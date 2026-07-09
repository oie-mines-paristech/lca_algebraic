import re
from collections import defaultdict
from copy import copy
from dataclasses import dataclass
from itertools import product
from types import FunctionType
from typing import Dict, List, Optional, Tuple, Union

import brightway2 as bw
import numpy as np
import pandas as pd
import sympy
from bw2data.backends.peewee import Activity
from pandas import DataFrame
from peewee import DoesNotExist
from pint import Quantity, Unit
from sympy import Add, Basic, Expr, IndexedBase, MatrixBase, Mul
from typing_extensions import deprecated

from lca_algebraic.activity import ActivityExtended, newActivity
from lca_algebraic.axis_dict import AxisDict
from lca_algebraic.base_utils import (
    MethodKey,
    TabbedDataframe,
    ValueOrExpression,
    _actName,
    _getDb,
    _isnumber,
    _isOutputExch,
)
from lca_algebraic.lambda_expression import LambdaExpr
from lca_algebraic.settings import Settings

from .cache import ExprCache, LCIACache
from .database import BIOSPHERE_PREFIX, DbContext, _isForeground, _setMeta
from .log import info, logger, warn
from .matrix import ActMatrix, invert
from .methods import method_name, method_unit
from .params import (
    FixedParamMode,
    _compute_param_length,
    _expanded_names_to_names,
    _fixed_params,
    _getAmountOrFormula,
    _param_registry,
    _toSymbolDict,
    all_params,
    freezeParams,
)
from .settings import PROXY_DB_FLAG, temp_settings
from .sympy_utils import replace_and_cleanup

# Symbol use in lambda expression to designate impact values
IMPACTS_SYMBOL = IndexedBase("impacts")


def register_user_function(sym, func):
    """Register a custom function with is python implementation
    Parameters
    ----------
    sym : the sympy function expression
    func : the implementation of the function

    Examples
    --------
    >>> def func_add(*args):
            returm sum(*args)
    >>>
    >>> func_add_sym = register_user_function(sympy.Function('func_add', real=True, imaginary=False), func_add)
    >>> e = sympy.Symbol('a') * func_add_sym(sympy.Symbol('b'), sympy.Symbol('c'))
    >>> sympy.srepr(e)
    "Mul(Symbol('a'), Function('func_add')(Symbol('b'), Symbol('c')))"
    """

    from lca_algebraic.base_utils import _user_functions

    global _user_functions
    _user_functions[sym.name] = (sym, func)
    return sym


def user_function(real=True, imaginary=False):
    """Function decorator to register custom Sympy user function.
    Beware that the implementation of the custom function should be compatible with numpy and work on vectors of values.

    Usage
    -----
    >>> @user_function()
    >>> def func_add(a, b):
    >>>      return a +b
    >>>
    >>> e = sympy.Symbol('a') * func_add(sympy.Symbol('b'), sympy.Symbol('c'))
    >>> sympy.srepr(e)
    "Mul(Symbol('a'), Function('func_add')(Symbol('b'), Symbol('c')))"
    """

    def lamb_f(func):
        sym_func = sympy.Function(func.__name__, real=real, imaginary=imaginary)
        return register_user_function(sym_func, func)

    return lamb_f


def _multiLCA(activities, methods):
    """Simple wrapper around brightway API"""
    bw.calculation_setups["process"] = {"inv": activities, "ia": methods}
    lca = bw.MultiLCA("process")
    cols = [_actName(act) for act_amount in activities for act, amount in act_amount.items()]
    return pd.DataFrame(lca.results.T, index=[method_name(method) for method in methods], columns=cols)


def multiLCA(models, methods, **params):
    """Compute LCA for a single activity and a set of methods,
    after settings the parameters and updating exchange amounts.
    This function does not use algebraic factorization and just calls the Brightway2 vanilla code.

    Parameters
    ----------
    model : Single activity (root model) or list of activities
    methods : Impact methods to consider
    params : Other parameters of the model
    """

    if not isinstance(models, list):
        models = [models]

    # Freeze params
    dbs = set(model[0] for model in models)
    for db in dbs:
        if _isForeground(db):
            freezeParams(db, **params)

    activities = [{act: 1} for act in models]
    return _multiLCA(activities, methods).transpose()


""" Compute LCA and return (act, method) => value """


def _group_acts_by_db(acts: list[ActivityExtended]):
    res = defaultdict(list)
    for act in acts:
        res[act["database"]].append(act)
    return res


def _multiLCAWithProxies(acts: list[ActivityExtended], methods: list[MethodKey]):
    proxy_acts = _createTechProxysForBio(acts)
    return _multiLCA([{act: 1} for act in proxy_acts.values()], methods)


def _multiLCAWithCache(all_acts, methods) -> Dict[Tuple[ActivityExtended, MethodKey], float]:
    res = dict()

    # Split activities by db_name
    for db_name, acts in _group_acts_by_db(all_acts).items():
        with LCIACache(db_name) as cache:
            # List activities with at least one missing value
            remaining_acts = list(act for act in all_acts if any(method for method in methods if (act, method) not in cache.data))

            # list methods with at least one missing value
            remaining_methods = list(
                method for method in methods if any(act for act in all_acts if (act, method) not in cache.data)
            )

            if len(remaining_acts) > 0 and len(remaining_methods) > 0:
                info(f"Computing LCA for {len(remaining_acts)} background acts on methods {remaining_methods}")

                lca = _multiLCAWithProxies(remaining_acts, remaining_methods)

                # Set output from dataframe
                for imethod, method in enumerate(remaining_methods):
                    for iact, act in enumerate(remaining_acts):
                        cache.data[(act, method)] = lca.iloc[imethod, iact]

            # Update res with a copy of the cache for selected impacts and activities
            res.update({(act, method): cache.data[(act, method)] for act in acts for method in methods})

    return res


def _actToLambdaExpr(model: Activity, axis=None, alpha=1) -> LambdaExpr:
    """
    For a given activity return a dict of bg _activity => expression
    """

    db_name = model["database"]

    if not _isForeground(db_name):
        # Already Bg activity ?
        res = LambdaExpr(alpha * IMPACTS_SYMBOL[0])
        res.background_activities = [model]
        return res

    # Key
    base_key = (axis, "factorized") if Settings.factorize_static_bg else (axis)

    # Cache around matrix computation
    with ExprCache(db_name) as cache:
        matrix_key = (base_key, "matrices")

        if matrix_key not in cache.data:
            logger.debug(f"{db_name} matrices were not in expression cache, computing...")
            matrices = _computeMatrices(db_name=db_name, axis=axis)
            cache.data[matrix_key] = matrices

        matrices: MatricesResult = cache.data[matrix_key]

        # Compute expression for a single fg model
        model_key = (base_key, "model", str(model))
        if model_key not in cache.data:
            cache.data[model_key] = matrices.compute_expr_for_model(model=model, with_axis=axis is not None, alpha=alpha)

        return cache.data[model_key]


def lambdify_expr(expr):
    return LambdaExpr(expr, params=[param.name for param in _param_registry().values()])


def _preMultiLCAAlgebric(
    model: ActivityExtended, methods: MethodKey, alpha: ValueOrExpression = 1, axis=None
) -> list[LambdaExpr]:
    """
    This method transforms an activity into a set of functions ready to compute LCA very fast on a set on methods.
    You may use is and pass the result to postMultiLCAAlgebric for fast computation on a model that does not change.

    This method is used by multiLCAAlgebric
    """

    with DbContext(model):
        if isinstance(alpha, Quantity):
            alpha = alpha.magnitude

        def _key(method):
            return (model, axis, method, alpha)

        with ExprCache(model["database"]) as cache:
            missing_methods = [method for method in methods if not _key(method) in cache.data]
            if len(missing_methods) > 0:
                exprs = _modelToExpr(model, methods=missing_methods, axis=axis, alpha=alpha)
                for method, expr in zip(missing_methods, exprs):
                    cache.data[_key(method)] = expr

            # At this point, everything is in cache
            # REturn the list in order
            return list(cache.data[_key(method)] for method in methods)


def _modelToExpr(model: Activity, methods: List[MethodKey], axis=None, alpha: ValueOrExpression = 1) -> List[LambdaExpr]:
    """
    Compute expressions corresponding to a model for each impact, replacing activities by the value of its impact

    Return
    ------
    <list of expressions (one per impact)>, <list of required param names>
    """

    generic_lambda_expr = _actToLambdaExpr(model, axis=axis, alpha=alpha)

    if len(generic_lambda_expr.background_activities) == 0:
        zero = LambdaExpr(AxisDict({"_all_": 0.0}))
        zero.background_activities = []
        zero.impacts = []
        return [zero] * len(methods)

    # Compute LCA for background activities
    all_impacts = _multiLCAWithCache(generic_lambda_expr.background_activities, methods)

    res = list()

    for method in methods:
        impacts = {bg_act: all_impacts[bg_act, method] for bg_act in generic_lambda_expr.background_activities}
        res.append(generic_lambda_expr.with_impacts(impacts))

    return res


def _slugify(str):
    return re.sub("[^0-9a-zA-Z]+", "_", str)


@dataclass
class ResultsWithParams:
    """Holds bith the result with context parameters"""

    dataframe: pd.DataFrame
    params: Dict


def _postMultiLCAAlgebric(methods, lambdas: List[LambdaExpr], with_params=False, unit: Unit = None, **params):
    """
    Compute LCA for a given set of parameters and pre-compiled lambda functions.
    This function is used by **multiLCAAlgebric**

    Parameters
    ----------
    methodAndLambdas : Output of preMultiLCAAlgebric
    **params : Parameters of the model
    """

    param_length = _compute_param_length(params)

    # lambda are SymDict ?
    # If use them as number of params
    if lambdas[0].has_axis and len(lambdas[0].expr) > 1:
        if param_length > 1:
            raise Exception("Multi params cannot be used together with 'axis'")
        param_length = len(lambdas[0].axis_keys)

    # Init output
    res = np.zeros((len(methods), param_length), float)

    # All params
    context_params = dict()

    # Compute result on whole vectors of parameter samples at a time : lambdas use numpy for vector computation
    def process(lambd):
        nonlocal context_params

        value_context = lambd.compute(**params)

        # Update the param values used
        context_params.update(value_context.context)

        value = value_context.value
        # Expand axis values as a list, to fit into the result numpy array
        if isinstance(value, dict):
            if len(value) > 1:
                # Ensure the values are in the same order as the value

                # XXX We use the order of the first lambda as each one might have different order
                axes = lambdas[0].axis_keys
                xvalue = [0.0] * len(axes)
                for i, axis_tag in enumerate(axes):
                    if axis_tag not in value:
                        continue
                    if axis_tag == "_all_":
                        xvalue[i] = float(value[axis_tag])
                    else:
                        xvalue[i] = float(value["_all_"]) - float(value[axis_tag])
                value = xvalue
            else:
                value = value["_all_"]
        return value

    # Use multithread for that
    for imethod, lambd in enumerate(lambdas):
        value = process(lambd)
        res[imethod, :] = value

    result = pd.DataFrame(
        res,
        index=[method_name(method) + "[%s]" % method_unit(method, fu_unit=unit) for method in methods],
    ).transpose()

    if with_params:
        return ResultsWithParams(dataframe=result, params=context_params)
    else:
        return result


# Add default values for issing parameters or warn about extra params
def _filter_params(params, expected_names, model):
    res = params.copy()

    expected_params_names = _expanded_names_to_names(expected_names)
    for expected_name in expected_params_names:
        if expected_name not in params:
            default = _param_registry()[expected_name].default
            res[expected_name] = default
            warn("Missing parameter %s, replaced by default value %s" % (expected_name, default))

    for key, value in params.items():
        if key not in expected_params_names:
            del res[key]
            if model:
                warn("Param %s not required for model %s" % (key, model))
    return res


def compute_value(formula, **params):
    """Compute actual value for a given formula, with possible parameters (or default ones)"""
    if isinstance(formula, float) or isinstance(formula, int):
        return formula

    lambd = LambdaExpr(formula)

    value_context = lambd.compute(**params)

    # TODO: add an option to keep axes ?
    value = value_context.value
    if isinstance(value, dict):
        return value["_all_"]
    return value


@deprecated("multiLCAAlgebric is deprecated, use compute_impacts instead")
def multiLCAAlgebric(*args, **kwargs):
    """deprecated. `compute_impacts()` instead"""
    warn("multiLCAAlgebric is deprecated, use compute_impacts instead")
    return compute_impacts(*args, **kwargs)


def _params_dataframe(param_values: Dict[str, float]):
    """Create a DataFrame, ordered by group, showing param values"""
    params_by_name = all_params()

    records = []

    plen = _compute_param_length(param_values)

    for param_name, value in param_values.items():
        param = params_by_name[param_name]
        record = {
            "group": param.group if param.group is not None else "",
            "name": param.name,
            "min": param.min,
            "max": param.max,
            "default": param.default,
        }

        if plen == 1:
            record["value"] = value
        else:
            if isinstance(value, (list, np.ndarray)):
                record.update({f"value_{i}": value for i, value in enumerate(value, 1)})
            else:
                # Repeat single value
                record.update({f"value_{i}": value for i in range(1, plen + 1)})

        records.append(record)

    df = pd.DataFrame.from_records(records).set_index(["group", "name"]).sort_index()

    return df


SingleOrMultipleFloat = Union[float, List[float], np.ndarray]


def compute_inventory(
    model: ActivityExtended,
    functional_unit=1,
    as_dict=False,
    impact_method=None,
    fields=["database", "name", "location", "unit"],
    **params,
):
    """

    This method computes the inventory of background activities for a given scenario / values of parameters.

    Parameters
    ----------
    model:
        Root activity

    functional_unit:
        Quantitity to divide the inventory. 1 by default

    as_dict:
        If true, returns a dict of act => value. If false (default) returns a dataframe

    fields:
        List of fields to be added in the ouput dataframe

    impact_method:
        If provided, return the impact for each activity, rather that its quantity

    params:
        All other attributes are treated as values of lca_algebraic parameters.
        If not specified, each parameters takes its default value.

    Returns
    -------
    Dataframe or Dict of act => value

    """

    with temp_settings(factorize_static_bg=False):
        lambda_expr = _actToLambdaExpr(model, alpha=1 / functional_unit)

        # Transform to dict of act => value
        val_by_act = dict()
        for bg_act in lambda_expr.background_activities:
            # Dummy impact set to 1 only for current bg act, zero for others
            impacts = {act: 1.0 if act == bg_act else 0.0 for act in lambda_expr.background_activities}

            bg_expr = lambda_expr.with_impacts(impacts)
            val = bg_expr.compute(**params).value
            if isinstance(val, dict):
                val = val["_all_"]
            val_by_act[bg_act] = val

        if impact_method is not None:
            # Compute LCA of background activities
            impact_by_act = _multiLCAWithCache(val_by_act.keys(), [impact_method])

            val_by_act: {act: value * impact_by_act[(act, impact_method)] for act, value in val_by_act.items()}

        if as_dict:
            return val_by_act

        # Transform to dataframe
        items = []
        for act, value in val_by_act.items():
            item = dict()

            for field in fields:
                if field in act:
                    item[field] = act[field]

            item["value"] = value
            items.append(item)

    return DataFrame(items)


def compute_impacts(
    models,
    methods,
    axis=None,
    functional_unit=1,
    return_params=False,
    description=None,
    **params,
):
    """
    Main parametric LCIA method :
    Computes LCA by expressing the foreground model as symbolic expression of background activities and parameters.
    Then, compute 'static' inventory of the referenced background activities.
    This enables a very fast recomputation of LCA with different parameters, \
    useful for stochastic evaluation of parametrized model

    Parameters
    ----------
    models :
        Single model or
        List of model or
        List of (model, alpha)
        or Dict of model:amount
        In case of several models, you cannot use list of parameters

    methods :
        List of methods / impacts to consider

    axis:
        Designates the name of a custom attribute of foreground activities.
        You may set this attribute using the method `myActivity.updateMeta(your_custom_attr="some_value")`

        The impacts will be ventilated by this attribute.
        This is useful to get impact by phase or sub-modules.

    params:
        Any other argument passed to this function is considered as a value of a parameter of the model :
        Values can be either single float values, list or ndarray of values.
        In the later case, all parameters should have the same number of values.
        Paremeters that are not provided will have their default value set.

    functional_unit:
        Quantity (static or Sympy formula) by which to divide impacts. Optional, 1 by default.

    return_params:
        If true, also returns the value of all parameters in as tabbed DataFrame

    description:
        Optional description/metadata to be added in output when using "return params" Dataframe

    Returns
    -------
    A dataframe with the results. If *return_params* is true, it returns `TabbedDataframe`,
    including all parameters values, that can be saved as a multi sheet excel file.

    Examples
    --------
    >>> compute_impacts(
    >>>    mainAct1, # The root activity of the foreground model
    >>>    [climate_change], # climate_change is the key (tuple) of the impact method
    >>>    functional_unit=energy_expression, # energy expression is a Sympy expression computing the energy in kWh
    >>>    axis="phase", # Split results by phase
    >>>    return_params=True, # Return all parameter values
    >>>
    >>>    # Parameter values
    >>>    p1=2.0,
    >>>    p2=3.0)


    """
    dfs = dict()

    if isinstance(models, list):

        def to_tuple(item):
            if isinstance(item, tuple):
                return item
            else:
                return (item, 1)

        models = dict(to_tuple(item) for item in models)
    elif not isinstance(models, dict):
        models = {models: 1}

    # Gather all param values (even default and computed)
    params_all = dict()

    # Single method provided ?
    if isinstance(methods, tuple):
        methods = [methods]

    for model, alpha in models.items():
        if type(model) is tuple:
            model, alpha = model

        alpha = float(alpha)

        dbname = model.key[0]
        with DbContext(dbname):
            # Check no params are passed for FixedParams
            for key in params:
                if key in _fixed_params():
                    warn("Param '%s' is marked as FIXED, but passed in parameters : ignored" % key)

            if functional_unit != 1:
                alpha = alpha / functional_unit

            lambdas = _preMultiLCAAlgebric(model, methods, alpha=alpha, axis=axis)

            unit: Optional[Unit] = functional_unit.units if isinstance(functional_unit, Quantity) else None

            res = _postMultiLCAAlgebric(methods, lambdas, with_params=return_params, unit=unit, **params)

            if return_params:
                df = res.dataframe
                params_all.update(res.params)
            else:
                df = res

            model_name = _actName(model)
            while model_name in dfs:
                model_name += "'"

            # param with several values
            list_params = {k: vals for k, vals in params.items() if isinstance(vals, list)}

            # Shapes the output / index according to the axis or multi param entry
            if axis:
                df[axis] = lambdas[0].axis_keys
                df = df.set_index(axis)
                df.index.set_names([axis])

                # Filter out line with zero output
                df = df.loc[
                    df.apply(
                        lambda row: not (row.name is None and row.values[0] == 0.0),
                        axis=1,
                    )
                ]

                # Rename "None" to others
                df = df.rename(index={None: "_other_"})

                # Sort index
                df.sort_index(inplace=True)

                # compute impacts that do not belong axis
                df.loc["*other*"] = df.loc["_all_"] - df.loc[list(set(df.index.to_list()) - set(["_all_"]))].sum()

                # Add "total" line
                df = df.rename(index={"_all_": "*all*"})

            elif len(list_params) > 0:
                # Use param values as index
                df.index = pd.MultiIndex.from_frame(pd.DataFrame(list_params))
                if df.index.nlevels == 1:
                    df.index = df.index.get_level_values(0)

            else:
                # Single output ? => give the single row the name of the model activity
                df = df.rename(index={0: model_name})

            dfs[model_name] = df

    if len(dfs) == 1:
        df = list(dfs.values())[0]
    else:
        # Concat several dataframes for several models
        df = pd.concat(list(dfs.values()))

    if return_params:
        metadata = {"Models": str(models), "Functional unit": functional_unit}
        if description:
            metadata["Description"] = description

        return TabbedDataframe(metadata=metadata, Results=df, Parameters=_params_dataframe(params_all))
    else:
        return df


def _isBioAct(act: ActivityExtended):
    db_name = act["database"]
    return (BIOSPHERE_PREFIX in db_name) or ("type" in act and act["type"] in ["emission", "natural resource"])


def _getOrCreateProxyDb(db_name):
    """Init proxy db to biosphere if not done yet"""
    proxyname = db_name + "-proxy"
    if proxyname not in bw.databases:
        db = bw.Database(proxyname)
        db.write(dict())
        _setMeta(proxyname, PROXY_DB_FLAG, True)
    return proxyname


def _getOrCreateProxy(act: ActivityExtended, exchanges: dict[ActivityExtended, float]):
    proxy_db = _getOrCreateProxyDb(act["database"])
    proxy_code = act["code"] + "#proxy"
    with temp_settings(strict_mode=False):
        try:
            proxy = _getDb(proxy_db).get(proxy_code)

            # Check exchanges
            existing_exchanges = {ex[1]: ex[2] for ex in proxy.listExchanges()}

            if existing_exchanges != exchanges:
                info(f"Proxy exchanges differ in {proxy['name']}, updating them")

                if len(existing_exchanges) > 0:
                    proxy.deleteExchanges()
                proxy.addExchanges(exchanges)

        except DoesNotExist:
            name = act["name"] + " # proxy"

            # Create biosphere proxy in User Db
            proxy = newActivity(
                db_name=proxy_db,
                name=name,
                code=proxy_code,
                switchActivity=True,
                isProxy=True,
                unit=act["unit"],
                exchanges=exchanges,
            )

    return proxy


def _createTechProxysForBio(acts: List[ActivityExtended]) -> Dict[ActivityExtended, ActivityExtended]:
    """
    Potentially create tech proxys for bio activity (brightway cannot proces LCIA on them
    Returns a dict of [OriginalAct -> OriginalOrProxyAct]
    """
    res = dict()
    for act in acts:
        res[act] = act if not _isBioAct(act) else _getOrCreateProxy(act, {act: 1})
    return res


def _replace_fixed_params(expr: Expr, fixed_params, fixed_mode=FixedParamMode.DEFAULT):
    """Replace fixed params with their value."""
    if not isinstance(expr, Basic):
        return expr

    sub = {key: val for param in fixed_params for key, val in param.expandParams(param.stat_value(fixed_mode)).items()}
    if len(sub) == 0:
        return expr
    sub = _toSymbolDict(sub)
    return replace_and_cleanup(expr, sub)


def _get_axis(act, axis_name: str):
    """Safe"""
    tag = act.get(axis_name, None)

    if tag is None:
        return None
    if tag.isalnum():
        return tag
    else:
        return re.sub("[^0-9a-zA-Z]+", "_", tag)


def _force_reduce(expr):
    """Force reduction of sum and multiplication : usefull for AxisDict"""
    if isinstance(expr, AxisDict):
        return AxisDict({key: _force_reduce(val) for key, val in expr.items()})
    if isinstance(expr, dict):
        return _force_reduce(AxisDict(expr))
    if isinstance(expr, Add):
        res = 0.0
        for arg in expr.args:
            res += _force_reduce(arg)
        return res
    if isinstance(expr, Mul):
        res = 1.0
        for arg in expr.args:
            res *= _force_reduce(arg)
        return res
    return expr


def replace_fixed_params(expr: Expr):
    return _replace_fixed_params(expr, _fixed_params().values())


def _clean_expr(expr):
    expr = _force_reduce(expr)
    return _replace_fixed_params(expr, _fixed_params().values())


def _solve_expression(fg_matrix: ActMatrix, bg_matrix: ActMatrix) -> ActMatrix:
    """Solve the foreground inventory."""

    # Case of empty matrix
    if len(bg_matrix.cols_acts()) == 0 or len(fg_matrix.row_acts()) == 0:
        ret = ActMatrix()
        ret._col_acts = copy(bg_matrix._col_acts)
        ret._row_acts = copy(fg_matrix._row_acts)
        return ret

    A = fg_matrix.to_sympy()
    B = bg_matrix.to_sympy()

    # Use fast inversion of sparse matrices
    inv_A = invert(A)
    res_mat = inv_A @ B

    ret = ActMatrix()
    ret._col_acts = copy(bg_matrix._col_acts)
    ret._row_acts = copy(fg_matrix._row_acts)
    for (i0, k0), (i1, k1) in product(enumerate(ret.row_acts()), enumerate(ret.cols_acts())):
        if (val := res_mat[i0, i1]) == 0:
            continue
        super(defaultdict, ret).__setitem__((k0, k1), _clean_expr(val))

    return ret


@dataclass
class MatricesResult:
    inv_fg: MatrixBase
    fg_acts: List[Activity]
    bg_acts: List[Activity]
    bg_matrix: Optional[MatrixBase] = None
    bg_cells: Optional[Dict[Tuple[Activity, Activity], ValueOrExpression]] = None

    def compute_expr_for_model(self, model: Activity, with_axis: bool, alpha=1) -> LambdaExpr:
        """Compute dict of bg_act => model for foreground model"""

        model_idx = self.fg_acts.index(model)

        axis_terms = defaultdict(list)

        if self.bg_cells is not None:
            for j, bg_act in enumerate(self.bg_acts):
                val = self.bg_cells.get((model, bg_act), 0)
                if val == 0:
                    continue
                if isinstance(val, AxisDict):
                    for k, v in val.items():
                        if v != 0:
                            axis_terms[k].append(Mul(v, IMPACTS_SYMBOL[j]))
                else:
                    axis_terms["_all_"].append(Mul(val, IMPACTS_SYMBOL[j]))
        else:
            fg_col = self.inv_fg.row(model_idx)
            row = fg_col * self.bg_matrix
            for j, bg_act in enumerate(self.bg_acts):
                val = row[j]
                if val == 0:
                    continue
                axis_terms["_all_"].append(Mul(val, IMPACTS_SYMBOL[j]))

        if with_axis:
            expr = AxisDict({k: Mul(Add(*terms), alpha) for k, terms in axis_terms.items() if terms})
            expr = _force_reduce(expr)
        else:
            terms = axis_terms.get("_all_", [])
            inner = Mul(Add(*terms), alpha) if terms else alpha * 0.0
            expr = AxisDict({"_all_": inner})

        res = LambdaExpr(expr)
        res.background_activities = self.bg_acts

        return res


def _computeMatrices(db_name, axis=None) -> MatricesResult:
    """
    Compute all expressions for a given DB
    Returns Dict[FgAct => Dict[BGAct => Expressions]]
    """

    if not _isForeground(db_name):
        raise ValueError(f"Can only compute expression on foreground activities. {db_name} is background")

    # Square technospere matrix dict of <act1, act2> => float
    fg_matrix = ActMatrix()

    # Rectangle matrix
    bg_matrix = ActMatrix()

    # Gather [axis_name] -> activity list
    axis_acts = defaultdict(list)

    # Fill the fg_matrix and bg_matrix, collect axis data
    # This function walk within the activity tree, it's needed because
    # activities may be spread across severals databases other than db_name
    def walk_activities(act: Activity):
        nonlocal axis, fg_matrix, bg_matrix, axis_acts

        # Update axis values
        if axis is not None:
            new_axis_val = _get_axis(act, axis)
            if new_axis_val is not None:
                axis_acts[new_axis_val].append(act)

        # Should not happen ?
        if not _isForeground(act["database"]):
            # We reached a background DB ? => stop developping and create reference to activity
            return

        if act in fg_matrix.row_acts():
            # Already explored
            return

        # Add current act to axes matrices, to keep correct shape
        fg_matrix.add_row(act)
        fg_matrix.add_col(act)
        bg_matrix.add_row(act)

        static_bg_amounts = dict()

        for exch in act.exchanges():
            amount = _getAmountOrFormula(exch)

            if isinstance(amount, FunctionType):
                # Some amounts in EIDB are functions ... we ignore them
                continue

            # Fetch activity referenced by the exchange
            input_db, input_code = exch["input"]
            sub_act = _getDb(input_db).get(input_code)

            # Fill the appropriate matrix
            if _isForeground(input_db):
                #  Production exchange
                if not _isOutputExch(exch):
                    amount = -amount

                fg_matrix[act, sub_act] += amount

                # Recursively explore the rest
                walk_activities(sub_act)

            else:
                if Settings.factorize_static_bg and _isnumber(amount):
                    static_bg_amounts[sub_act] = amount
                else:
                    bg_matrix[act, sub_act] += amount

        # Static bg amount not empty ? create and reference it
        if len(static_bg_amounts) > 0:
            proxy_act = _getOrCreateProxy(act, static_bg_amounts)
            bg_matrix[act, proxy_act] = 1.0

    # Fill the matrices exploring everything
    for act in bw.Database(db_name):
        walk_activities(act)

    A = fg_matrix.to_sympy()
    A = replace_fixed_params(A)

    inv_A = invert(A)

    if axis is None:
        B = bg_matrix.to_sympy()
        B = replace_fixed_params(B)
        return MatricesResult(
            bg_matrix=B,
            inv_fg=inv_A,
            fg_acts=fg_matrix.cols_acts(),
            bg_acts=bg_matrix.cols_acts(),
        )

    # solve the linear equation algebically for the axis _all_
    data = _solve_expression(fg_matrix, bg_matrix)

    tmp = defaultdict(AxisDict)
    for k, v in data.items():
        tmp[k] += AxisDict({"_all_": v})

    # Solve LCA equation for remaining axis
    for axis_tag, acts in axis_acts.items():
        xfg = copy(fg_matrix)
        xbg = copy(bg_matrix)

        # Set all activities belong axis activities to 0
        for a in acts:
            for o in xfg.cols_acts():
                xfg[a, o] = 0.0
            for o in xbg.cols_acts():
                xbg[a, o] = 0.0
            xfg[a, a] = 1.0

        out = _solve_expression(xfg, xbg)
        for k, v in out.items():
            tmp[k] += AxisDict({axis_tag: v})

    return MatricesResult(
        bg_cells=dict(tmp),
        inv_fg=inv_A,
        fg_acts=fg_matrix.cols_acts(),
        bg_acts=bg_matrix.cols_acts(),
    )


def _reverse_dict(dic):
    return {v: k for k, v in dic.items()}
