import concurrent.futures
import re
from dataclasses import dataclass
from types import FunctionType
from typing import Dict, List, Optional, Tuple, Union

import brightway2 as bw
import numpy as np
import pandas as pd
from black.trans import defaultdict
from bw2data.backends.peewee import Activity
from pandas import DataFrame
from peewee import DoesNotExist
from pint import Quantity, Unit
from sympy import Add, Basic, Expr, ImmutableMatrix, Mul, Symbol, lambdify, parse_expr
from sympy.printing.numpy import NumPyPrinter
from typing_extensions import deprecated

from . import newActivity
from .activity import ActivityExtended
from .axis_dict import AxisDict
from .base_utils import (
    MethodKey,
    TabbedDataframe,
    ValueOrExpression,
    _actName,
    _getDb,
    _isOutputExch,
    _user_functions,
)
from .cache import ExprCache, LCIACache
from .database import BIOSPHERE_PREFIX, DbContext, _isForeground
from .log import logger, warn
from .methods import method_name, method_unit
from .params import (
    FixedParamMode,
    _complete_params,
    _compute_param_length,
    _expand_param_names,
    _expand_params,
    _expanded_names_to_names,
    _fixed_params,
    _getAmountOrFormula,
    _param_registry,
    _toSymbolDict,
    all_params,
    freezeParams,
)


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
    global _user_functions
    _user_functions[sym.name] = (sym, func)
    return sym


def user_function(sym):
    """Function decorator to register user function

    Usage
    -----
    >>> @user_function(sympy.Function('func_add', real=True, imaginary=False))
    >>> def func_add(*args):
            returm sum(*args)
    >>>
    >>> e = sympy.Symbol('a') * func_add(sympy.Symbol('b'), sympy.Symbol('c'))
    >>> sympy.srepr(e)
    "Mul(Symbol('a'), Function('func_add')(Symbol('b'), Symbol('c')))"
    """
    return lambda func: register_user_function(sym, func)


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


def _multiLCAWithCache(acts, methods) -> Dict[Tuple[ActivityExtended, MethodKey], float]:
    # Create proxys for bio activities
    proxy_acts = _createTechProxysForBio(acts)

    with LCIACache() as cache:
        # List activities with at least one missing value
        remaining_acts = list(
            proxy_act
            for proxy_act in proxy_acts.values()
            if any(method for method in methods if (proxy_act, method) not in cache.data)
        )

        if len(remaining_acts) > 0:
            lca = _multiLCA([{act: 1} for act in remaining_acts], methods)

            # Set output from dataframe
            for imethod, method in enumerate(methods):
                for iact, act in enumerate(remaining_acts):
                    cache.data[(act, method)] = lca.iloc[imethod, iact]

        # Return a copy of the cache for selected impacts and activities
        return {(act, method): cache.data[(proxy_acts[act], method)] for act in acts for method in methods}


def _replace_symbols_with_params_in_exp(expr: Basic):
    """After unpickling, the Symbols in the expression with same named are not the same object than Parameters.
    This is problematic for xreplace wich relies on it :
    Here we replace them.
    """

    if not isinstance(expr, Basic):
        return expr

    all_params = _param_registry().as_dict()
    subs = {symb: all_params[symb.name] for symb in expr.free_symbols if isinstance(symb, Symbol) and symb.name in all_params}

    logger.debug("Replace: %s", subs)

    expr = expr.xreplace(subs)
    return expr


def _actToExpressionDict(model: Activity, axis=None) -> Dict[Activity, ValueOrExpression]:
    """
    For a given activity return a dict of bg _activity => expression
    """

    db_name = model["database"]

    # Try to load from cache
    with ExprCache() as cache:
        key = (db_name, axis)

        if key not in cache.data:
            logger.debug(f"{db_name} was not in expression cache, computing...")

            res_by_act = _computeDbExpressions(db_name=db_name, axis=axis)

            cache.data[key] = res_by_act
        else:
            logger.debug(f"{model} found in expression cache")

        res_by_act = cache.data[key]

    res = res_by_act[model]

    return res


def _modelToExpr(model: Activity, methods: List[MethodKey], axis=None, alpha: ValueOrExpression = 1):
    """
    Compute expressions corresponding to a model for each impact, replacing activities by the value of its impact

    Return
    ------
    <list of expressions (one per impact)>, <list of required param names>

    """

    expr_by_bg_act = _actToExpressionDict(model, axis=axis)

    bg_acts = list(expr_by_bg_act.keys())

    # Compute LCA for background activities
    impacts = _multiLCAWithCache(bg_acts, methods)

    # For each method, compute an algebric expression with activities replaced by their values
    exprs = []
    for method in methods:
        expr = 0.0
        for bg_act, value in expr_by_bg_act.items():
            impact = impacts[(bg_act, method)]
            expr += impact * value

        # Ensure symbols are params
        expr = _replace_symbols_with_params_in_exp(expr)

        if isinstance(alpha, Quantity):
            alpha = alpha.magnitude
        expr *= alpha

        exprs.append(expr)

    return exprs


def _filter_param_values(params, expanded_param_names):
    return {key: val for key, val in params.items() if key in expanded_param_names}


def _free_symbols(expr: Basic):
    if isinstance(expr, Basic):
        return set([str(symb) for symb in expr.free_symbols])
    else:
        # Static value
        return set()


def _lambdify(expr: Basic, expanded_params):
    """Lambdify, handling manually the case of SymDict (for impacts by axis)"""

    printer = NumPyPrinter(
        {"fully_qualified_modules": False, "inline": True, "allow_unknown_functions": True, "user_functions": dict()}
    )

    modules = [{x[0].name: x[1] for x in _user_functions.values()}, "numpy"]

    if isinstance(expr, Basic):
        lambd = lambdify(expanded_params, expr, modules, printer=printer, cse=LambdaWithParamNames._use_sympy_cse)

        def func(*arg, **kwargs):
            res = lambd(*arg, **kwargs)
            if isinstance(res, dict):
                # Transform key symbols into Str
                return {str(k): v for k, v in res.items()}
            else:
                return res

        return func

    else:
        # Not an expression : return static func
        def static_func(*args, **kargs):
            return expr

        return static_func


@dataclass
class ValueContext:
    """Represents a result value, with all parameters values used in context"""

    value: float
    context: Dict[str, float]


class LambdaWithParamNames:
    """
    This class represents a compiled (lambdified) expression together
    with the list of requirement parameters and the source expression
    """

    _use_sympy_cse = False

    def __init__(self, expr: Expr, expanded_params=None, params=None, sobols=None):
        """Computes a lamdda function from expression and list of expected parameters.
        you can provide either the list pf expanded parameters (full vars for enums) for the 'user' param names
        """

        if isinstance(expr, dict):
            # Come from JSON serialization
            obj = expr
            # LIst of required params for this lambda
            self.params: List[str] = obj["params"]

            # Full names
            self.expanded_params = _expand_param_names(self.params)
            local_dict = {x[0].name: x[0] for x in _user_functions.values()}
            self.expr = parse_expr(obj["expr"], local_dict=local_dict)
            self.lambd = _lambdify(self.expr, self.expanded_params)
            self.sobols = obj["sobols"]

        else:
            self.expr = expr
            self.params = params

            if expanded_params is None:
                if params is None:
                    expanded_params = _free_symbols(expr)
                    params = _expanded_names_to_names(expanded_params)
                    self.params = params

                # We expand again the parameters
                # If we expect an enum param name, we also expect the other ones :
                # enumparam_val1 => enumparam_val1, enumparam_val2, ...
                expanded_params = _expand_param_names(params)

            elif self.params is None:
                self.params = _expanded_names_to_names(expanded_params)

            self.lambd = _lambdify(expr, expanded_params)
            self.expanded_params = expanded_params
            self.sobols = sobols

    @property
    def has_axis(self):
        return isinstance(self.expr, AxisDict)

    @property
    def axis_keys(self):
        if self.has_axis:
            return self.expr.str_keys()
        else:
            return None

    def compute(self, **params) -> ValueContext:
        """Compute result value based of input parameters"""

        # Add default or computed values
        completed_params = _complete_params(params, self.params)

        # Expand enums
        expanded_params = _expand_params(completed_params)

        # Remove parameters that are not required
        expanded_params = _filter_param_values(expanded_params, self.expanded_params)

        value = self.lambd(**expanded_params)

        return ValueContext(value=value, context=completed_params)

    def serialize(self):
        expr = str(self.expr)
        return dict(params=self.params, expr=expr, sobols=self.sobols)

    @staticmethod
    def use_sympy_cse(b=True):
        LambdaWithParamNames._use_sympy_cse = b

    def __repr__(self):
        return repr(self.expr)

    def _repr_latex_(self):
        return self.expr._repr_latex_()


def lambdify_expr(expr):
    return LambdaWithParamNames(expr, params=[param.name for param in _param_registry().values()])


def _preMultiLCAAlgebric(model: ActivityExtended, methods: MethodKey, alpha: ValueOrExpression = 1, axis=None):
    """
    This method transforms an activity into a set of functions ready to compute LCA very fast on a set on methods.
    You may use is and pass the result to postMultiLCAAlgebric for fast computation on a model that does not change.

    This method is used by multiLCAAlgebric
    """

    with DbContext(model):
        exprs = _modelToExpr(model, methods, alpha=alpha, axis=axis)

        # Lambdify (compile) expressions
        return [LambdaWithParamNames(expr) for expr in exprs]


def _slugify(str):
    return re.sub("[^0-9a-zA-Z]+", "_", str)


@dataclass
class ResultsWithParams:
    """Holds bith the result with context parameters"""

    dataframe: pd.DataFrame
    params: Dict


def _postMultiLCAAlgebric(methods, lambdas: List[LambdaWithParamNames], with_params=False, unit: Unit = None, **params):
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
    if lambdas[0].has_axis:
        if param_length > 1:
            raise Exception("Multi params cannot be used together with 'axis'")
        param_length = len(lambdas[0].axis_keys)

    # Init output
    res = np.zeros((len(methods), param_length), float)

    # All params
    context_params = dict()

    # Compute result on whole vectors of parameter samples at a time : lambdas use numpy for vector computation
    def process(args):
        nonlocal context_params

        imethod = args[0]
        lambd: LambdaWithParamNames = args[1]

        value_context = lambd.compute(**params)

        # Update the param values used
        context_params.update(value_context.context)

        value = value_context.value

        # Expand axis values as a list, to fit into the result numpy array
        if isinstance(value, dict):
            # Ensure the values are in the same order as the value

            # XXX We use the order of the first lambda as each one might have different order
            axes = lambdas[0].axis_keys
            value = list(float(value[axis]) if axis in value else 0.0 for axis in axes)

        return (imethod, value)

    # Use multithread for that
    with concurrent.futures.ThreadPoolExecutor() as exec:
        for imethod, value in exec.map(process, enumerate(lambdas)):
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

    lambd = LambdaWithParamNames(formula)

    value_context = lambd.compute(**params)

    return value_context.value


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
    model: ActivityExtended, functional_unit=1, as_dict=False, fields=["database", "name", "location", "unit"], **params
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

    params:
        All other attributes are treated as values of lca_algebraic parameters.
        If not specified, each parameters takes its default value.

    Returns
    -------
    Dataframe or Dict of act => value

    """

    exprs_by_bg_act = _actToExpressionDict(model)

    # Transform to dict of act => value
    val_by_act = dict()
    for bg_act, expr in exprs_by_bg_act.items():
        val_by_act[bg_act] = compute_value(expr, **params) / functional_unit

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

                # Add "total" line
                df.loc["*sum*"] = df.sum(numeric_only=True)

            elif len(list_params) > 0:
                for k, vals in list_params.items():
                    df[k] = vals
                df = df.set_index(list(list_params.keys()))

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


def _createTechProxysForBio(acts: List[ActivityExtended]) -> Dict[ActivityExtended, ActivityExtended]:
    """
    Potentially create tech proxys for bio activity (brightway cannot proces LCIA on them
    Returns a dict of [OriginalAct -> OriginalOrProxyAct]
    """
    res = dict()
    for act in acts:
        if not _isBioAct(act):
            proxy_act = act
        else:
            proxy_code = act["code"] + "#asTech"
            try:
                proxy_act = _getDb(act["database"]).get(proxy_code)
            except DoesNotExist:
                name = act["name"] + " # asTech"

                # Create biosphere proxy in User Db
                proxy_act = newActivity(
                    db_name=act["database"],
                    name=name,
                    code=proxy_code,
                    switchActivity=True,
                    isProxy=True,
                    unit=act["unit"],
                    exchanges={act: 1},
                )
        res[act] = proxy_act

    return res


def _replace_fixed_params(expr, fixed_params, fixed_mode=FixedParamMode.DEFAULT):
    """Replace fixed params with their value."""
    if not isinstance(expr, Basic):
        return expr

    sub = {key: val for param in fixed_params for key, val in param.expandParams(param.stat_value(fixed_mode)).items()}
    sub = _toSymbolDict(sub)
    return expr.xreplace(sub)


def _get_axis(act, axis_name: str):
    """Safe"""
    tag = act.get(axis_name, None)

    if tag is None:
        return None
    if tag.isalnum():
        return tag
    else:
        return re.sub("[^0-9a-zA-Z]+", "_", tag)


class ActMatrix(defaultdict):
    def __init__(self):
        super().__init__(lambda: 0.0)
        self._col_acts = list()
        self._row_acts = list()

    def __getitem__(self, key):
        row_act, col_act = key
        self.add_row(row_act)
        self.add_col(col_act)
        return super().__getitem__(key)

    def __setitem__(self, key, value):
        row_act, col_act = key
        self.add_row(row_act)
        self.add_col(col_act)
        return super().__setitem__(key, value)

    def add_col(self, col_act):
        if col_act not in self._col_acts:
            self._col_acts.append(col_act)
            self._col_acts.sort()

    def add_row(self, row_act):
        if row_act not in self._row_acts:
            self._row_acts.append(row_act)
            self._row_acts.sort()

    def demand_vector(self, act, value=1.0):
        """Generate vector of len (rows) with zeros and 1 only at the index of act"""
        act_idx = self.row_acts().index(act)
        res = [0.0] * len(self._row_acts)
        res[act_idx] = value
        return ImmutableMatrix([res])

    def row_acts(self):
        return self._row_acts

    def cols_acts(self):
        return self._col_acts

    def shape(self):
        return len(self._row_acts), len(self._col_acts)

    def to_sympy(self):
        """Return an immutable sympy matrix"""
        rows = list()
        for row_act in self.row_acts():
            row = list()
            rows.append(row)
            for col_act in self.cols_acts():
                row.append(self.get((row_act, col_act), 0.0))
        return ImmutableMatrix(rows)

    def to_dataframe(self):
        res = dict()
        for row_act in self.row_acts():
            res[str(row_act)] = {str(col_act): self[(row_act, col_act)] for col_act in self.cols_acts()}
        return pd.DataFrame(res)

    def __repr__(self):
        return self.to_dataframe().__repr__()


def _force_reduce(expr):
    """Force reduction of sum and multiplication : usefull for AxisDict"""
    if isinstance(expr, dict):
        return AxisDict(expr)
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


def _clean_expr(expr):
    expr = _force_reduce(expr)
    return _replace_fixed_params(expr, _fixed_params().values())


def _solve_expression(
    fg_matrix: ActMatrix, bg_matrix: ActMatrix
) -> Dict[ActivityExtended, Dict[ActivityExtended, ValueOrExpression]]:
    """solve the foreground inventory. Retrurn dict o"""

    A = fg_matrix.to_sympy()
    B = bg_matrix.to_sympy()

    # BG
    if len(bg_matrix.cols_acts()) == 0:
        # Case of empty matrix
        res_mat = ImmutableMatrix([[]])
    else:
        res_mat = (A**-1) * B

    # Transform to dict of dict
    res = dict()
    for i_fg, fg_act in enumerate(fg_matrix.cols_acts()):
        row = dict()
        res[fg_act] = row
        for i_bg, bg_act in enumerate(bg_matrix.cols_acts()):
            val = res_mat[i_fg, i_bg]
            if val == 0:
                continue
            row[bg_act] = _clean_expr(val)

    return res


def _computeDbExpressions(db_name, axis=None) -> Dict[Activity, Dict[Activity, ValueOrExpression]]:
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

    # Fill the matrices
    def fill_matrices_rec(act: Activity, axis_tag=None):
        # Update axis values
        if axis is not None:
            new_axis_val = _get_axis(act, axis)
            if new_axis_val is not None:
                axis_tag = new_axis_val

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
                fill_matrices_rec(sub_act, axis_tag)

            else:
                # Only tag bg values
                if axis_tag is not None:
                    amount = AxisDict({axis_tag: amount})

                bg_matrix[act, sub_act] += amount

    # Fill the matrices exploring everything
    for act in bw.Database(db_name):
        fill_matrices_rec(act)

    # solve the linear equation algebically
    return _solve_expression(fg_matrix, bg_matrix)


def _reverse_dict(dic):
    return {v: k for k, v in dic.items()}
