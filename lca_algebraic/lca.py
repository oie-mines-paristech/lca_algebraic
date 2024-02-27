import builtins
import concurrent.futures
from collections import OrderedDict
from dataclasses import dataclass
from typing import List
from warnings import warn

import numpy as np
import pandas as pd
from peewee import DoesNotExist
from sympy import Symbol, lambdify, parse_expr

from . import TabbedDataframe
from .axis_dict import AxisDict
from .base_utils import _actName, _getDb, _method_unit
from .cache import ExprCache, LCIACache
from .helpers import (
    BIOSPHERE_PREFIX,
    Activity,
    ActivityExtended,
    Basic,
    DbContext,
    Dict,
    _actDesc,
    _getAmountOrFormula,
    _isForeground,
    bw,
    error,
    newActivity,
    re,
    simplify,
    symbols,
    types,
    with_db_context,
)
from .log import logger
from .params import (
    FixedParamMode,
    _complete_params,
    _compute_param_length,
    _expand_param_names,
    _expand_params,
    _expanded_names_to_names,
    _fixed_params,
    _param_registry,
    _toSymbolDict,
    all_params,
    compute_expr_value,
    freezeParams,
)


def _impact_labels():
    """Dictionnary of custom impact names
    Dict of "method tuple" => string
    """
    # Prevent reset upon auto reload in jupyter notebook
    if "_impact_labels" not in builtins.__dict__:
        builtins._impact_labels = dict()

    return builtins._impact_labels


def set_custom_impact_labels(impact_labels: Dict):
    """Global function to override name of impact method in graphs"""
    _impact_labels().update(impact_labels)


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


def _multiLCAWithCache(acts, methods):
    with LCIACache() as cache:
        # List activities with at least one missing value
        remaining_acts = list(act for act in acts if any(method for method in methods if (act, method) not in cache.data))

        if len(remaining_acts) > 0:
            lca = _multiLCA([{act: 1} for act in remaining_acts], methods)

            # Set output from dataframe
            for imethod, method in enumerate(methods):
                for iact, act in enumerate(remaining_acts):
                    cache.data[(act, method)] = lca.iloc[imethod, iact]

        # Return a copy of the cache for selected impacts and activities
        return {(act, method): cache.data[(act, method)] for act in acts for method in methods}


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


def _modelToExpr(model: ActivityExtended, methods, alpha=1, axis=None):
    """
    Compute expressions corresponding to a model for each impact, replacing activities by the value of its impact

    Return
    ------
    <list of expressions (one per impact)>, <list of required param names>

    """

    # Try to load from cache
    with ExprCache() as cache:
        key = (model, axis)

        if key not in cache.data:
            logger.debug(f"{model} was not in exrepssion cache, computing...")

            expr, actBySymbolName = actToExpression(model, axis=axis)

            cache.data[key] = (expr, actBySymbolName)
        else:
            logger.debug(f"{model} found in exrepssion cache")

        expr, actBySymbolName = cache.data[key]

    logger.debug("Alpha passed %s", alpha)

    # logger.debug("Raw expression for %s/%s : '%s'", model, str(methods), expr)
    # logger.debug("Act by symbol : %s", actBySymbolName)
    # if logger.isEnabledFor("DEBUG") :
    # logger.debug(f"Length of expression : {len(str(expr))}")

    expr = expr * alpha

    # Create dummy reference to biosphere
    # We cannot run LCA to biosphere activities
    # We create a technosphere activity mapping exactly to 1 biosphere item
    pureTechActBySymbol = OrderedDict()
    for act_name, act in actBySymbolName.items():
        pureTechActBySymbol[act_name] = _createTechProxyForBio(act, model.key[0])

    # Compute LCA for background activities
    lcas = _multiLCAWithCache(pureTechActBySymbol.values(), methods)

    # For each method, compute an algebric expression with activities replaced by their values
    exprs = []
    for method in methods:
        # Replace activities by their value in expression for this method
        sub = dict({symbol: lcas[(act, method)] for symbol, act in pureTechActBySymbol.items()})

        expr_curr = expr.xreplace(sub)

        # Ensure symbols are params
        expr_curr = _replace_symbols_with_params_in_exp(expr_curr)

        exprs.append(expr_curr)

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
    if isinstance(expr, Basic):
        lambd = lambdify(expanded_params, expr, "numpy")

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


def compute_param_formulas(param_values: Dict, required_params: List[str]):
    # Compute the values params computed from formulas (and not specified in input)
    pass


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

    def __init__(self, exprOrDict, expanded_params=None, params=None, sobols=None):
        """Computes a lamdda function from expression and list of expected parameters.
        you can provide either the list pf expanded parameters (full vars for enums) for the 'user' param names
        """

        if isinstance(exprOrDict, dict):
            # Come from JSON serialization
            obj = exprOrDict
            # LIst of required params for this lambda
            self.params: List[str] = obj["params"]

            # Full names
            self.expanded_params = _expand_param_names(self.params)
            self.expr = parse_expr(obj["expr"])
            self.lambd = _lambdify(self.expr, self.expanded_params)
            self.sobols = obj["sobols"]

        else:
            self.expr = exprOrDict
            self.params = params

            if expanded_params is None:
                if params is None:
                    expanded_params = _free_symbols(exprOrDict)
                    params = _expanded_names_to_names(expanded_params)
                    self.params = params

                # We expand again the parameters
                # If we expect an enum param name, we also expect the other ones :
                # enumparam_val1 => enumparam_val1, enumparam_val2, ...
                expanded_params = _expand_param_names(params)

            elif self.params is None:
                self.params = _expanded_names_to_names(expanded_params)

            self.lambd = _lambdify(exprOrDict, expanded_params)
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

    def __repr__(self):
        return repr(self.expr)

    def _repr_latex_(self):
        return self.expr._repr_latex_()


def _preMultiLCAAlgebric(model: ActivityExtended, methods, alpha=1, axis=None):
    """
    This method transforms an activity into a set of functions ready to compute LCA very fast on a set on methods.
    You may use is and pass the result to postMultiLCAAlgebric for fast computation on a model that does not change.

    This method is used by multiLCAAlgebric
    """
    with DbContext(model):
        exprs = _modelToExpr(model, methods, alpha=alpha, axis=axis)

        # Lambdify (compile) expressions
        return [LambdaWithParamNames(expr) for expr in exprs]


def method_name(method):
    """Return name of method, taking into account custom label set via set_custom_impact_labels(...)"""
    if method in _impact_labels():
        return _impact_labels()[method]
    return method[1] + " - " + method[2]


def _slugify(str):
    return re.sub("[^0-9a-zA-Z]+", "_", str)


@dataclass
class ResultsWithParams:
    """Holds bith the result with context parameters"""

    dataframe: pd.DataFrame
    params: Dict


def _postMultiLCAAlgebric(methods, lambdas: List[LambdaWithParamNames], with_params=False, **params):
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
            value = list(float(val) for val in value.values())

        return (imethod, value)

    # Use multithread for that
    with concurrent.futures.ThreadPoolExecutor() as exec:
        for imethod, value in exec.map(process, enumerate(lambdas)):
            res[imethod, :] = value

    result = pd.DataFrame(
        res,
        index=[method_name(method) + "[%s]" % _method_unit(method) for method in methods],
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
            error("Missing parameter %s, replaced by default value %s" % (expected_name, default))

    for key, value in params.items():
        if key not in expected_params_names:
            del res[key]
            if model:
                error("Param %s not required for model %s" % (key, model))
    return res


def compute_value(formula, **params):
    """Compute actual value for a given formula, with possible parameters (or default ones)"""
    if isinstance(formula, float) or isinstance(formula, int):
        return formula

    lambd = LambdaWithParamNames(formula)

    value_context = lambd.compute(**params)

    return value_context.value


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
    methods : List of methods / impacts to consider
    params : You should provide named values of all the parameters declared in the model. \
             Values can be single value or list of samples, all of the same size
    axis: Designates the name of an attribute of user activities to split impacts by their value. \
        This is useful to get impact by phase or sub modules
    functional_unit: quantity (static or Sypy formula) by which to divide impacts
    return_params: If true, also returns the value of all parameters in as tabbed DataFrame
    description: Optional description/metadata to be added in output when using "return params" Dataframe
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

    for model, alpha in models.items():
        if type(model) is tuple:
            model, alpha = model

        alpha = float(alpha)

        dbname = model.key[0]
        with DbContext(dbname):
            # Check no params are passed for FixedParams
            for key in params:
                if key in _fixed_params():
                    error("Param '%s' is marked as FIXED, but passed in parameters : ignored" % key)

            if functional_unit != 1:
                alpha = alpha / functional_unit

            lambdas = _preMultiLCAAlgebric(model, methods, alpha=alpha, axis=axis)

            res = _postMultiLCAAlgebric(methods, lambdas, with_params=return_params, **params)
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
                df = df.rename(index={None: "*other*"})

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


def _createTechProxyForBio(act_key, target_db):
    """
    We cannot reference directly biosphere in the model, since LCA can only be applied to products
    We create a dummy activity in our DB, with same code, and single exchange of amount '1'
    """
    dbname, code = act_key
    act = _getDb(dbname).get(code)

    # Biosphere ?
    if (BIOSPHERE_PREFIX in dbname) or ("type" in act and act["type"] in ["emission", "natural resource"]):
        code_to_find = code + "#asTech"

        try:
            # Already created ?
            return _getDb(target_db).get(code_to_find)
        except DoesNotExist:
            name = act["name"] + " # asTech"

            # Create biosphere proxy in User Db
            res = newActivity(
                target_db, name, act["unit"], {act: 1}, code=code_to_find, isProxy=True
            )  # add a this flag to distinguish this dummy activity from others
            return res
    else:
        return act


def _replace_fixed_params(expr, fixed_params, fixed_mode=FixedParamMode.DEFAULT):
    """Replace fixed params with their value."""

    sub = {key: val for param in fixed_params for key, val in param.expandParams(param.stat_value(fixed_mode)).items()}
    sub = _toSymbolDict(sub)
    return expr.xreplace(sub)


def _safe_axis(axis_name: str):
    if axis_name.isalnum():
        return axis_name
    else:
        return re.sub("[^0-9a-zA-Z]+", "_", axis_name)


def _tag_expr(expr, act, axis):
    """Tag expression for one axe. Check the child expression is not already tagged with different values"""
    axis_tag = act.get(axis, None)

    if axis_tag is None:
        return expr

    axis_tag = _safe_axis(axis_tag)

    if isinstance(expr, AxisDict):
        res = 0
        for key, val in expr._dict.items():
            if key is not None and str(key) != axis_tag:
                raise ValueError(
                    "Inconsistent axis for one change of  '%s' : attempt to tag as '%s'. "
                    "Already tagged as '%s'. Value of the exchange : %s" % (act["name"], axis_tag, key, str(val))
                )
            res += val
    else:
        res = expr

    return AxisDict({axis_tag: res})


@with_db_context(arg="act")
def actToExpression(act: Activity, axis=None):
    """Computes a symbolic expression of the model, referencing background activities and model parameters as symbols

    Returns
    -------
        (sympy_expr, dict of symbol => activity)
    """

    act_symbols: Dict[Symbol] = dict()  # Cache of  act = > symbol

    def act_to_symbol(sub_act):
        """Transform an activity to a named symbol and keep cache of it"""

        db_name, code = sub_act.key

        # Look in cache
        if not (db_name, code) in act_symbols:
            act = _getDb(db_name).get(code)
            name = act["name"]
            base_slug = _slugify(name)

            slug = base_slug
            i = 1
            while symbols(slug) in act_symbols.values():
                slug = f"{base_slug}{i}"
                i += 1

            act_symbols[(db_name, code)] = symbols(slug)

        return act_symbols[(db_name, code)]

    # Local cache of expressions

    def actToExpressionRec(act: ActivityExtended, parents=[]):
        res = 0
        outputAmount = act.getOutputAmount()

        if not _isForeground(act["database"]):
            # We reached a background DB ? => stop developping and create reference to activity
            return act_to_symbol(act)

        for exch in act.exchanges():
            formula = _getAmountOrFormula(exch)

            if isinstance(formula, types.FunctionType):
                # Some amounts in EIDB are functions ... we ignore them
                continue

            #  Production exchange
            if exch["input"] == exch["output"]:
                continue

            input_db, input_code = exch["input"]
            sub_act = _getDb(input_db).get(input_code)

            # Background DB => reference it as a symbol
            if not _isForeground(input_db):
                act_expr = act_to_symbol(sub_act)

            # Our model : recursively it to a symbolic expression
            else:
                parents = parents + [act]
                if sub_act in parents:
                    raise Exception("Found recursive activities : " + ", ".join(_actName(act) for act in (parents + [sub_act])))

                act_expr = actToExpressionRec(sub_act, parents)

            avoidedBurden = 1

            if exch.get("type") == "production" and not exch.get("input") == exch.get("output"):
                avoidedBurden = -1

            res += formula * act_expr * avoidedBurden

        res = res / outputAmount

        # Axis ? transforms this to a dict with the correct Tag
        if axis:
            res = _tag_expr(res, act, axis)

        return res

    expr = actToExpressionRec(act)

    if isinstance(expr, float):
        expr = simplify(expr)
    else:
        # Replace fixed params with their default value
        expr = _replace_fixed_params(expr, _fixed_params().values())

    return (expr, _reverse_dict(act_symbols))


def _reverse_dict(dic):
    return {v: k for k, v in dic.items()}


def exploreImpacts(impact, *activities: ActivityExtended, **params):
    """
    Advanced version of #printAct()

    Displays all exchanges of one or several activities and their impacts.
    If parameter values are provided, formulas will be evaluated accordingly.
    If two activities are provided, they will be shown side by side and compared.
    """
    tables = []
    names = []

    diffOnly = params.pop("diffOnly", False)
    withImpactPerUnit = params.pop("withImpactPerUnit", False)

    for main_act in activities:
        inputs_by_ex_name = dict()
        data = dict()

        for i, (name, input, amount) in enumerate(main_act.listExchanges()):
            # Params provided ? Evaluate formulas
            if len(params) > 0 and isinstance(amount, Basic):
                amount = compute_expr_value(amount, params)

            i = 1
            ex_name = name
            while ex_name in data:
                ex_name = "%s#%d" % (name, i)
                i += 1

            inputs_by_ex_name[ex_name] = _createTechProxyForBio(input.key, main_act.key[0])

            input_name = _actName(input)

            if _isForeground(input.key[0]):
                input_name += "{user-db}"

            data[ex_name] = dict(input=input_name, amount=amount)

        # Provide impact calculation if impact provided
        all_acts = list(set(inputs_by_ex_name.values()))
        res = multiLCAAlgebric(all_acts, [impact], **params)
        impacts = res[res.columns.tolist()[0]].to_list()

        impact_by_act = {act: value for act, value in zip(all_acts, impacts)}

        # Add impacts to data
        for key, vals in data.items():
            amount = vals["amount"]
            act = inputs_by_ex_name[key]
            impact = impact_by_act[act]

            if withImpactPerUnit:
                vals["impact_per_unit"] = impact
            vals["impact"] = amount * impact

        # To dataframe
        df = pd.DataFrame(data)

        tables.append(df.T)
        names.append(_actDesc(main_act))

    full = pd.concat(tables, axis=1, keys=names, sort=True)

    # Highlight differences in case two activites are provided
    if len(activities) == 2:
        YELLOW = "background-color:NavajoWhite"
        diff_cols = ["amount", "input", "impact"]
        cols1 = dict()
        cols2 = dict()
        for col_name in diff_cols:
            cols1[col_name] = full.columns.get_loc((names[0], col_name))
            cols2[col_name] = full.columns.get_loc((names[1], col_name))

        if diffOnly:

            def isDiff(row):
                return row[cols1["impact"]] != row[cols2["impact"]]

            full = full.loc[isDiff]

        def same_amount(row):
            res = [""] * len(row)
            for col_name in diff_cols:
                if row[cols1[col_name]] != row[cols2[col_name]]:
                    res[cols1[col_name]] = YELLOW
                    res[cols2[col_name]] = YELLOW
            return res

        full = full.style.apply(same_amount, axis=1)

    return full
