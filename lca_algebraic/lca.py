import concurrent.futures
from collections import OrderedDict
from typing import Tuple, Dict

from pandas import DataFrame
from sympy import lambdify, simplify

from .base_utils import _actName, _getDb, _method_unit
from .base_utils import _getAmountOrFormula
from .helpers import *
from .helpers import _isForeground
from .params import _param_registry, _completeParamValues, _fixed_params, _expanded_names_to_names, _expand_param_names
from ipywidgets import Output


def _impact_labels():
    """Dictionnary of custom impact names
        Dict of "method tuple" => string
    """
    # Prevent reset upon auto reload in jupyter notebook
    if not '_impact_labels' in builtins.__dict__:
        builtins._impact_labels = dict()

    return builtins._impact_labels

def set_custom_impact_labels(impact_labels:Dict) :
    """ Global function to override name of impact method in graphs """
    _impact_labels().update(impact_labels)

def _multiLCA(activities, methods):
    """Simple wrapper around brightway API"""
    bw.calculation_setups['process'] = {'inv': activities, 'ia': methods}
    lca = bw.MultiLCA('process')
    cols = [_actName(act) for act_amount in activities for act, amount in act_amount.items()]
    return pd.DataFrame(lca.results.T, index=[method_name(method) for method in methods], columns=cols)


def multiLCA(models, methods, **params):
    """Compute LCA for a single activity and a set of methods, after settings the parameters and updating exchange amounts.
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
    for db in dbs :
        if _isForeground(db) :
            freezeParams(db, **params)

    activities = [{act: 1} for act in models]
    return _multiLCA(activities, methods).transpose()



# Cache of (act, method) => values
_BG_IMPACTS_CACHE = dict()

def _clearLCACache() :
    _BG_IMPACTS_CACHE.clear()

""" Compute LCA and return (act, method) => value """
def _multiLCAWithCache(acts, methods) :

    # List activities with at least one missing value
    remaining_acts = list(act for act in acts if any(method for method in methods if (act, method) not in _BG_IMPACTS_CACHE))

    if (len(remaining_acts) > 0) :
        debug("remaining act", remaining_acts, methods)
        lca = _multiLCA(
            [{act: 1} for act in remaining_acts],
            methods)

        # Set output from dataframe
        for imethod, method in enumerate(methods) :
            for iact, act in enumerate(remaining_acts) :
                _BG_IMPACTS_CACHE[(act, method)] = lca.iloc[imethod, iact]

    return _BG_IMPACTS_CACHE





def _filter_param_values(params, expanded_param_names) :
    return {key : val for key, val in params.items() if key in expanded_param_names}

class LambdaWithParamNames :
    """
    This class represents a compiled (lambdified) expression together with the list of requirement parameters and the source expression
    """
    def __init__(self, exprOrDict, expanded_params=None, params=None, sobols=None):
        """ Computes a lamdda function from expression and list of expected parameters.
        you can provide either the list pf expanded parameters (full vars for enums) for the 'user' param names
        """

        if isinstance(exprOrDict, dict) :
            # Come from JSON serialization
            obj = exprOrDict
            # Parse
            self.params = obj["params"]

            # Full names
            self.expanded_params = _expand_param_names(self.params)
            self.expr = parse_expr(obj["expr"])
            self.lambd = lambdify(self.expanded_params, self.expr, 'numpy')
            self.sobols = obj["sobols"]

        else :
            self.expr = exprOrDict
            self.params = params

            if expanded_params is None :
                expanded_params = _expand_param_names(params)
            if self.params is None :
                self.params = _expanded_names_to_names(expanded_params)

            self.lambd = lambdify(expanded_params, exprOrDict, 'numpy')
            self.expanded_params = expanded_params
            self.sobols = sobols

    def complete_params(self, params):

        param_length = _compute_param_length(params)

        # Check and expand enum params, add default values
        res = _completeParamValues(params, self.params)

        # Expand single params and transform them to np.array
        for key in res.keys():
            val = res[key]
            if not isinstance(val, list):
                val = list([val] * param_length)
            res[key] = np.array(val, float)
        return res

    def compute(self, **params):
        """Compute result value based of input parameters """
        # Filter on required parameters
        params = _filter_param_values(params, self.expanded_params)
        return self.lambd(**params)

    def serialize(self) :
        return dict(
            params=self.params,
            expr=str(self.expr),
            sobols=self.sobols)

    def __repr__(self):
        return repr(self.expr)

    def _repr_latex_(self):
        return self.expr._repr_latex_()

MethodId = Tuple[str, str, str]

def _modelsToLambdas(models : List[ActivityExtended], methods : List[MethodId]) -> Dict[Tuple[ActivityExtended, MethodId], LambdaWithParamNames]:
    """
    Compile list of models (=activities) and methods (=impacts) to SymPy expression,
    replacing the background activities by te values of the impacts

    :param models: List of activitiies
    :param methods: List of impact methods (tuples)
    :return: dict of (model, method) => LambdaWithParamNames
    """

    # Cache of (db, key) => symbol name
    symbol_by_act = dict()
    expr_by_act = dict()
    res = dict()

    # Generate expressions for all models
    for model in models :
        expr_by_act[model] = actToExpression(model, symbol_by_act)

    # Create tech proxys if necessary (for bio flux)
    acts_by_name = {symbol: _createTechProxyForBio(act, act[0]) for act, symbol in symbol_by_act.items()}

    # Compute LCA for all background activities
    lcas = _multiLCAWithCache(
        acts_by_name.values(),
        methods)

    # Loop on models
    for model in models :

        with DbContext(model) :

            expr = expr_by_act[model]

            # Compute the required params
            free_names = set([str(symb) for symb in expr.free_symbols])
            act_names = set([str(symb) for symb in symbol_by_act.values()])
            expected_names = free_names - act_names

            # If we expect an enum param name, we also expect the other ones : enumparam_val1 => enumparam_val1, enumparam_val2, ...
            expected_names = _expand_param_names(_expanded_names_to_names(expected_names))

            # Loop on impact methods and replace background activities by their values
            for method in methods:

                # Replace activities by their value in expression for this method
                sub = dict({symbol: lcas[(act, method)] for symbol, act in acts_by_name.items()})

                res[(model, method)] = LambdaWithParamNames(
                    expr.xreplace(sub),
                    expected_names)

    return res


def _modelToLambdas(model: ActivityExtended, methods):
    """
    Wraps _modelsToLambdas for a single model, returning a list of lambda function in the same order as methods.
    """
    with DbContext(model) :
        res = _modelsToLambdas([model], methods)
        return list(res[(model, method)] for method in methods)

def method_name(method):
    """Return name of method, taking into account custom label set via set_custom_impact_labels(...) """
    if method in _impact_labels() :
        return _impact_labels()[method]
    return method[1] + " - " + method[2]

def _slugify(str) :
    return re.sub('[^0-9a-zA-Z]+', '_', str)

def _compute_param_length(params) :
    # Check length of parameter values
    param_length = 1
    for key, val in params.items():
        if isinstance(val, list):
            if param_length == 1:
                param_length = len(val)
            elif param_length != len(val):
                raise Exception("Parameters should be a single value or a list of same number of values")
    return param_length


def _applyParams(lambd, params) :
    """
        Apply given params to lambda function.
        return either a single value of numpy vector
    """
    completed_params = lambd.complete_params(params)
    res =  lambd.compute(**completed_params)
    if isinstance(res, (list, np.ndarray)) and len(res) == 1:
        return res[0]
    return res


def _postMultiLCAAlgebric(methods, lambdas, **params):
    '''
        Compute LCA for a given set of parameters and pre-compiled lambda functions.
        This function is used by **multiLCAAlgebric**.

        Returns a dataframe.

        Parameters
        ----------
        methods: List of impact methods
        lambdas: dict of lambda:alpha
        **params : Parameters of the model
    '''

    param_length = _compute_param_length(params)

    # Init output
    res = np.zeros((len(methods), param_length), float)

    # Compute result on whole vectors of parameter samples at a time : lambdas use numpy for vector computation
    def process(args):
        imethod, lambd = args
        return (imethod, _applyParams(lambd, params))

    # Use multithread for that
    with concurrent.futures.ThreadPoolExecutor() as exec:
        for imethod, value in exec.map(process, enumerate(lambdas)):
            res[imethod, :] = value

    return pd.DataFrame(
        res,
        index=[method_name(method) + "[%s]" % _method_unit(method) for method in methods]).transpose()


# Add default values for issing parameters or warn about extra params
def _filter_params(params, expected_names, model) :
    res = params.copy()

    expected_params_names = _expanded_names_to_names(expected_names)
    for expected_name in expected_params_names:
        if expected_name not in params:
            default = _param_registry()[expected_name].default
            res[expected_name] = default
            error("Missing parameter %s, replaced by default value %s" % (expected_name, default))

    for key, value in params.items():
        if not key in expected_params_names:
            del res[key]
            if model :
                error("Param %s not required for model %s" % (key, model))
    return res

def multiLCAAlgebric(
        models,
        methods,
        raw=False, # If raw if true, returns dict of [(model, method)] -> value(s). Otherwize, returns Dataframe
        **params):
    """
    Main parametric LCIA method : Computes LCA by expressing the foreground model as symbolic expression of background activities and parameters.
    Then, compute 'static' inventory of the referenced background activities.
    This enables a very fast recomputation of LCA with different parameters, useful for stochastic evaluation of parametrized model

    Parameters
    ----------
    models :
        Single model or
        List of model or
        List of (model, alpha)
        or Dict of model:amount
        In case of several models, you cannot use list of parameters
    methods : List of methods / impacts to consider
    extract_activities : Optionnal : list of foregound or background activities. If provided, the result only integrate their contribution
    params : You should provide named values of all the parameters declared in the model. \
             Values can be single value or list of samples, all of the same size
    """
    if not isinstance(methods, list) :
        methods = [methods]

    if isinstance(models, list):
        def to_tuple(item) :
            if isinstance(item, tuple) :
                return item
            else:
                return (item, 1)
        models = dict(to_tuple(item) for item in models)
    elif not isinstance(models, dict):
        models = {models:1}

    for model in models.keys():
        dbname = model.key[0]
        with DbContext(dbname):

            # Check no params are passed for FixedParams
            for key in params:
                if key in _fixed_params() :
                    error("Param '%s' is marked as FIXED, but passed in parameters : ignored" % key)

            debug("computing : ", model)

    lambdas_per_model_method = _modelsToLambdas(models.keys(), methods)

    param_length = _compute_param_length(params)

    res = defaultdict(dict)

    # Choose shape of the dataframe depending of size of each axe
    if param_length > 1 :
        row_axe = "params"
        if len(models) > 1 :
            if len(methods) > 1 :
                raise Exception("Only two of three dimensions (models, methods, param list) can be larger than 1")
            else:
                col_axe = "models"
        else:
            col_axe = "methods"
    else :
        col_axe = "methods"
        row_axe = "models"


    for method in methods :

        method_n = method_name(method) + "[%s]" % _method_unit(method)

        for model, alpha in models.items() :
            with DbContext(model) :

                model_n = _actName(model)

                lambd = lambdas_per_model_method[(model, method)]
                values = _applyParams(lambd, params) * alpha

                if raw :
                    res[(model, method)] = values
                else:
                    col = method_n if col_axe == "methods" else model_n

                    if row_axe == "params" :
                        res[col] = values
                    else:
                        row = method_n if row_axe == "methods" else model_n
                        res[col][row] = values

    if raw :
        return res
    else:
        return DataFrame.from_dict(res)


def _createTechProxyForBio(act_key, target_db):
    """
        We cannot reference directly biosphere in the model, since LCA can only be applied to products
        We create a dummy activity in our DB, with same code, and single exchange of amount '1'
    """
    dbname, code = act_key
    act = _getDb(dbname).get(code)

    # Biosphere ?
    if (dbname == BIOSPHERE3_DB_NAME) or ("type" in act and act["type"] in ["emission", "natural resource"]) :

        code_to_find = code + "#asTech"

        try:
            # Already created ?
            return _getDb(target_db).get(code_to_find)
        except:
            name = act['name'] + ' # asTech'

            # Create biosphere proxy in User Db
            res = newActivity(target_db, name, act['unit'], {act: 1},
                              code=code_to_find,
                              isProxy=True) # add a this flag to distinguish this dummy activity from others
            return res
    else :
        return act


def _replace_fixed_params(expr, fixed_params, fixed_mode=FixedParamMode.DEFAULT) :
    """Replace fixed params with their value."""
    sub = {symbols(key): val for param in fixed_params for key, val in param.expandParams(param.stat_value(fixed_mode)).items()}
    return expr.xreplace(sub)

@with_db_context(arg="act")
def actToExpression(
        act: Activity,
        act_symbols : Dict[str, Symbol] = dict()):

    """
    Computes a symbolic expression of the model, referencing background activities and model parameters as symbols
    (sympy_expr, dict of symbol => activity)
    :param act: Activity
    :param act_symbols: Cache of activity -> Symbol (filled during this process)
    :return: Sympy expression
    """

    def act_to_symbol(sub_act):
        """ Transform an activity to a named symbol and keep cache of it """

        db_name, code = sub_act.key

        # Look in cache
        if not (db_name, code) in act_symbols:

            act = _getDb(db_name).get(code)
            name = act['name']
            base_slug = _slugify(name)

            slug = base_slug
            i = 1
            while symbols(slug) in act_symbols.values():
                slug = f"{base_slug}{i}"
                i += 1

            act_symbols[(db_name, code)] = symbols(slug)

        return act_symbols[(db_name, code)]

    def rec_func(act: Activity, parents):

        res = 0
        outputAmount = act.getOutputAmount()

        if not _isForeground(act["database"]) :
            # We reached a background DB ? => stop developping and create reference to activity
            return act_to_symbol(act)

        for exch in act.exchanges():

            formula = _getAmountOrFormula(exch)

            if isinstance(formula, types.FunctionType):
                # Some amounts in EIDB are functions ... we ignore them
                continue

            #  Production exchange
            if exch['input'] == exch['output']:
                continue

            input_db, input_code = exch['input']
            sub_act = _getDb(input_db).get(input_code)


            # Background DB => reference it as a symbol
            if not _isForeground(input_db) :

                # Add to dict of background symbols
                act_expr = act_to_symbol(sub_act)

            # Our model : recursively it to a symbolic expression
            else:

                parents = parents + [act]
                if sub_act in parents :
                    raise Exception("Found recursive activities : " + ", ".join(_actName(act) for act in (parents + [sub_act])))

                act_expr = rec_func(sub_act, parents)

            avoidedBurden = 1

            if exch.get('type') == 'production' and not exch.get('input') == exch.get('output') :
                debug("Avoided burden", exch[name])
                avoidedBurden = -1

            #debug("adding sub act : ", sub_act, formula, act_expr)

            res += formula * act_expr * avoidedBurden

        return res / outputAmount

    expr = rec_func(act, [])

    if isinstance(expr, float) :
        expr = simplify(expr)
    else:
        # Replace fixed params with their default value
        expr = _replace_fixed_params(expr, _fixed_params().values())

    return expr


def _reverse_dict(dic):
    return {v: k for k, v in dic.items()}



