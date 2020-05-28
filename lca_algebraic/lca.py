from collections import OrderedDict
import concurrent.futures

from sympy import lambdify, simplify

from .base_utils import _actName, _eprint, _getDb
from .base_utils import _getAmountOrFormula
from .helpers import *
from .params import _param_registry, _completeParamValues, _fixed_params


def _multiLCA(activities, methods):
    """Simple wrapper around brightway API"""
    bw.calculation_setups['process'] = {'inv': activities, 'ia': methods}
    lca = bw.MultiLCA('process')
    cols = [_actName(act) for act_amount in activities for act, amount in act_amount.items()]
    return pd.DataFrame(lca.results.T, index=[method_name(method) for method in methods], columns=cols)


def multiLCA(model, methods, **params):
    """Compute LCA for a single activity and a set of methods, after settings the parameters and updating exchange amounts.

    Parameters
    ----------
    model : Single activity (root model) or list of activities
    methods : Impact methods to consider
    params : Other parameters of the model
    """

    # Check and expand params
    params = _completeParamValues(params)

    # Update brightway parameters
    bwParams = [dict(name=key, amount=value) for key, value in params.items()]
    bw.parameters.new_project_parameters(bwParams)

    # ActivityParameter.recalculate_exchanges(DEFAULT_PARAM_GROUP)
    bw.parameters.recalculate()

    if isinstance(model, list):
        activities = [{act: 1} for act in model]
    else:
        activities = [{model: 1}]
    return _multiLCA(activities, methods).transpose()


def _modelToExpr(model: ActivityExtended, methods, extraFixedParams=None, fixed_mode=FixedParamMode.DEFAULT) :
    '''
    Compute expressions corresponding to a model for each impact, replacing activities by the value of its impact


    Return
    ------
    <list of expressions (one per impact)>, <list of required param names>

    '''
    dbname = model.key[0]

    # print("computing model to expression for %s" % model)
    expr, actBySymbolName = actToExpression(model, extraFixedParams=extraFixedParams, fixed_mode=fixed_mode)

    # Required params
    free_names = set([str(symb) for symb in expr.free_symbols])
    act_names = set([str(symb) for symb in actBySymbolName.keys()])
    expected_names = free_names - act_names

    # If we expect an enu param name, we also expect the other ones : enumparam_val1 => enumparam_val1, enumparam_val2, ...
    expected_names = _expand_param_names(_expanded_names_to_names(expected_names))

    # Create dummy reference to biosphere
    # We cannot run LCA to biosphere activities
    # We create a technosphere activity mapping exactly to 1 biosphere item
    pureTechActBySymbol = OrderedDict()
    for name, act in actBySymbolName.items():
        if act[0] == BIOSPHERE3_DB_NAME:
            act = _getOrCreateDummyBiosphereActCopy(dbname, act[1])
        else:
            act = getActByCode(*act)
        pureTechActBySymbol[name] = act

    # List of activities, ordered
    acts = pureTechActBySymbol.values()

    # Transform to [{act1:1], {act2:1}, etc] for MultiLCA
    actsWithAmount = [{act: 1} for act in acts]

    # Compute LCA for all background activities and methods
    lca = _multiLCA(actsWithAmount, methods)

    # For each method, compute an algebric expression with activities replaced by their values
    exprs = []
    for imethod, method in enumerate(methods):
        # print("Generating lamba function for %s / %s" % (model, method))

        # Replace activities by their value in expression for this method
        sub = dict({symbol: lca.iloc[imethod, iact] for iact, symbol in enumerate(pureTechActBySymbol.keys())})
        exprs.append(expr.xreplace(sub))

    return exprs, expected_names

def simplifiedModel(model, impacts, extraFixedParams=None, fixed_mode=FixedParamMode.DEFAULT) :
    '''
    Computes simplified expressions corresponding to a model for each impact, replacing activities by the value of its impact.

    Return
    ------
    list of sympy expressions (one per impact)

    Parameters
    ----------
    extraFixedParams : List of extra parameters to fix
    '''
    exprs, expected_names = _modelToExpr(model, impacts, extraFixedParams=extraFixedParams, fixed_mode=fixed_mode)

    return [simplify(ex) for ex in exprs]



class LambdaWithParamNames :
    """
    This class represents a compiled (lambdified) expression together with the list of requirement parameters and the source expression
    """
    def __init__(self, expr, expanded_params=None, params=None):
        """ Computes a lamdda function from expression and list of expected parameters.
        you can provide either the list pf expanded parameters (full vars for enums) for the 'user' param names
        """
        self.expr = expr
        self.params = params
        if expanded_params is None :
            expanded_params = _expand_param_names(params)
        if self.params is None :
            self.params = _expanded_names_to_names(expanded_params)

        self.lambd = lambdify(expanded_params, expr, 'numpy')
        self.expanded_params = expanded_params

    def compute(self, **params):
        """Compute result value based of input parameters """
        # Filter on required parameters
        params = {key : val for key, val in params.items() if key in self.expanded_params}
        return self.lambd(**params)

def preMultiLCAAlgebric(model: ActivityExtended, methods):
    '''
        This method transforms an activity into a set of functions ready to compute LCA very fast on a set on methods.
        You may use is and pass the result to postMultiLCAAlgebric for fast computation on a model that does not change.

        This method is used by multiLCAAlgebric
    '''
    exprs, expected_names = _modelToExpr(model, methods)

    # Lambdify (compile) expressions
    return [LambdaWithParamNames(expr, expected_names) for expr in exprs]


def method_name(method):
    return method[1] + " - " + method[2]

def _slugify(str) :
    return re.sub('[^0-9a-zA-Z]+', '_', str)

def postMultiLCAAlgebric(methods, lambdas, alpha=1, **params):
    '''
        Compute LCA for a given set of parameters and pre-compiled lambda functions.
        This function is used by **multiLCAAlgebric**

        Parameters
        ----------
        methodAndLambdas : Output of preMultiLCAAlgebric
        **params : Parameters of the model
    '''

    # Check length of parameter values
    param_length = 1
    for key, val in params.items():
        if isinstance(val, list):
            if param_length == 1:
                param_length = len(val)
            elif param_length != len(val):
                raise Exception("Parameters should be a single value or a list of same number of values")

    def complete_params(required_param_names) :

        # Check and expand enum params, add default values
        res = _completeParamValues(params, required_param_names)

        # Expand single params and transform them to np.array
        for key in res.keys():
            val = res[key]
            if not isinstance(val, list):
                val = list([val] * param_length)
            res[key] = np.array(val, float)
        return res

    # Init output
    res = np.zeros((len(methods), param_length), float)

    # Compute result on whole vectors of parameter samples at a time : lambdas use numpy for vector computation
    def process(args) :
        imethod = args[0]
        lambd_with_params : LambdaWithParamNames = args[1]

        completed_params = complete_params(lambd_with_params.params)

        value = alpha * lambd_with_params.compute(**completed_params)
        return (imethod, value)

    # Use multithread for that
    with concurrent.futures.ThreadPoolExecutor() as exec:
        for imethod, value in exec.map(process, enumerate(lambdas)):
            res[imethod, :] = value

    return pd.DataFrame(res, index=[method_name(method) for method in methods]).transpose()


def _expand_param_names(param_names):
    '''Expand parameters names (with enum params) '''
    return [name for key in param_names for name in _param_registry()[key].names()]


def _expanded_names_to_names(param_names):
    """Find params corresponding to expanded names, including enums."""
    param_names = set(param_names)
    res = dict()
    for param in _param_registry().values():
        for name in param.names():
            if name in param_names:
                res[name] = param

    missing = param_names - set(res.keys())
    if len(missing) > 0:
        raise Exception("Unkown params : %s" % missing)

    return {param.name for param in res.values()}

# Add default values for issing parameters or warn about extra params
def _filter_params(params, expected_names, model) :
    res = params.copy()

    expected_params_names = _expanded_names_to_names(expected_names)
    for expected_name in expected_params_names:
        if expected_name not in params:
            default = _param_registry()[expected_name].default
            res[expected_name] = default
            _eprint("Missing parameter %s, replaced by default value %s" % (expected_name, default))

    for key, value in params.items():
        if not key in expected_params_names:
            del res[key]
            if model :
                _eprint("Param %s not required for model %s" % (key, model))
    return res

def multiLCAAlgebric(models, methods, **params):
    """Compute LCA by expressing the foreground model as symbolic expression of background activities and parameters.
    Then, compute 'static' inventory of the referenced background activities.
    This enables a very fast recomputation of LCA with different parameters, useful for stochastic evaluation of parametrized model

    Parameters
    ----------
    models : Single model or list of models or dict of model:amount : if list of models, you cannot use param lists
    methods : List of methods / impacts to consider
    params : You should provide named values of all the parameters declared in the model. \
             Values can be single value or list of samples, all of the same size
    """
    dfs = dict()

    if not isinstance(models, list):
        models = [models]

    for model in models:

        alpha = 1
        if type(model) is tuple:
            model, alpha = model

        # Fill default values
        lambdas = preMultiLCAAlgebric(model, methods)

        # Filter on required parameters
        df = postMultiLCAAlgebric(methods, lambdas, alpha=alpha, **params)

        model_name = _actName(model)

        # Single params ? => give the single row the name of the model activity
        if df.shape[0] == 1:
            df = df.rename(index={0: model_name})

        dfs[model_name] = df

    if len(dfs) == 1:
        df = list(dfs.values())[0]
        return df
    else:
        # Concat several dataframes for several models
        return pd.concat(list(dfs.values()))


def _getOrCreateDummyBiosphereActCopy(dbname, code):
    """
        We cannot reference directly biosphere in the model, since LCA can only be applied to products
        We create a dummy activity in our DB, with same code, and single exchange of amount '1'
    """

    code_to_find = code + "#asTech"
    try:
        return _getDb(dbname).get(code_to_find)
    except:
        bioAct = _getDb(BIOSPHERE3_DB_NAME).get(code)
        name = bioAct['name'] + ' # asTech'
        res = newActivity(dbname, name, bioAct['unit'], {bioAct: 1}, code=code_to_find)
        return res


def actToExpression(act: Activity, extraFixedParams=None, fixed_mode=FixedParamMode.DEFAULT):

    """Computes a symbolic expression of the model, referencing background activities and model parameters as symbols

    Returns
    -------
        (sympy_expr, dict of symbol => activity)
    """

    act_symbols = dict()  # Dict of  act = > symbol

    def act_to_symbol(db_name, code):

        act = _getDb(db_name).get(code)
        name = act['name']
        base_slug = _slugify(name)

        slug = base_slug
        i = 1
        while symbols(slug) in act_symbols.values():
            slug = f"{base_slug}{i}"
            i += 1

        return symbols(slug)

    def rec_func(act: Activity):

        res = 0
        outputAmount = 1

        for exch in act.exchanges():

            formula = _getAmountOrFormula(exch)

            if isinstance(formula, types.FunctionType):
                # Some amounts in EIDB are functions ... we ignore them
                continue

            input_db, input_code = exch['input']

            #  Different output ?
            if exch['input'] == exch['output']:
                if exch['amount'] != 1:
                    outputAmount = exch['amount']
                continue

            # Background DB => reference it as a symbol
            if input_db in [BIOSPHERE3_DB_NAME, ECOINVENT_DB_NAME()]:
                if not (input_db, input_code) in act_symbols:
                    act_symbols[(input_db, input_code)] = act_to_symbol(input_db, input_code)
                act_expr = act_symbols[(input_db, input_code)]

            # Our model : recursively it to a symbolic expression
            else:

                if input_db == act['database'] and input_code == act['code']:
                    raise Exception("Recursive exchange : %s" % (act.__dict__))

                sub_act = _getDb(input_db).get(input_code)
                act_expr = rec_func(sub_act)

            res += formula * act_expr

        return res / outputAmount

    expr = rec_func(act)

    fixed_params = set(_fixed_params().values())
    if extraFixedParams :
        fixed_params |= set(extraFixedParams)

    # Replace fixed values
    sub = {symbols(key): val for param in fixed_params for key, val in param.expandParams(param.stat_value(fixed_mode)).items()}


    expr = expr.xreplace(sub)

    return (expr, _reverse_dict(act_symbols))


def _reverse_dict(dic):
    return {v: k for k, v in dic.items()}
