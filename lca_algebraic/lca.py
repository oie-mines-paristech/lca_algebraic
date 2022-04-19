import concurrent.futures
from collections import OrderedDict
from typing import Dict, List

from sympy import lambdify, simplify

from .base_utils import _actName, error, _getDb, _method_unit
from .base_utils import _getAmountOrFormula
from .helpers import *
from .helpers import _actDesc, _isForeground
from .params import _param_registry, _completeParamValues, _fixed_params, _expanded_names_to_names, _expand_param_names


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
        lca = _multiLCA(
            [{act: 1} for act in remaining_acts],
            methods)

        # Set output from dataframe
        for imethod, method in enumerate(methods) :
            for iact, act in enumerate(remaining_acts) :
                _BG_IMPACTS_CACHE[(act, method)] = lca.iloc[imethod, iact]

    return _BG_IMPACTS_CACHE


def _modelToExpr(
        model: ActivityExtended,
        methods,
        extract_activities=None) :
    '''
    Compute expressions corresponding to a model for each impact, replacing activities by the value of its impact

    Return
    ------
    <list of expressions (one per impact)>, <list of required param names>

    '''
    # print("computing model to expression for %s" % model)
    expr, actBySymbolName = actToExpression(
        model,
        extract_activities=extract_activities)

    # Required params
    free_names = set([str(symb) for symb in expr.free_symbols])
    act_names = set([str(symb) for symb in actBySymbolName.keys()])
    expected_names = free_names - act_names

    # If we expect an enum param name, we also expect the other ones : enumparam_val1 => enumparam_val1, enumparam_val2, ...
    expected_names = _expand_param_names(_expanded_names_to_names(expected_names))

    # Create dummy reference to biosphere
    # We cannot run LCA to biosphere activities
    # We create a technosphere activity mapping exactly to 1 biosphere item
    pureTechActBySymbol = OrderedDict()
    for name, act in actBySymbolName.items():
        pureTechActBySymbol[name] = _createTechProxyForBio(act, model.key[0])

    # Compute LCA for background activities
    lcas = _multiLCAWithCache(pureTechActBySymbol.values(), methods)

    # For each method, compute an algebric expression with activities replaced by their values
    exprs = []
    for method in methods:
        # Replace activities by their value in expression for this method
        sub = dict({symbol: lcas[(act, method)] for symbol, act in pureTechActBySymbol.items()})
        exprs.append(expr.xreplace(sub))

    return exprs, expected_names


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

def _preMultiLCAAlgebric(model: ActivityExtended, methods, extract_activities:List[Activity]=None):
    '''
        This method transforms an activity into a set of functions ready to compute LCA very fast on a set on methods.
        You may use is and pass the result to postMultiLCAAlgebric for fast computation on a model that does not change.

        This method is used by multiLCAAlgebric
    '''
    with DbContext(model) :
        exprs, expected_names = _modelToExpr(model, methods, extract_activities=extract_activities)

        # Lambdify (compile) expressions
        return [LambdaWithParamNames(expr, expected_names) for expr in exprs]


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

def _postMultiLCAAlgebric(methods, lambdas, alpha=1, **params):
    '''
        Compute LCA for a given set of parameters and pre-compiled lambda functions.
        This function is used by **multiLCAAlgebric**

        Parameters
        ----------
        methodAndLambdas : Output of preMultiLCAAlgebric
        **params : Parameters of the model
    '''

    param_length = _compute_param_length(params)

    # Init output
    res = np.zeros((len(methods), param_length), float)

    # Compute result on whole vectors of parameter samples at a time : lambdas use numpy for vector computation
    def process(args) :
        imethod = args[0]
        lambd_with_params : LambdaWithParamNames = args[1]

        completed_params = lambd_with_params.complete_params(params)

        value = alpha * lambd_with_params.compute(**completed_params)
        return (imethod, value)

    # Use multithread for that
    with concurrent.futures.ThreadPoolExecutor() as exec:
        for imethod, value in exec.map(process, enumerate(lambdas)):
            res[imethod, :] = value

    return pd.DataFrame(res, index=[method_name(method) + "[%s]" % _method_unit(method) for method in methods]).transpose()


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

def multiLCAAlgebric(models, methods, extract_activities:List[Activity]=None, **params):
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
    dfs = dict()

    if isinstance(models, list):
        def to_tuple(item) :
            if isinstance(item, tuple) :
                return item
            else:
                return (item, 1)
        models = dict(to_tuple(item) for item in models)
    elif not isinstance(models, dict):
        models = {models:1}

    for model, alpha in models.items():

        if type(model) is tuple:
            model, alpha = model

        dbname = model.key[0]
        with DbContext(dbname):

            # Check no params are passed for FixedParams
            for key in params:
                if key in _fixed_params() :
                    error("Param '%s' is marked as FIXED, but passed in parameters : ignored" % key)


            lambdas = _preMultiLCAAlgebric(model, methods, extract_activities=extract_activities)

            df = _postMultiLCAAlgebric(methods, lambdas, alpha=alpha, **params)

            model_name = _actName(model)
            while model_name in dfs :
                model_name += "'"

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
def actToExpression(act: Activity, extract_activities=None):

    """Computes a symbolic expression of the model, referencing background activities and model parameters as symbols

    Returns
    -------
        (sympy_expr, dict of symbol => activity)
    """

    act_symbols : Dict[Symbol]= dict()  # Cache of  act = > symbol

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

    def rec_func(act: Activity, in_extract_path, parents=[]):

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


            # If list of extract activites requested, we only integrate activites below a tracked one
            exch_in_path = in_extract_path
            if extract_activities is not None:
                if sub_act in extract_activities :
                    exch_in_path = in_extract_path or (sub_act in extract_activities)

            # Background DB => reference it as a symbol
            if not _isForeground(input_db) :

                if exch_in_path :
                    # Add to dict of background symbols
                    act_expr = act_to_symbol(sub_act)
                else:
                    continue

            # Our model : recursively it to a symbolic expression
            else:

                parents = parents + [act]
                if sub_act in parents :
                    raise Exception("Found recursive activities : " + ", ".join(_actName(act) for act in (parents + [sub_act])))

                act_expr = rec_func(sub_act, exch_in_path, parents)

            avoidedBurden = 1

            if exch.get('type') == 'production' and not exch.get('input') == exch.get('output') :
                debug("Avoided burden", exch[name])
                avoidedBurden = -1

            #debug("adding sub act : ", sub_act, formula, act_expr)

            res += formula * act_expr * avoidedBurden

        return res / outputAmount

    expr = rec_func(act, extract_activities is None)

    if isinstance(expr, float) :
        expr = simplify(expr)
    else:
        # Replace fixed params with their default value
        expr = _replace_fixed_params(expr, _fixed_params().values())

    return (expr, _reverse_dict(act_symbols))


def _reverse_dict(dic):
    return {v: k for k, v in dic.items()}


def exploreImpacts(impact, *activities : ActivityExtended, **params):
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
                new_params = [(name, value) for name, value in _completeParamValues(params).items()]
                amount = amount.subs(new_params)

            i = 1
            ex_name = name
            while ex_name in data:
                ex_name = "%s#%d" % (name, i)
                i += 1

            inputs_by_ex_name[ex_name] = _createTechProxyForBio(input.key, main_act.key[0])

            input_name = _actName(input)

            if _isForeground(input.key[0]) :
                input_name += "{user-db}"

            data[ex_name] = dict(input=input_name, amount=amount)


        # Provide impact calculation if impact provided
        all_acts = list(set(inputs_by_ex_name.values()))
        res = multiLCAAlgebric(all_acts, [impact], **params)
        impacts = res[res.columns.tolist()[0]].to_list()

        impact_by_act = {act : value for act, value in zip(all_acts, impacts)}

        # Add impacts to data
        for key, vals in data.items() :
            amount = vals["amount"]
            act = inputs_by_ex_name[key]
            impact = impact_by_act[act]

            if withImpactPerUnit :
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
        for col_name in diff_cols :
            cols1[col_name] = full.columns.get_loc((names[0], col_name))
            cols2[col_name] = full.columns.get_loc((names[1], col_name))

        if diffOnly:
            def isDiff(row):
                return row[cols1['impact']] != row[cols2['impact']]

            full = full.loc[isDiff]

        def same_amount(row):
            res = [""] * len(row)
            for col_name in diff_cols :
                if row[cols1[col_name]] != row[cols2[col_name]] :
                    res[cols1[col_name]] = YELLOW
                    res[cols2[col_name]] = YELLOW
            return res
        full = full.style.apply(same_amount, axis=1)

    return full