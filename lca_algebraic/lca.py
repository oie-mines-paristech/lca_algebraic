import concurrent.futures
from collections import OrderedDict
from functools import partial
from typing import Dict, List

from IPython import display
from sympy import lambdify, simplify

from .base_utils import _actName, error, _getDb, _method_unit
from .base_utils import _getAmountOrFormula
from .helpers import *
from .helpers import _actDesc, _isForeground
from .params import _param_registry, _completeParamValues, _fixed_params, _expanded_names_to_names, _expand_param_names
from ipywidgets import HTML, Label, Button, HBox, VBox, Layout, Combobox, Output, Accordion, ToggleButton
from ipyevents import Event


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

            debug("computing : ", model)
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


def _uniq_name(name, names) :
    i = 1
    res = name
    while res in names:
        res = "%s#%d" % (name, i)
        i += 1
    return res

def _ellipsis(value, max=40) :
    if len(value) < max :
        return value
    else :
        short = value[0:max]
        return "<span title='%s'>%s…</span>" % (value, short)

def _sci_repr(float_val) :
    if isinstance(float_val, float):
        return "%02.3g" % float_val
    else:
        return float_val


def _none_float(float_val) :
    return 0.0 if float_val is None else float_val

_SMALL_DIFF_COLOR = "LemonChiffon"
_HIGH_DIFF_COLOR = "Coral"
_RIGHT_BORDER = {'border-right-width':'2px', 'border-right-color':'black'}
_DB_COLORS = ["lightgreen", "lightyellow", "lightsteelblue", "lightpink", "NavajoWhite", "khaki"]

def _db_shortname(db_name) :
    """Return 'shortname' metadata if Db has one. Otherwize, return DB name"""
    db = bw.Database(db_name)
    if "shortname" in db.metadata :
        return db.metadata["shortname"]
    return db_name

def _db_idx(db) :
    """Return index if DB, in list of all DBs"""
    dbs = list(bw.databases.keys())
    return dbs.index(db)

def _explore_impacts(
        activities : ActivityExtended,
        method,
        withImpactPerUnit=False,
        diff_impact=0.1,
        act_callback=None,
        **params):
    """
    Advanced version of #printAct()

    Displays all exchanges of one or several activities and their impacts.
    If parameter values are provided, formulas will be evaluated accordingly.
    If two activities are provided, they will be shown side by side and compared.
    """

    if not isinstance(activities, (list, tuple)):
        activities = [activities]

    # List of dict[exchange=>attributes]
    # As many items as source activities
    datas = list()
    nact = len(activities)

    for main_act in activities:

        inputs_by_ex_name = dict()

        # Dict of exchange name => attributes
        data = dict()
        datas.append(data)

        for ex in main_act.listExchanges():

            amount = ex.amount

            # Params provided ? Evaluate formulas
            if len(params) > 0 and isinstance(ex.amount, Basic):
                new_params = [(name, value) for name, value in _completeParamValues(params).items()]
                amount = ex.amount.subs(new_params)

            ex_name = _uniq_name(ex.name, data)
            inputs_by_ex_name[ex_name] = _createTechProxyForBio(ex.input.key, main_act.key[0])

            data[ex_name] = SimpleNamespace(
                input=ex.input["name"],
                amount=amount,
                db= ex.input.key[0],
                code= ex.input.key[1],
                loc="GLO" if not "location" in ex.input else ex.input["location"],
                unit=ex.unit)

        # Provide impact calculation if impact provided
        all_acts = list(set(inputs_by_ex_name.values()))

        res = multiLCAAlgebric(all_acts, [method], **params)
        impacts = res[res.columns.tolist()[0]].to_list()

        impact_by_act = {act : value for act, value in zip(all_acts, impacts)}

        # Add impacts to data
        for key, attrs in data.items() :
            amount = attrs.amount
            act = inputs_by_ex_name[key]
            impact = impact_by_act[act]

            if withImpactPerUnit :
                attrs.impact_per_unit = impact
            attrs.impact = amount * impact

    # All exch
    all_ex_names = list(set(ex for data in datas for ex in data.keys()))

    headers = ["exchange", "unit"]

    # Exchange, Unit, Inputs, Amounts, Impacts
    widths = [170, 50] \
             + [210] * nact \
             + [70] * nact \
             + [70] * nact

    def _numerate(str_val) :
        return [str_val] + ["%s #%d" % (str_val, i) for i in range(1, nact)]

    headers += _numerate("input") + _numerate("amount") + _numerate("impact")

    # Compute sum of impacts for first activity
    sum_impact = sum(attrs.impact for attrs in datas[0].values())

    # To sheet
    sheet  = ipysheet.sheet(
        rows=len(all_ex_names),
        columns=len(headers),
        column_headers=headers,
        row_headers=False,
        stretch_headers="none",
        column_width=widths)

    # Loop on exchanges
    for i_ex, ex_name in enumerate(all_ex_names) :

        ipysheet.cell(row=i_ex, column=0, value=HTML("<b class='trunc' title='%s'>%s</b>" % (ex_name, ex_name)))

        # Loop on type of columns
        for column in ["input", "amount", "impact"] :

            # On or two cells, one per activity
            values = []
            cells = []

            # Loop on activity (1 or 2)
            for i_act, (act, data) in enumerate(zip(activities, datas)) :

                # Empty attr for this exchange ?
                attrs = data[ex_name] if ex_name in data else SimpleNamespace(
                    unit="-",
                    input="_",
                    amount="_",
                    impact=None,
                    db="-",
                    loc="-")

                # Override unit (common column)
                if attrs.unit != "-" :
                    ipysheet.cell(row=i_ex, column=1, value=attrs.unit, style=_RIGHT_BORDER, read_only=True)

                # Switch on column type
                if column == "input" :

                    if attrs.input == "_":
                        input = "-"
                        values.append(input)
                    else:
                        input_name = "<a href='#'>%s</a>" % attrs.input
                        loc_html = "<span class='loc'> %s </span>" % attrs.loc
                        db_html = "&nbsp;<span class='db%d'> %s </span>" % (_db_idx(attrs.db), _db_shortname(attrs.db))
                        html = input_name + "<br/>" + loc_html + db_html

                        values.append(html)

                        input = HTML(html)

                        # Link the callback to click on it
                        if act_callback is not None :

                            sub_act = getActByCode(attrs.db, attrs.code)

                            def cb(sub_act, event):
                                act_callback(sub_act)

                            event = Event(source=input, watched_events=['click'])
                            event.on_dom_event(partial(cb, sub_act))

                    cells.append(input)

                elif column == "amount":

                    values.append(attrs.amount)
                    cells.append(_sci_repr(attrs.amount))

                else:

                    values.append(attrs.impact)
                    cells.append(_sci_repr(attrs.impact))

            # Colorize background according to Diff
            color = None

            if (nact > 1) and values[0] != values[1] :

                # For impact, highlight differently difference higher than a given share of total impact
                if column == "impact":
                    diff = abs(
                        _none_float(values[1]) -
                        _none_float(values[0]))
                    rel_diff = diff / sum_impact

                    color = _HIGH_DIFF_COLOR if rel_diff > diff_impact else _SMALL_DIFF_COLOR

                else:
                    color = _SMALL_DIFF_COLOR


            # Display cells for this column
            for icell, cell in enumerate(cells) :
                icol = 2 + nact * (0 if column == "input" else 1 if column == "amount" else 2)


                style= _RIGHT_BORDER if (icell == nact - 1) else None

                ipysheet.cell(
                    row=i_ex,
                    column=icol + icell,
                    style=style.copy() if style else None,
                    background_color=color,
                    value=cell,
                    read_only=True)

    # Styling

    style = """
    .trunc {
           white-space: nowrap;
           overflow: hidden;
           text-overflow: ellipsis;
    }"""

    # Add colors for each DB number
    style += "".join(""".db%d {
            background-color : %s    
        }""" % (i, color) for i, color in enumerate(_DB_COLORS))

    style = HTML("""<style>%s</style>""" % style)
    display(style)

    return sheet



def _acts_by_name(db_name) :
    return {_actName(act):act for act in bw.Database(db_name)}

def explore_impacts(main_db, base_db, method):
    """
    Interactive interface to explore activities, compare them to other activities of other DB and compare impacts

    :param main_db: Name of main DB
    :param base_db: Name of base DB
    :param method: Method for impact calculation
    :return: Interactive widget ready to be displayed
    """
    message = HTML()
    table_container = VBox()
    history = []

    home_button = Button(icon="home")
    back_button = Button(icon="arrow-left")

    output = Output()

    main_act = None
    base_act = None

    # List activities by name in other DB
    main_acts_by_name = _acts_by_name(main_db)
    base_acts_by_name = _acts_by_name(base_db)

    main_combo = Combobox(
        description=_db_shortname(main_db),
        layout=widgets.Layout(width='700px'),
        options=list(main_acts_by_name.keys()),
        placeholder="Select activity",
        ensure_option=True)

    base_combo = Combobox(
        description=_db_shortname(base_db),
        layout=widgets.Layout(width='700px'),
        options=list(base_acts_by_name.keys()),
        placeholder="Select activity",
        ensure_options=True)

    # Comments
    main_comment = HTML(layout=Layout(display="none"))
    base_comment = HTML(layout=Layout(display="none"))

    main_button_switch = ToggleButton(
        value=False,
        icon="chevron-down",
        description="Comment")

    base_button_switch = ToggleButton(
        value=False,
        icon="chevron-down",
        description="Comment")

    with output :

        def update_main(act) :
            """Update main activity. Try to update base activity with same name"""

            nonlocal main_act, base_act

            if main_act == act :
                return

            if main_act is not None or base_act is not None :
                history.append((main_act, base_act))

            main_act = act


            name = _actName(act)

            if act.key[0] == main_db and name in base_acts_by_name:
                base_act = base_acts_by_name[name]
            else:
                base_act = None


            update_view()

        def update_base(act) :

            nonlocal base_act

            if base_act == act:
                return

            base_act = act

            history.append((main_act, base_act))
            update_view()

        def pop_and_update(event) :
            nonlocal main_act, base_act
            if len(history) > 0 :
                main_act, base_act = history.pop()
                update_view()

        def reset(event) :
            nonlocal main_act, base_act
            if len(history) > 0 :
                main_act, base_act = history[0]
                history.clear()
                update_view()

        def update_view():

            # Update names in combo
            main_combo.value = "" if main_act is None else _actName(main_act)
            base_combo.value = "" if base_act is None else _actName(base_act)

            # Update comments
            main_comment.value = "" if main_act is None else main_act["comment"]
            base_comment.value = "" if base_act is None else base_act["comment"]

            message.value = "<h4>Computing ...</h4>"

            try:
                acts = main_act if base_act is None else [main_act, base_act]
                table = _explore_impacts(acts, method, act_callback=update_main)
                table_container.children = [table]

            finally:
                message.value = ""

        def main_combo_change(change):
            val = change.new
            if not val in main_acts_by_name :
                return
            act = main_acts_by_name[val]
            update_main(act)

        def base_combo_change(change):
            val = change.new
            if not val in base_acts_by_name :
                return
            act = base_acts_by_name[val]
            update_base(act)

        def comment_visibility_switch(container, change) :
            val = change.new
            container.layout.display = "block" if val else "none"


    # Bind events
    home_button.on_click(reset)
    back_button.on_click(pop_and_update)

    main_combo.observe(main_combo_change, "value")
    base_combo.observe(base_combo_change, "value")

    main_button_switch.observe(partial(comment_visibility_switch, main_comment), "value")
    base_button_switch.observe(partial(comment_visibility_switch, base_comment), "value")

    display(VBox([
        HBox([home_button, back_button]),
        HBox([main_combo, main_button_switch]),
        main_comment,
        HBox([base_combo, base_button_switch]),
        base_comment,
        message,
        table_container,
        output]))
