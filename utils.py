#
# This file defines several utility functions above brightway2 to be used by notebooks
#
import re
import sys
from collections import defaultdict
from typing import Dict, Union, List, Any

import brightway2 as bw
import numpy as np
import pandas as pd
from IPython.display import display
from bw2calc.lca import LCA
from bw2data.backends.peewee import Activity, ActivityDataset
from bw2data.parameters import ActivityParameter, ProjectParameter, DatabaseParameter, Group, ExchangeDataset
from sympy import Symbol, Basic, simplify, symbols
from sympy.parsing.sympy_parser import parse_expr
from slugify import slugify
from sympy.utilities.lambdify import lambdify
from collections import OrderedDict
from bw2data.backends.peewee.utils import dict_as_exchangedataset
from copy import deepcopy
from itertools import chain
import builtins

# -- Constants

# DB names
ECOINVENT_DB_NAME='ecoinvent 3.4 cut off'
BIOSPHERE3_DB_NAME='biosphere3'
ACV_DB_NAME='incer-acv'

DEFAULT_PARAM_GROUP="acv"

# Global
def param_registry() :
    if not 'param_registry' in builtins.__dict__ :
        builtins.param_registry = dict()
    return builtins.param_registry

# Sympy symbols
old_amount = symbols("old_amount") # Can be used in epxression of amount for updateExchanges, in order to reference the previous value

# Improved API for activity
class BetterActivity(Activity) :

    # Add method to get exchange by name
    def get_exchange(self, name=None, input=None, single=True) :

        exchs = list([exch for exch in self.exchanges_np() if name and name == exch['name'] or input and input == exch['input']])
        if single and len(exchs) != 1:
            raise Exception("Expected 1 exchange with name '%s' for '%s', found %d" % (name, self, len(exchs)))
        if len(exchs) == 1 :
            return exchs[0]
        elif len(exchs) == 0 :
            return None
        else:
            return exchs
    # Return list of exchanges, except prodution
    def exchanges_np(self) :
        for exch in self.exchanges() :
            if exch['input'] != exch['output'] :
                yield  exch

# Amend the existing class so that all instances of Activity benefit it
setattr(Activity, "get_exchange", BetterActivity.get_exchange)
setattr(Activity, "exchanges_np", BetterActivity.exchanges_np)

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


def isnumber(value) :
    return isinstance(value, int) or isinstance(value, float)

def print_act(*activities, **params):
    '''Print activities and their exchanges.
    If parameter values are provided, formulas will be evaluated accordingly'''
    tables = []
    names = []
    for act in activities :
        df = pd.DataFrame(index = ['input', 'amount', 'unit', 'type'])
        for (i,exc) in enumerate(act.exchanges()):
            input = bw.get_activity(exc.input.key)
            amount = getAmountOrFormula(exc)

            # Params provided ? Evaluate formulas
            if len(params) > 0 and isinstance(amount, Basic) :
                new_params = [(name, value) for name, value in completeParamValues(params).items()]
                amount = amount.subs(new_params)

            name = exc['name']
            if 'location' in input and input['location'] != "GLO" :
                name += "[%s]" % input['location']

            df[name] = [str(input), amount, exc.unit, exc['type']]
        tables.append(df.T)
        names.append(act['name'])

    display(pd.concat(tables, axis = 1, keys = names, sort=True))

# Create or cleanup DB
def resetDb(db_name=ACV_DB_NAME):
    if db_name in bw.databases :
        eprint("Db %s was here. Reseting it" % db_name)
        del bw.databases[db_name]
    db = bw.Database(ACV_DB_NAME)
    db.write(dict())

def importDb(dbname, path) :
    if dbname in bw.databases:
        eprint("Database '%s' has already been imported " % dbname)
    else:
        ei34 = bw.SingleOutputEcospold2Importer(path, dbname)
        ei34.apply_strategies()
        ei34.statistics()
        ei34.write_database()

dbs = dict()
def getDb(dbname)  -> bw.Database :
    '''Pool of Database instances '''
    if not dbname in dbs :
        dbs[dbname] = bw.Database(dbname)
    return dbs[dbname]

def resetParams(db_name=ACV_DB_NAME) :
    ''' Reset project and activity parameters '''
    param_registry().clear()
    ProjectParameter.delete().execute()
    ActivityParameter.delete().execute()
    DatabaseParameter.delete().execute()
    Group.delete().execute()

# Index of activities per name, for fast search dict[db_name][activity_word] => list of activitites
db_index = dict()

def split_words(name) :

    clean = re.sub('[^0-9a-zA-Z]+', ' ', name)
    clean = re.sub(' +', ' ', clean)
    clean = clean.lower()

    return clean.split(' ')


def build_index(db) :
    res = defaultdict(list)
    for act in db :
        words = split_words(act['name'])
        for word in words :
            res[word].append(act)
    return res

def get_indexed_db(db_name) :
    if not db_name in db_index :
        db_index[db_name] = build_index(getDb(db_name))
    return db_index[db_name]

def find_candidates(db_name, name) :
    '''Return the shortest list of candidates among all the words in name'''
    res = []
    index = get_indexed_db(db_name)
    words = split_words(name)
    for word in words :
        candidates = index[word]
        if len(res) == 0 or ( 0 < len(candidates) < len(res)) :
            res = candidates
    return res

def getActByCode(db_name, code) :
    return getDb(db_name).get(code)

def findActivity(name, location=None, categories=None, category=None, db_name=None, single=True) :
    '''
        Find single activity by name.
        Uses index for fast fetching

        :param name: Part of the name to search for
        :param location: [optional] 'GLO' or other
        :param category: if provided, activity should have at least this category
        :param categories: if provided, activity should have this exact list of categories
    '''
    def act_filter(act) :
        if not name == act['name'] :
            return False
        if location and not location == act['location'] :
            return False
        if category and not category in act['categories'] :
            return False
        if categories and not tuple(categories) == act['categories'] :
            return False
        return True

    # Find candidates via index
    candidates = find_candidates(db_name, name)

    # Exact match
    acts = list(filter(act_filter, candidates))
    if single and len(acts) == 0 :
        raise Exception("No activity found in '%s' with name '%s' and location '%s'" % (db_name, name, location))
    if single and len(acts) > 1 :
        raise Exception("Several activity found in '%s' with name '%s' and location '%s':\n%s" % (db_name, name, location, str(acts)))
    if len(acts) == 1 :
        return acts[0]
    else :
        return acts

def findBioAct(name, location=None, **kwargs):
    ''' Alias for findActivity(name, ... db_name=BIOSPHERE3_DB_NAME) '''
    return findActivity(name, location, db_name=BIOSPHERE3_DB_NAME, **kwargs)

def findTechAct(name, location=None, **kwargs):
    ''' Alias for findActivity(name, ... db_name=BIOSPHERE3_DB_NAME) '''
    return findActivity(name, location, db_name=ECOINVENT_DB_NAME, **kwargs)

def interpolate(x, x1, x2, y1, y2):
    return y1 + (y2 - y1) * (x - x1) / (x2 - x1)


def interpolated_act(name:str, act1:BetterActivity, act2:BetterActivity, x1, x2, x, alpha1=1, alpha2=1, **kwargs) :

    res = copyActivity(act1, name, withExchanges=False, **kwargs)

    exch1_by_input = dict({exch['input'] : exch for exch in act1.exchanges_np()})
    exch2_by_input = dict({exch['input']: exch for exch in act2.exchanges_np()})

    inputs = set(chain(exch1_by_input.keys(), exch2_by_input.keys()))

    for input in inputs :

        exch1 = exch1_by_input.get(input)
        exch2 = exch2_by_input.get(input)
        exch = exch1 if exch1 else exch2

        amount1 = exch1['amount'] if exch1 else 0
        amount2 = exch2['amount'] if exch2 else 0

        if exch1 and exch2 and exch1['name'] != exch2['name'] :
            raise Exception("Input %s refer two different names : %s, %s" % (input,  exch1['name'],  exch2['name']))

        amount = interpolate(x, x1, x2, amount1 * alpha1, amount2 * alpha2)
        act = getActByCode(*input)
        addExchanges(res, {act: dict(amount=amount, name=exch['name'])})
    return res

# Type of parameters
class ParamType :
    ENUM = "enum"
    BOOL = "bool"
    FLOAT = "float"

class ParamDef(Symbol) :
    '''
    Generic definition of a parameter, with name, bound, type, distribution
    This definition will serve both to generate brightway2 parameters and to evaluate.

    This class inherits sympy Symbol, making it possible to use in standard arithmetic python
    while keeping it as a symbolic expression (delayed evaluation).
    '''

    def __new__(cls, name, *karg, **kargs):
        return Symbol.__new__(cls, name)

    def __init__(self, name, type : str, default, description=""):

        self.name = name
        self.type = type
        self.default = default
        self.description = description

    # Expand parameter (usefull for enum param)
    def expandParams(self, value=None) -> Dict[str, float] :
        if value == None:
            value = self.default
        return {self.name : value}

    def __repr__(self):
        return self.name

class EnumParam(ParamDef) :
    '''
    Enum param is a facility representing a choice / switch as many 0/1 parameters.
    It is not itself a Sympy symbol. use #symbol("value") to access it
    '''
    def __init__(self, name, values: List[str], **argv):
        super(EnumParam, self).__init__(name, ParamType.ENUM, **argv)

        self.values = values

    def expandParams(self, value=None):
        return {"%s_%s" % (self.name, enum_val) : 1.0 if enum_val == value else 0.0 for enum_val in self.values}

    def symbol(self, enumValue):
        '''
            Access parameter for each enum value :
            <paramName>_<paramValue>
        '''
        if not enumValue in self.values :
            raise Exception("enumValue should be one of %s. Was %s" % (str(self.values), enumValue))
        return Symbol(self.name + '_' + enumValue)



# Creates a param and register it into a global registry and as a brightway parameter
def newParamDef(name, type, **kwargs) :


    if type == ParamType.ENUM :
        param = EnumParam(name, **kwargs)
    else :
        param = ParamDef(name, type=type, **kwargs)

    # Put it in local registry (in memory)
    if name in param_registry():
        eprint("Param %s was already defined : overriding" % name)
    param_registry()[name] = param

    # Save in brightway2 project
    bwParams = [dict(name=key, amount=value) for key, value in param.expandParams().items()]
    bw.parameters.new_project_parameters(bwParams)

    return param

def newFloatParam(name, default, **kwargs) :
    return newParamDef(name, ParamType.FLOAT, default=default, **kwargs)

# Transform amount to either simple amount or formula
def amountToFormula(amount : Union[float, str, Basic], currentAmount=None) :
    res = dict()
    if isinstance(amount, Basic):

        if currentAmount != None :
            amount = amount.subs(old_amount, currentAmount)

        # Check the expression does not reference undefined params
        all_symbols = list([key for param in param_registry().values() for key, val in param.expandParams().items()])
        for symbol in amount.free_symbols :
            if not str(symbol) in all_symbols :
                raise Exception("Symbol '%s' not found in params : %s" % (symbol, all_symbols))

        res['formula'] = str(amount)
        res['amount'] = 0
    elif isinstance(amount, float) or isinstance(amount, int):
        res['amount'] = amount
    else:
        raise Exception(
            "Amount should be either a constant number or a Sympy expression (expression of ParamDef). Was : %s" % type(
                amount))
    return res

NumOrExpression = Union[float, Basic]

def addExchanges(act, exchanges : Dict[Activity, Union[NumOrExpression, dict]] = dict()) :
    ''' Add exchanges to an existing activity, with a compact syntax :
    :param  exchanges Dict of activity => amount or activity => attributes_dict.
            Amount being either a fixed value or Sympy expression (arithmetic expression of Sympy symbols)
    '''
    parametrized = False
    for sub_act, attrs in exchanges.items():

        if isinstance(attrs, dict):
            amount = attrs.pop('amount')
        else :
            attrs = dict()

        exch = act.new_exchange(
            input=sub_act.key,
            name=sub_act['name'],
            unit=sub_act['unit'] if 'unit' in sub_act else None,
            type='biosphere' if sub_act['database'] == BIOSPHERE3_DB_NAME else 'technosphere')

        exch.update(attrs)
        exch.update(amountToFormula(amount))
        if 'formula' in exch :
            parametrized = True

        exch.save()
    act.save()
    if parametrized:
        bw.parameters.add_exchanges_to_group(DEFAULT_PARAM_GROUP, act)

def updateExchanges(activity: Activity, updates: Dict[str, any] = dict()):
    '''
        Update existing exchanges, by name.
        :param updates Dict of exchange name => either single value (float or SympPy expression) for updating only amount,
        or dict of attributes, for updating several at a time. The sympy expression can reference the symbol 'old_amount'
        that will be replaced with the current value.
    '''

    # Update exchanges
    for name, attrs in updates.items():
        exch = activity.get_exchange(name)

        # Single value ? => amount
        if not isinstance(attrs, dict):
            attrs = dict(amount=attrs)

        if 'amount' in attrs:
            attrs.update(amountToFormula(attrs['amount'], exch['amount']))

        exch.update(attrs)
        exch.save()

        # We have a formula now ? => register it to parametrized exchange
        if 'formula' in attrs:
            bw.parameters.add_exchanges_to_group(DEFAULT_PARAM_GROUP, activity)



def getAmountOrFormula(ex:ExchangeDataset) -> Union[Basic, float] :
    '''Return either static amount or Sympy formula'''
    if 'formula' in ex:
        try :
            return parse_expr(ex['formula'])
        except :
            eprint("Error while parsing formula '%s' : backing to amount" % ex['formula'])

    return ex['amount']

def substituteWithDefault(act: Activity, exchange_name: str, switch_act : Activity) :
    '''
        Substitutes one exchange with a switch on other activities, or fallback to the current one as default (parameter set to None)
        For this purpose, we create a new exchange referencing the activity switch, and we multiply current activity by (1-sum(enum_params)),
        making it null as soon as one enum value is set.
        This is useful for changing electricty mix, leaving the default one if needed
    '''
    sum = 0
    for exch in switch_act.exchanges() :
        if exch['input'] != exch['output'] :
            sum += getAmountOrFormula(exch)

    current_exch = act.get_exchange(exchange_name)
    addExchanges(act, {switch_act : getAmountOrFormula(current_exch)})
    updateExchanges(act, {exchange_name : (1-sum)* old_amount})



def _newAct(code, db_name=ACV_DB_NAME) :
    db = getDb(db_name)
    # Already present : delete it ?
    for act in db:
        if act['code'] == code:
            eprint("Activity '%s' was already in '%s'. Overwriting it" % (code, db_name))
            act.delete()

    return db.new_activity(code)

# Creates a new activity
# @param db_name Destination DB : ACV DB by default
# @param exchanges Dict of activity => amount. If amount is a string, is it considered as a formula with parameters
# @param argv : extra params passed as properties of the new activity
def newActivity(
        name, unit=None,
        exchanges : Dict[Activity, Union[float, str]] = dict(),
        db_name=ACV_DB_NAME, **argv) :

    act = _newAct(name, db_name)
    act['name'] = name
    act['unit'] = unit
    act['type'] = 'process'
    act.update(argv)

    # Add exchanges
    addExchanges(act, exchanges)

    return act





def copyActivity(activity: BetterActivity, code=None, db_name=ACV_DB_NAME, withExchanges=True, **kwargs):
    """
        Copy activity into a new DB
    """

    res = _newAct(code, db_name)

    for key, value in activity.items():
        if key not in ['database', 'code'] :
            res[key] = value
    for k, v in kwargs.items():
        res._data[k] = v
    res._data[u'code'] = code
    res['name'] = code
    res.save()

    if withExchanges :
        for exc in activity.exchanges():
            data = deepcopy(exc._data)
            data['output'] = res.key
            # Change `input` for production exchanges
            if exc['input'] == exc['output']:
                data['input'] = res.key
            ExchangeDataset.create(**dict_as_exchangedataset(data))

    return res






# Create a new parametrized, virtual activity, made of a map of other activities, controlled by an enum parameter.
# This enables to implement a "Switch" with brightway parameters
# Internally, this will create a linear sum of other activities controlled by <param_name>_<enum_value> : 0 or 1
# @param paramDef : parameter definition of type enum
# @param acts_dict : dict of <enumValue> => activity
def switch(name, paramDef: ParamDef, acts_dict: Dict[str, Activity]) :

    # Transform map of enum values to correspoding formulas <param_name>_<enum_value>
    exch = {act: paramDef.symbol(key) for key, act in acts_dict.items()}
    res = newActivity(
        name,
        exchanges=exch)

    # Unit of switch activity is the one of the children
    for key, act in acts_dict.items() :
        if 'unit' in act :
            res['unit'] = act['unit']
            res.save()
    return res


# Compte LCA for several methods and several param values : useful for stochastic evaluation
def multi_params_lca(methods, func_unit, nb_iter, update_params) :


    lca = LCA(demand=func_unit, method=methods[0])
    lca.lci(factorize=True)
    method_matrices = []

    results = np.zeros((nb_iter, len(methods)))

    for method in methods:
        lca.switch_method(method)
        method_matrices.append(lca.characterization_matrix)

    for iParam in range(0, nb_iter):
        update_params(iParam)

        lca.redo_lci(func_unit)
        lca.rebuild_biosphere_matrix()
        lca.rebuild_technosphere_matrix()

        for col, cf_matrix in enumerate(method_matrices):
            lca.characterization_matrix = cf_matrix
            lca.lcia_calculation()
            results[iParam, col] = lca.score

    return results

def act_name(act:Activity) :
    '''Generate pretty name for activity, appending location if not "GLO"'''
    res = act['name']
    if act['location'] != 'GLO' :
        res += "[%s]" % act["location"]
    return res

def _multiLCA(activities, methods) :
    ''' Simple wrapper around brightway API'''
    bw.calculation_setups['process'] = {'inv': activities, 'ia': methods}
    lca = bw.MultiLCA('process')
    cols = [act_name(act) for act_amount in activities for act, amount in act_amount.items()]
    return pd.DataFrame(lca.results.T, index=[method[2] for method in methods], columns=cols)

def listOfDictToDictOflist(LD) :
    return {k: [dic[k] for dic in LD] for k in LD[0]}

def completeParamValues(params) :
    '''Check parameters and expand enum params.
    :return Dict of param_name => float value '''

    #undef_params = param_registry.keys() - params.keys()
    #if undef_params :
    #    raise Exception("Some model parameters are not set : %s" % undef_params)

    res = dict()
    for key, val in params.items():
        param = param_registry()[key]

        if isinstance(val, list) :
            newvals = [param.expandParams(val) for val in val]
            res.update(listOfDictToDictOflist(newvals))
        else :
            res.update(param.expandParams(val))
    return res


def multiLCA(model, methods, **params) :

    '''
        Compute LCA for a single activity and a set of methods, after settings the parameters and updating exchange amounts.
    '''

    # Check and expand params
    params = completeParamValues(params)

    # Update brightway parameters
    bwParams = [dict(name=key, amount=value) for key, value in params.items()]
    bw.parameters.new_project_parameters(bwParams)

    # ActivityParameter.recalculate_exchanges(DEFAULT_PARAM_GROUP)
    bw.parameters.recalculate()

    if isinstance(model, list) :
        activities = [{act:1} for act in model]
    else :
        activities = [{model: 1}]
    return _multiLCA(activities, methods).transpose()




def multiLCAAlgebric(models, methods, **params) :
    '''
    Compute LCA by expressing the foreground model as symbolic expression of background activities and parameters.
    Then, compute 'static' inventory of the referenced background activities.
    This enables a very fast recomputation of LCA with different parameters, useful for stochastic evaluation of parametrized model
    :param: Single model or list of models : if list of models, you cannot use param lists
    :param: methods List of methods / impacts to consider
    :param: params You should provide named values of all the parameters declared in the model.
            Values can be single value or list of samples, all of the same size
    '''
    dfs = dict()

    if not isinstance(models, list) :
        models = [models]

    # Check and expand params
    params = completeParamValues(params)

    for model in models :

        expr, actBySymbolName = actToExpression(model)

        # Create dummy reference to biosphere
        # We cannot run LCA to biosphere activities
        # We create a technosphere activity mapping exactly to 1 biosphere item
        pureTechActBySymbol = OrderedDict()
        for name, act in actBySymbolName.items() :
            if act[0] == BIOSPHERE3_DB_NAME :
                act = getOrCreateDummyBiosphereActCopy(act[1])
            else :
                act = getActByCode(*act)
            pureTechActBySymbol[name] = act

        # List of activities, ordered
        acts = pureTechActBySymbol.values()

        # Transform to [{act1:1], {act2:1}, etc] for MultiLCA
        actsWithAmount = [{act:1} for act in acts]

        # Compute LCA for all background activities and methods
        lca = _multiLCA(actsWithAmount, methods)

        # For each method, compute an algebric expression with activities replaced by their values
        lambdas = []
        for imethod, method in enumerate(methods) :

            # Replace activities by their value for this method
            sub = [(symbol, lca.iloc[imethod, iact]) for iact, symbol in enumerate(pureTechActBySymbol.keys())]
            method_expr = expr.subs(sub)

            # Tranform Sympy expression to lambda function, based on numpy to fast vectorial evaluation
            lambd = lambdify(params.keys(), method_expr, 'numpy')
            lambdas.append(lambd)

        # Expand parameters as list of parameters
        param_length = 1

        for key, val in params.items() :
            if isinstance(val, list) :
                if param_length == 1 :
                    param_length = len(val)
                elif param_length != len(val) :
                    raise Exception("Parameters should be a single value or a list of same number of values")

        # Expand params and transform lists to np.array for vector computation
        for key in params.keys() :
            val = params[key]
            if not isinstance(val, list) :
                val = list([val] * param_length)
            params[key] = np.array(val)

        res = np.zeros((len(methods), param_length))


        for imethod, lambd in enumerate(lambdas) :
            # Compute result on whole vectors of parameter samples at a time : lambdas use numpy for vector computation
            res[imethod, :] = lambd(**params)

        name = act_name(act)
        df = pd.DataFrame(res, index=[method[2] for method in methods]).transpose()

        # Single params ? => give the single row the name of the model activity
        if df.shape[0] == 1 :
            df = df.rename(index={0:act_name(model)})

        dfs[name] = df

    if len(dfs) == 1 :
        display(list(dfs.values())[0])
    else :
        # Concat several dataframes for several models
        display(pd.concat(list(dfs.values())))


def getOrCreateDummyBiosphereActCopy(code) :
    '''
        We cannot reference directly biosphere in the model, since LCA can only be applied to products
        We create a dummy activity in our DB, with same code, and single exchange of '1'
    '''
    try:
        return getDb(ACV_DB_NAME).get(code)
    except:
        bioAct = getDb(BIOSPHERE3_DB_NAME).get(code)
        res = newActivity(bioAct['name'] + '#copy', bioAct['unit'], exchanges={bioAct:1})
        return res



def actToExpression(act : Activity):
    '''
        Computes a symbolic expression of the model, referencing background activities and model parameters as symbols
        :return (sympy_expr, dict of symbol => activity)
    '''

    act_symbols = dict() # Dict of  act = > symbol

    def act_to_symbol(db_name, code) :
        act = getDb(db_name).get(code)
        name = act['name']
        base_slug = slugify(name, separator='_')

        slug = base_slug
        i=1
        while symbols(slug) in act_symbols.values():
            slug = f"{base_slug}{i}"
            i += 1

        return symbols(slug)

    def rec_func(act:Activity) :
        res = 0

        for exch in act.exchanges():

            formula = getAmountOrFormula(exch)

            input_db, input_code = exch['input']

            #  Ignore output
            if exch['input'] == exch['output'] :
                if exch['amount'] != 1 :
                    raise Exception("Output quantity not 1")
                continue

            # Background DB => reference it as a symbol
            if input_db != ACV_DB_NAME :
                if not (input_db, input_code) in act_symbols :
                    act_symbols[(input_db, input_code)] = act_to_symbol(input_db, input_code)
                act_expr = act_symbols[(input_db, input_code)]

            # Our model : recursively transform it to a symbolic expression
            else :

                if input_db == act['database'] and input_code == act['code']:
                    raise Exception("Recursive exchange : %s" % (act.__dict__))

                sub_act = getDb(input_db).get(input_code)
                act_expr = rec_func(sub_act)

            res += formula * act_expr

        return res

    expr =  rec_func(act)

    return (expr, reverse_dict(act_symbols))

def reverse_dict(dic) :
    return {v: k for k, v in dic.items()}