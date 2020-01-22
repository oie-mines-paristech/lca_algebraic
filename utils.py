#
# This file defines several utility functions above brightway2 to be used by notebooks
#
import re
import sys
from collections import defaultdict
from typing import Dict, Union, List

import brightway2 as bw
import numpy as np
import pandas as pd
from IPython.display import display
from bw2calc.lca import LCA
from bw2data.backends.peewee import Activity, ActivityDataset
from bw2data.parameters import ActivityParameter, ProjectParameter, DatabaseParameter, Group
from sympy import Symbol, Basic, simplify, symbols
from sympy.parsing.sympy_parser import parse_expr
from slugify import slugify
from sympy.utilities.lambdify import lambdify
from collections import OrderedDict


# Constants
ECOINVENT_DB_NAME='ecoinvent 3.4 cut off'
BIOSPHERE3_DB_NAME='biosphere3'
ACV_DB_NAME='incer-acv'
DEFAULT_PARAM_GROUP="acv"

# Global
param_registry = dict()

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

# Print activity
def print_act(act):
    df = pd.DataFrame(index = ['amount /formula','unit','type'])
    for (i,exc) in enumerate(act.exchanges()):
        df[str(bw.get_activity(exc.input.key))] = [exc['formula'] if 'formula' in exc else exc.amount, exc.unit,exc['type']]
    display(df.T)

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
    global param_registry
    param_registry = dict()
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

def findActivity(name, location=None, categories=None, category=None, db_name=None) :
    '''
        Find single activity by name.
        Uses index for fast fetching

        :param name: Part of the name to search for
        :param location: [optional] 'GLO' or other
        :param category: if provided, activity should have at least this category
        :param categories: if provided, activity should have this exact list of categories
    '''
    def act_filter(act) :
        if not name in act['name'] :
            return False
        if location and not location in act['location'] :
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
    if len(acts) == 0 :
        raise Exception("No activity found in '%s' with name '%s' and location '%s'" % (db_name, name, location))
    if len(acts) > 1 :
        raise Exception("Several activity found in '%s' with name '%s' and location '%s':\n%s" % (db_name, name, location, str(acts)))
    return acts[0]

def findBioAct(name, location=None, **kwargs):
    ''' Alias for findActivity(name, ... db_name=BIOSPHERE3_DB_NAME) '''
    return findActivity(name, location, db_name=BIOSPHERE3_DB_NAME, **kwargs)

def findTechAct(name, location=None, **kwargs):
    ''' Alias for findActivity(name, ... db_name=BIOSPHERE3_DB_NAME) '''
    return findActivity(name, location, db_name=ECOINVENT_DB_NAME, **kwargs)

# Type of parameters
class ParamType :
    ENUM = "enum"
    BOOL = "bool"
    NUMBER = "number"

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

    # Transform to bw2 params (several params for enum type, one for others)
    def toBwParams(self, value=None):
        if value == None:
            value = self.default
        return [dict(
            name=self.name,
            amount=value)]

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

    def toBwParams(self, value=None):
        res = []
        for enum_val in self.values:
            res.append(dict(
                name="%s_%s" % (self.name, enum_val),
                amount=1.0 if enum_val == value else 0.0))
        return res

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
    global param_registry

    if type == ParamType.ENUM :
        param = EnumParam(name, **kwargs)
    else :
        param = ParamDef(name, type=type, **kwargs)

    # Put it in local registry (in memory)
    if name in param_registry:
        raise Exception("Param %s already defined" % name)
    param_registry[name] = param

    # Save in brightway2 project
    bw.parameters.new_project_parameters(param.toBwParams())

    return param

# Creates a new activity
# @param db_name Destination DB : ACV DB by default
# @param exchanges Dict of activity => amount. If amount is a string, is it considered as a formula with parameters
# @param argv : extra params passed as properties of the new activity
def newActivity(
        name, unit=None,
        exchanges : Dict[Activity, Union[float, str]] = dict(),
        db_name=ACV_DB_NAME, **argv) :

    db = getDb(db_name)

    # Already present : delete it ?
    for act in db :
        if act['code'] == name :
            eprint("Activity '%s' was already in '%s'. Overwriting it" % (name, db_name))
            act.delete()

    #code = binascii.hexlify(os.urandom(16))
    act = db.new_activity(name)
    act['name'] = name
    act['unit'] = unit
    act['type'] = 'process'
    act.update(argv)

    # Add exchanges
    parametrized = False
    for sub_act, amount in exchanges.items() :
        exch = act.new_exchange(
            input=sub_act.key,
            name=sub_act['name'],
            unit=sub_act['unit'] if unit in sub_act else None,
            type='biosphere' if sub_act['database'] == BIOSPHERE3_DB_NAME else 'technosphere')
        if isinstance(amount, Basic) :
            parametrized = True
            exch.expr = amount
            exch['formula'] = str(amount)
            exch['amount'] = 0
        elif isinstance(amount, float) or isinstance(amount, int) :
            exch['amount'] = amount
        else :
            raise Exception("Amount should be either a constant number or a Sympy expression (expression of ParamDef). Was : %s" % type(amount))
        exch.save()

    act.save()
    if parametrized :
       bw.parameters.add_exchanges_to_group(DEFAULT_PARAM_GROUP, act)
    return act

# Create a new parametrized, virtual activity, made of a map of other activities, controlled by an enum parameter.
# This enables to implement a "Switch" with brightway parameters
# Internally, this will create a linear sum of other activities controlled by <param_name>_<enum_value> : 0 or 1
# @param paramDef : parameter definition of type enum
# @param acts_dict : dict of <enumValue> => activity
def switch(paramDef: ParamDef, acts_dict: Dict[str, Activity]) :

    # Transform map of enum values to correspoding formulas <param_name>_<enum_value>
    exch = {act: paramDef.symbol(key) for key, act in acts_dict.items()}
    res = newActivity(
        paramDef.name + " switch" ,
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


def multiLCA(activities, methods) :
    ''' Simple wrapper around brightway API'''
    bw.calculation_setups['process'] = {'inv': activities, 'ia': methods}
    lca = bw.MultiLCA('process')
    cols = [act['code'] for act_amount in activities for act, amount in act_amount.items()]
    return pd.DataFrame(lca.results.T, index=methods, columns=cols)

def multiLCAWithParams(activity, methods, amount=1, **params) :

    '''
        Compute LCA for a single activity and a set of methods, after settings the parameters and updating exchange amounts.
    '''

    # Update parameters in brigthway
    if params:
        bw_params = []
        for key, val in params.items():
            param = param_registry[key]
            bw_params.extend(param.toBwParams(val))

        bw.parameters.new_project_parameters(bw_params)
        # ActivityParameter.recalculate_exchanges(DEFAULT_PARAM_GROUP)
        bw.parameters.recalculate()

    #undef_params = param_registry.keys() - params.keys()
    #if undef_params :
    #    raise Exception("Some model parameters are not set : %s" % undef_params)

    return multiLCA([{activity: amount}], methods)


def multiLCAWithParamsAlgebric(activity, methods, amount=1, **params) :

    '''
        Compute LCA by expressing the foregrounf model as symbolic expression of background
        activities and parameters. Then, compute 'static' inventory of the referenced background activities.
        This enables a very fast recomputation of LCA with different parameters, usefull for stochastic evaluation of parametrized models
    '''

    expr, actBySymbolName = actToExpression(activity)

    # Create dummy reference to biosphere
    # We cannot run LCA to biosphere
    # We create technosphere activity mappaing exactly to 1 biosphere item
    actBySymbolNameNoBio = OrderedDict()
    for name, act in actBySymbolName.items() :
        if act[0] == BIOSPHERE3_DB_NAME :
            act = getOrCreateDummyBiosphereActCopy(act[1])
        else :
            act = getActByCode(*act)
        actBySymbolNameNoBio[name] = act

    # List of activities, ordered
    acts = actBySymbolNameNoBio.values()

    # Transform to [{act1:1], {act2:1}, etc] for MultiLCA
    actsWithAmount = [{act:1} for act in acts]

    # Compute LCA for all background activities and methods
    lca =  multiLCA(actsWithAmount, methods)

    # For each method, compute an algebric expression with activities replaced by their values
    lambdas = []
    for imethod, method in enumerate(methods) :

        # Replace activities by their value for this method
        sub = [(symbol, lca.iloc[imethod, iact]) for iact, symbol in enumerate(actBySymbolNameNoBio.keys())]
        method_expr = expr.subs(sub)

        # Tranform Sympy expression to lambda function to fast evaluation
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
        # Compute result on whole vertors at a time : lambdas use numpy for vetor computation
        res[imethod, :] = lambd(**params)

    return pd.DataFrame(res, index=[method[2] for method in methods]).transpose()



def getOrCreateDummyBiosphereActCopy(code) :
    '''
        We cannot reference directly biosphere in the model, since LCA can only be applied to products
        We create a dummy activity in our DB, with same code, and single exchange of '1'
    '''
    try:
        return getDb(ACV_DB_NAME).get(code)
    except:
        bioAct = getDb(BIOSPHERE3_DB_NAME).get(code)
        res = newActivity(bioAct['name'] + 'copy', bioAct['unit'], exchanges={bioAct:1})
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
        i=0
        while slug in act_symbols.values():
            slug =  f"{base_slug}{i}"

        return symbols(slug)

    def rec_func(act:Activity) :

        res = 0

        for exch in act.exchanges() :

            if 'formula' in exch :
                amount = parse_expr(exch['formula'], param_registry)
            else :
                amount = exch['amount']

            db_name, act_key = exch['input']

            # Background DB => reference it as a symbol
            if db_name != ACV_DB_NAME :
                if not (db_name, act_key) in act_symbols :
                    act_symbols[(db_name, act_key)] = act_to_symbol(db_name, act_key)
                act_expr = act_symbols[(db_name, act_key)]

            # Our model : recursively transform it to a symbolic expression
            else :
                if db_name == act['database'] and act_key == act['code']:
                    raise Exception("Recursive exchange : %s" % (act.__dict__))

                sub_act = getDb(db_name).get(act_key)
                act_expr = rec_func(sub_act)

            res += amount * act_expr

        return simplify(res)

    expr =  rec_func(act)

    return (expr, reverse_dict(act_symbols))

def reverse_dict(dic) :
    return {v: k for k, v in dic.items()}