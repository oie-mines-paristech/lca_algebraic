#
# This file defines several utility functions above brightway2 to be used by notebooks
#
import brightway2 as bw
from bw2data.backends.peewee import Activity
from bw2data.parameters import ActivityParameter, ProjectParameter, DatabaseParameter, Group
from typing import Dict, Union, List
from bw2calc.lca import LCA
import pandas as pd
import sys
from sympy import Symbol, Basic
import numpy as np

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
    return df.T

# Create or cleanup DB
def resetDb(db_name=ACV_DB_NAME):
    if db_name in bw.databases :
        print("Db %s was here. Deleting it" % db_name)
        del bw.databases[db_name]
    db = bw.Database(ACV_DB_NAME)
    db.write(dict())


def resetParams(db_name=ACV_DB_NAME) :
    ''' Reset project and activity parameters '''
    global param_registry
    param_registry = dict()
    ProjectParameter.delete().execute()
    ActivityParameter.delete().execute()
    DatabaseParameter.delete().execute()
    Group.delete().execute()


def findActivity(name, location=None, categories=None, category=None, db_name=ECOINVENT_DB_NAME) :
    '''
        Find single activity by name

        :param name: Part of the name to search for
        :param location: [optional] 'GLO' or other
        :param category: if provided, activity should have at least this category
        :param categories: if provided, activities should be exactly this
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
    db = bw.Database(db_name)
    acts = list(filter(act_filter, db))
    if len(acts) == 0 :
        raise Exception("No activity found in '%s' with name '%s' and location '%s'" % (db_name, name, location))
    if len(acts) > 1 :
        raise Exception("Several activity found in '%s' with name '%s' and location '%s':\n%s" % (db_name, name, location, str(acts)))
    return acts[0]

def findBiosphere(name, location=None, **kwargs):
    ''' Alias for findActivity(name, ... db_name=BIOSPHERE3_DB_NAME) '''
    return findActivity(name, location, db_name=BIOSPHERE3_DB_NAME, **kwargs)


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



# Creates a param and register it into global registry and as brightway parameter
def newParamDef(name, type, **kwargs) :
    global param_registry

    if type == ParamType.ENUM :
        param = EnumParam(name, **kwargs)
    else :
        param = ParamDef(name, type=type, **kwargs)

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

    db = bw.Database(db_name)

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

        print(lca.tech_params)

        lca.redo_lci(func_unit)
        lca.rebuild_biosphere_matrix()
        lca.rebuild_technosphere_matrix()

        for col, cf_matrix in enumerate(method_matrices):
            lca.characterization_matrix = cf_matrix
            lca.lcia_calculation()
            results[iParam, col] = lca.score

    return results


def computeLCA(activity, methods, amount=1, **params) :

    '''
        Compute LCA for a single activity and a set of methods, after settings the parameters and updateing exchange amounts
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

    bw.calculation_setups['process'] = {'inv': [{activity: amount}], 'ia': methods}
    lca = bw.MultiLCA('process')
    df = pd.DataFrame(index=methods)
    df[0] = lca.results.T
    return df
