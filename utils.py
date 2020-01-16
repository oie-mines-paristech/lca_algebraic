#
# This file defines several utility functions above brightway2 to be used by notebooks
#
import brightway2 as bw
from bw2data.backends.peewee import Activity
from bw2data.parameters import ActivityParameter
from typing import Dict, Union
import pandas as pd
import sys

# Constants
ECOINVENT_DB_NAME='ecoinvent 3.4 cut off'
BIOSPHERE3_DB_NAME='biosphere3'
ACV_DB_NAME='incer-acv'

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

# Print activity
def print_act(act):
    df = pd.DataFrame(index = ['amount /formula','unit','type'])
    for (i,exc) in enumerate(act.exchanges()):
        df[str(bw.get_activity(exc.input.key))] = [exc['formula'] if 'formula' in exc else exc.amount, exc.unit,exc['type']]
    return df.T

# Create or cleanup DB
def initDb(db_name=ACV_DB_NAME) :
    if db_name in bw.databases :
        del bw.databases[db_name]
    db = bw.Database(ACV_DB_NAME)
    db.write(dict())

# Find single activity by name and location
def findActivity(name, location=None, db_name=ECOINVENT_DB_NAME) :
    def act_filter(act) :
        if not name in act['name'] :
            return False
        if location and not location in act['location'] :
            return False
        return True
    db = bw.Database(db_name)
    acts = list(filter(act_filter, db))
    if len(acts) == 0 :
        raise Exception("No activity found in '%s' with name '%s' and location '%s'" % (db_name, name, location))
    if len(acts) > 1 :
        raise Exception("Several activity found in '%s' with name '%s' and location '%s':\n%s" % (db_name, name, location, str(acts)))
    return acts[0]

# Type of parameters
class ParamType :
    ENUM = "enum"
    BOOL = "bool"
    NUMBER = "number"

# Generic definition of a parameter, with name, bound, type, distribution
# This definition will serve both to generate brightway2 parameters and to evaluate
class ParamDef :
    def __init__(self, name, type : ParamType, default, description="", **argv):
        self.name = name
        self.type = type
        self.default = default
        self.description = description
        if type == ParamType.ENUM :
            if not 'values' in argv :
                raise Exception("You should provide a list of possible values for enum")
            self.values = argv['values']

    # Transform to bw2 params (several params for enum type, one for others)
    def toBwParams(self, value=None):
        if value == None :
            value = self.default
        if self.type == ParamType.ENUM :
            res = []
            for enum_val in self.values :
                res.append(dict(
                    name="%s_%s" % (self.name, enum_val),
                    amount = 1.0 if enum_val == value else 0.0))
            return res
        else :
            return [dict(
                name=self.name,
                amount=value)]

# Creates a param and register it into global registry and as bright way parameter
param_registry = dict()
def newParamDef(name, *args, **kwargs) :
    param = ParamDef(name, *args, **kwargs)
    if name in param_registry:
        raise Exception("Param %s already defined" % name)
    param_registry[name] = param

    # Save in brightway2 project
    bw.parameters.new_project_parameters(param.toBwParams())

    return param

# Creates a new activity
# @param db_name Destination DB : ACER DB by default
# @param exchanges Dict of activity => amount. If amount is a string, is it considered as a formula with parameters
def newActivity(
        name, unit=None,
        exchanges : Dict[Activity, Union[float, str]] = None,
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
    act.update(argv)

    # Add exchanges
    parametrized = False
    if exchanges :
        for sub_act, amount in exchanges.items() :
            ex = act.new_exchange(
                input=sub_act.key,
                name=sub_act['name'],
                unit=sub_act['unit'] if unit in sub_act else None,
                type='technosphere')
            if isinstance(amount, str) :
                parametrized = True
                ex['formula'] = amount
                ex['amount'] = 0
            else :
                ex['amount'] = amount
            ex.save()
    act.save()
    if parametrized :
        bw.parameters.add_exchanges_to_group(db_name, act)
    return act

# Create a new paramatrized, virtual activity, made of a map of other activities, controlled by an enum parameter.
# This enables to implement a "Switch" with brightway parameters
# Internally, this will create a linear sum of other activities controlled by <param_name>_<enum_value> : 0 or 1
# @param paramDef : parameter definition of type enum
# @param map : dict of <enumValue> => activity
def newSwitchActivity(paramDef: ParamDef, map: Dict[str, Activity]) :
    res = newActivity(paramDef.name + " switch")

    for key, act in map.items() :
        if not key in paramDef.values :
            raise Exception("Value '%s' not defined in '%s'" % (key, paramDef.name))

        exchange = res.new_exchange(
            amount = 0,
            type='technosphere',
            input = act.key,
            name = act['name'],
            formula='%s_%s' % (paramDef.name, key))
        if 'unit' in act :
            exchange['unit'] = act['unit']

            # Switch activity has same unit as child ones
            res['unit'] = act['unit']
        exchange.save()

    res.save()
    return res

# Compute LCA for a single activity and a set of methods, after settings the parameters
def computeLCA(activity, methods, amount=1, **params) :

    # Update parameters in brigthway
    if params:
        bw_params= []
        for key, val in params.items() :
            param = param_registry[key]
            bw_params.extend(param.toBwParams(val))
        print(bw_params)
        bw.parameters.new_project_parameters(bw_params)
        bw.parameters.recalculate()

    bw.calculation_setups['process'] = {'inv': [{activity:amount}], 'ia': methods}
    lca = bw.MultiLCA('process')
    df = pd.DataFrame(index=methods)
    df[0] = lca.results.T
    return df