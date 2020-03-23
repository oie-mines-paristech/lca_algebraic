#
# This file defines several utility functions above brightway2 to be used by notebooks
#
import re
import sys
from collections import defaultdict
from typing import Dict, Union, List, Any
import types
import brightway2 as bw
import numpy as np
import pandas as pd
from IPython.display import display
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
import matplotlib.pyplot as plt
from ipywidgets import interact, interactive, fixed, interact_manual
import ipywidgets as widgets
import warnings
from SALib.sample import saltelli
import math
from SALib.analyze import sobol
from scipy.stats import binned_statistic
import seaborn as sns
from sys import stderr
import math
from enum import Enum
from scipy.stats import triang
from scipy.stats import truncnorm
from math import sqrt

# -- Constants
DEBUG=False

def debug(*args, **kwargs) :
    if DEBUG :
        print(*args, **kwargs)

def error(*args, **kwargs):
    print(*args, **kwargs, file=stderr)

# DB names

# Set by detect_ecoinvent
ECOINVENT_DB_NAME = None
BIOSPHERE3_DB_NAME = 'biosphere3'

USER_DB_NAME = None

DEFAULT_PARAM_GROUP = "acv"

# Global
def _param_registry() :

    # Prevent reset upon auto reload in jupyter notebook
    if not 'param_registry' in builtins.__dict__:
        builtins.param_registry = dict()

    return builtins.param_registry


# Sympy symbols
old_amount = symbols("old_amount")  # Can be used in epxression of amount for updateExchanges, in order to reference the previous value

NumOrExpression = Union[float, Basic]



class ParamType:
    '''Type of parameters'''
    ENUM = "enum"
    BOOL = "bool"
    FLOAT = "float"

class DistributionType :
    '''Type of distribution'''
    LINEAR = "linear"
    NORMAL = "normal"
    TRIANGLE = "triangle"
    FIXED = "fixed"

class ParamDef(Symbol):
    '''Generic definition of a parameter, with name, bound, type, distribution
    This definition will serve both to generate brightway2 parameters and to evaluate.

    This class inherits sympy Symbol, making it possible to use in standard arithmetic python
    while keeping it as a symbolic expression (delayed evaluation).
    '''

    def __new__(cls, name, *karg, **kargs):
        return Symbol.__new__(cls, name)

    def __init__(self, name, type: str, default, min=None, max=None, unit="", description="", label=None, label_fr=None, group=None, distrib=DistributionType.LINEAR, std=None):
        self.name = name
        self.type = type
        self.default = default
        self.description = description
        self.min = min
        self.max = max
        self.unit = unit
        self.label = label
        self.label_fr = label_fr
        self.group=group
        self.distrib = distrib

        if type == ParamType.FLOAT and self.min is None :
            self.distrib = DistributionType.FIXED

        if distrib == DistributionType.NORMAL and std is None :
            raise Exception("Standard deviation is mandatory for normal distribution")
        self.std = std

    def label(self):
        if self.label is not None :
            return self.label
        else :
            return self.name.replace("_", " ")

    def range(self, n) :
        '''Used for parametric analysis'''
        step = (self.max - self.min) / (n - 1)
        return list(i * step + self.min for i in range(0, n))


    def rand(self, alpha):
        """Transforms a random number between 0 and 1 to valid value according to the distribution of probability of the parameter"""
        if self.distrib == DistributionType.LINEAR :
            return self.min + alpha * (self.max - self.min)

        elif self.distrib == DistributionType.TRIANGLE :
            if not hasattr(self, "_distrib") :
                scale = self.max - self.min
                c = (self.default - self.min) / scale
                self._distrib = triang(c, loc=self.min, scale=scale)

            return self._distrib.ppf(alpha)

        elif self.distrib == DistributionType.NORMAL :
            if not hasattr(self, "_distrib") :
                self._distrib = truncnorm(
                    (self.min - self.default) / self.std,
                    (self.max - self.min) / self.std,
                    loc=self.default,
                    scale=self.std)

            return self._distrib.ppf(alpha)

        else :
            raise Exception("Unknowk distribution type " + self.distrib)



    # Expand parameter (useful for enum param)
    def expandParams(self, value=None) -> Dict[str, float]:
        if value == None:
            value = self.default
        return {self.name: value}

    # Useful for enum param, having several names
    def names(self) :
        return [self.name]

    def __repr__(self):
        return self.name


class BooleanDef(ParamDef):
    """Parameter with discrete value 0 or 1"""

    def __init__(self, name, **argv):
        super(BooleanDef, self).__init__(name, ParamType.BOOL, min=0, max=1, **argv)

    def range(self, n):
        return [0, 1]

    def rand(self, alpha):
        return round(alpha)




class EnumParam(ParamDef):
    """Enum param is a facility representing a choice / switch as many boolean parameters.
    It is not itself a Sympy symbol. use #symbol("value") to access it"""

    def __init__(self, name, values: List[str], **argv):
        super(EnumParam, self).__init__(name, ParamType.ENUM, min=None, max=None, **argv)
        self.values = values

    def expandParams(self, currValue=None):
        values = self.values + [None]
        res = dict()
        for enum_val in values:
            var_name = "%s_%s" % (self.name, enum_val if enum_val is not None else "default")
            res[var_name] = 1.0 if enum_val == currValue else 0.0
        return res

    def symbol(self, enumValue):
        """Access parameter for each enum value : <paramName>_<paramValue>"""
        if enumValue is None:
            return Symbol(self.name + '_default')
        if not enumValue in self.values:
            raise Exception("enumValue should be one of %s. Was %s" % (str(self.values), enumValue))
        return Symbol(self.name + '_' + enumValue)

    def names(self) :
        return ["%s_%s" % (self.name, value) for value in (self.values + ["default"]) ]

    def rand(self, alpha):
        i = math.ceil(alpha * (len(self.values))) -1
        return self.values[int(i)]

    def range(self, n):
        return self.values


class ActivityExtended(Activity):
    """Improved API for activity : adding a few useful methods.
    Those methods are backported to #Activity in order to be directly available on all existing instances
    """

    def getExchange(self, name=None, input=None, single=True):
        """Get exchange by name or input

        Parameters
        ----------
        name : name of the exchange. Name can be suffixed with '#LOCATION' to distinguish several exchanges with same name. \
            It can also be suffised by '*' to match an exchange starting with this name. Location can be a negative match '!'
            Exampple : "Wood*#!RoW" matches any exchange with name  containing Wood, and location not "RoW"

        single :True if a single match is expected. Otherwize, a list of result is returned

        Returns
        -------
            Single exchange or list of exchanges (if _single is False or "name" contains a '*')
            raise Exception if not matching exchange found
        """

        def single_match(name, exch) :

            # Name can be "Elecricity#RER"
            if "#" in name:
                name, loc = name.split("#")
                negative = False
                if loc.startswith("!") :
                    negative = True
                    loc = loc[1:]
                act = getActByCode(*exch['input'])

                if not 'location' in act or (negative and act['location'] == loc) or (not negative and act['location'] != loc) :
                    return False

            if '*' in name :
                name = name.replace('*', '')
                return name in exch['name']
            else :
                return name == exch['name']


        def match(exch):
            if name :
                if isinstance(name, list):
                    return any(single_match(iname, exch) for iname in name)
                else:
                    return single_match(name, exch)

            if input:
                return input == exch['input']

        exchs = list(exch for exch in self.exchangesNp() if match(exch))
        if len(exchs) == 0:
            raise Exception("Found no exchange matching name : %s" % name)

        if single and len(exchs) != 1:
            raise Exception("Expected 1 exchange with name '%s' found %d" % (name, len(exchs)))
        if single:
            return exchs[0]
        else:
            return exchs

    def setOutputAmount(self, amount):
        '''Set the amount for the single output exchange (1 by default)'''
        self.addExchanges({self : amount})

    def updateExchanges(self, updates: Dict[str, any] = dict()):
        """Update existing exchanges, by name.

        Parameters
        ----------
        updates : Dict of "<exchange name>" => <new value>

            <exchange name> can be suffixed with '#LOCATION' to distinguish several exchanges with same name. \
            It can also be suffixed by '*' to match an exchange starting with this name. Location can be a negative match '!'
            Exampple : "Wood*#!RoW" matches any exchange with name  containing Wood, and location not "RoW"

            <New Value>  : either single value (float or SympPy expression) for updating only amount, or activity for updating only input,
            or dict of attributes, for updating both at once, or any other attribute.
            The amount can reference the symbol 'old_amount' that will be replaced with the current amount of the exchange.
        """

        # Update exchanges
        for name, attrs in updates.items():

            exchs = self.getExchange(name, single=not '*' in name)
            if not isinstance(exchs, list):
                exchs = [exchs]
            for exch in exchs:

                if attrs is None:
                    exch.delete()
                    exch.save()
                    continue

                # Single value ? => amount
                if not isinstance(attrs, dict):
                    if isinstance(attrs, Activity):
                        attrs = dict(input=attrs)
                    else :
                        attrs = dict(amount=attrs)

                if 'amount' in attrs:
                    attrs.update(_amountToFormula(attrs['amount'], exch['amount']))

                exch.update(attrs)
                exch.save()

                # We have a formula now ? => register it to parametrized exchange
                if 'formula' in attrs:
                    bw.parameters.add_exchanges_to_group(DEFAULT_PARAM_GROUP, self)

    def deleteExchanges(self, name, single=True):
        ''' Remove matching exchanges '''
        exchs = self.getExchange(name, single=single)
        if not isinstance(exchs, list):
            exchs = [exchs]
        if len(exchs) == 0 :
            raise Exception("No exchange found for '%s'" % name)
        for ex in exchs :
            ex.delete()
            ex.save()
        self.save()

    def substituteWithDefault(self, exchange_name: str, switch_act: Activity, paramSwitch: EnumParam, amount=None):

        """Substitutes one exchange with a switch on other activities, or fallback to the current one as default (parameter set to None)
        For this purpose, we create a new exchange referencing the activity switch, and we multiply current activity by '<param_name>_default',
        making it null as soon as one enum value is set.

        This is useful for changing electricty mix, leaving the default one if needed

        Parameters
        ----------
        act : Activity to update
        exchange_name : Name of the exchange to update
        switch_act : Activity to substitue as input
        amount : Amount of the input (uses previous amount by default)
        """

        current_exch = self.getExchange(exchange_name)

        prev_amount = amount if amount else _getAmountOrFormula(current_exch)

        self.addExchanges({switch_act: prev_amount})
        self.updateExchanges({exchange_name: paramSwitch.symbol(None) * prev_amount})

    def addExchanges(self, exchanges: Dict[Activity, Union[NumOrExpression, dict]] = dict()):
        """Add exchanges to an existing activity, with a compact syntax :

        Parameters
        ----------
        exchanges : Dict of activity => amount or activity => attributes_dict. \
            Amount being either a fixed value or Sympy expression (arithmetic expression of Sympy symbols)
        """
        parametrized = False
        for sub_act, attrs in exchanges.items():

            if isinstance(attrs, dict):
                amount = attrs.pop('amount')
            else:
                amount = attrs
                attrs = dict()

            exch = self.new_exchange(
                input=sub_act.key,
                name=sub_act['name'],
                unit=sub_act['unit'] if 'unit' in sub_act else None,
                type='production' if self == sub_act else 'biosphere' if sub_act['database'] == BIOSPHERE3_DB_NAME else 'technosphere')

            exch.update(attrs)
            exch.update(_amountToFormula(amount))
            if 'formula' in exch:
                parametrized = True

            exch.save()
        self.save()
        if parametrized:
            bw.parameters.add_exchanges_to_group(DEFAULT_PARAM_GROUP, self)

    def getAmount(self, *args, sum=False, **kargs):
        """
        Get the amount of one or several exchanges, selected by name or input. See #getExchange()
        """
        exchs = self.getExchange(*args, single=not sum, **kargs)
        if sum:
            res = 0
            if len(exchs) == 0:
                raise Exception("No exchange found")
            for exch in exchs:
                res += _getAmountOrFormula(exch)
            return res
        else:
            return _getAmountOrFormula(exchs)

    def exchangesNp(self):
        """ List of exchange, except production (output) one."""
        for exch in self.exchanges():
            if exch['input'] != exch['output']:
                yield exch


# Backport new methods to vanilla Activity class in order to benefit from it for all existing instances
for name, item in ActivityExtended.__dict__.items():
    if isinstance(item, types.FunctionType):
        setattr(Activity, name, item)


def _eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


def _isnumber(value):
    return isinstance(value, int) or isinstance(value, float)


def printAct(*activities,  **params):
    """
    Print activities and their exchanges.
    If parameter values are provided, formulas will be evaluated accordingly
    """
    tables = []
    names = []

    for act in activities:
        df = pd.DataFrame(index=['input', 'amount', 'unit', 'type'])
        data = dict()
        for (i, exc) in enumerate(act.exchanges()):
            input = bw.get_activity(exc.input.key)
            amount = _getAmountOrFormula(exc)

            # Params provided ? Evaluate formulas
            if len(params) > 0 and isinstance(amount, Basic):
                new_params = [(name, value) for name, value in _completeParamValues(params).items()]
                amount = amount.subs(new_params)

            name = exc['name']
            if 'location' in input and input['location'] != "GLO":
                name += "#%s" % input['location']
            if exc.input.key[0] not in [BIOSPHERE3_DB_NAME, ECOINVENT_DB_NAME] :
                name += " {user-db}"

            iname = name
            i=1
            while iname in data :
                iname = "%s#%d" % (name, i)
                i += 1

            data[iname] = [str(input), amount, exc.unit, exc['type']]

        for key, values in data.items() :
            df[key] = values

        tables.append(df.T)
        names.append(_actDesc(act))

    full = pd.concat(tables, axis=1, keys=names, sort=True)

    if len(activities) == 2 :
        yellow = "background-color:yellow"
        iamount1 = full.columns.get_loc((names[0], "amount"))
        iamount2 = full.columns.get_loc((names[1], "amount"))
        iact1 = full.columns.get_loc((names[0], "input"))
        iact2 = full.columns.get_loc((names[1], "input"))
        def same_amount(row) :
            res = [""] * len(row)

            if row[iamount1] != row[iamount2] :
                res[iamount1] = yellow
                res[iamount2] = yellow
            if row[iact1] != row[iact2]:
                res[iact1] = yellow
                res[iact2] = yellow
            return res

        full = full.style.apply(same_amount, axis=1)

    display(full)



def resetDb(db_name):
    """ Create or cleanup a user DB"""
    if db_name in bw.databases:
        _eprint("Db %s was here. Reseting it" % db_name)
        del bw.databases[db_name]
    db = bw.Database(db_name)
    db.write(dict())

def initDb(project_name) :
    global ECOINVENT_DB_NAME
    '''Init brightway and detect version of existing installation of ecoinvent'''
    bw.projects.set_current(project_name)
    bw.bw2setup()

    for key in bw.databases.keys() :
        if 'ecoinvent' in key :
            ECOINVENT_DB_NAME = key
            print("Found ecoinvent DB : %s" % key)
            return

    _eprint("No ecoinvent DB found")


def importDb(dbname, path):
    '''Import eco invent DB'''

    if dbname in bw.databases:
        _eprint("Database '%s' has already been imported " % dbname)
    else:
        ei34 = bw.SingleOutputEcospold2Importer(path, dbname)
        ei34.apply_strategies()
        ei34.statistics()
        ei34.write_database()


dbs = dict()

def _getDb(dbname) -> bw.Database:
    """Pool of Database instances"""
    if not dbname in dbs:
        dbs[dbname] = bw.Database(dbname)
    return dbs[dbname]


def resetParams(db_name):
    """Reset project and activity parameters"""
    _param_registry().clear()
    ProjectParameter.delete().execute()
    ActivityParameter.delete().execute()
    DatabaseParameter.delete().execute()
    Group.delete().execute()


# Index of activities per name, for fast search dict[db_name][activity_word] => list of activitites
db_index = dict()


def _split_words(name):
    clean = re.sub('[^0-9a-zA-Z]+', ' ', name)
    clean = re.sub(' +', ' ', clean)
    clean = clean.lower()

    return clean.split(' ')


def _build_index(db):
    res = defaultdict(set)
    for act in db:
        words = _split_words(act['name'])
        for word in words:
            res[word].add(act)
    return res


def _get_indexed_db(db_name):
    if not db_name in db_index:
        db_index[db_name] = _build_index(_getDb(db_name))
    return db_index[db_name]


def _find_candidates(db_name, name):

    res = []
    index = _get_indexed_db(db_name)
    words = _split_words(name)
    for word in words:
        candidates = index[word]
        if len(res) == 0 or (0 < len(candidates) < len(res)):
            res = list(candidates)
    return res


def getActByCode(db_name, code):
    """ Get activity by code """
    return _getDb(db_name).get(code)


def findActivity(name=None, loc=None, in_name=None, code=None, categories=None, category=None, db_name=None,
                 single=True, unit=None):
    """
        Find single activity by name & location
        Uses index for fast fetching
    """

    if name and '*' in name :
        in_name = name.replace("*", "")
        name = None

    def act_filter(act):
        if name and not name == act['name']:
            return False
        if in_name and not in_name in act['name']:
            return False
        if loc and not loc == act['location']:
            return False
        if unit and not unit == act['unit'] :
            return False
        if category and not category in act['categories']:
            return False
        if categories and not tuple(categories) == act['categories']:
            return False
        return True

    if code:
        acts = [getActByCode(db_name, code)]
    else:
        name_key = name if name else in_name

        # Find candidates via index
        candidates = _find_candidates(db_name, name_key)

        # Exact match
        acts = list(filter(act_filter, candidates))

    if single and len(acts) == 0:
        raise Exception("No activity found in '%s' with name '%s' and location '%s'" % (db_name, name, loc))
    if single and len(acts) > 1:
        raise Exception("Several activity found in '%s' with name '%s' and location '%s':\n%s" % (
        db_name, name, loc, str(acts)))
    if len(acts) == 1:
        return acts[0]
    else:
        return acts


def findBioAct(name=None, loc=None, **kwargs):
    """Alias for findActivity(name, ... db_name=BIOSPHERE3_DB_NAME)
    """
    return findActivity(name=name, loc=loc, db_name=BIOSPHERE3_DB_NAME, **kwargs)


def findTechAct(name=None, loc=None, **kwargs):
    """Alias for findActivity(name, ... db_name=ECOINVENT_DB_NAME)
    """
    return findActivity(name=name, loc=loc, db_name=ECOINVENT_DB_NAME, **kwargs)


def Max(a, b) :
    """Max define as algrebraic forumal with 'abs' for proper computation on vectors """
    return (a + b + abs(a - b)) / 2

def Min(a, b) :
    """Max define as algrebraic forumal with 'abs' for proper computation on vectors """
    return (a + b - abs(b - a)) / 2

def interpolate(x, x1, x2, y1, y2):
    """Build an expression for linear interpolation between two points.
    If x is not within [x1, x2] the corresponding bound Y values are returned"""
    x = Min(Max(x, x1), x2)
    return y1 + (y2 - y1) * (x - x1) / (x2 - x1)

def newInterpolatedAct(dbname: str, name: str, act1: ActivityExtended, act2: ActivityExtended, x1, x2, x, alpha1=1, alpha2=1, **kwargs):

    """Creates a new activity made of interpolation of two similar activities.
    For each exchange :
    amount = alpha1 * a1 + (x - X1) * (alpha2 * a2 - alpha1 * a1) / (x2 - x1)

    Parameters
    ----------
    name : Name of new activity
    act1 : Activity 1
    act2 : Activity 2
    x1 : X for act1
    x2 : X for act 2
    x : Should be a parameter symbol
    alpha1 : Ratio for act1 (Default value = 1)
    alpha2 : Ratio for act2 (Default value = 1)
    kwargs : Any other param will be added as attributes of new activity
    """
    res = copyActivity(dbname, act1, name, withExchanges=False, **kwargs)

    exch1_by_input = dict({exch['input']: exch for exch in act1.exchangesNp()})
    exch2_by_input = dict({exch['input']: exch for exch in act2.exchangesNp()})

    inputs = set(chain(exch1_by_input.keys(), exch2_by_input.keys()))

    for input in inputs:

        exch1 = exch1_by_input.get(input)
        exch2 = exch2_by_input.get(input)
        exch = exch1 if exch1 else exch2

        amount1 = exch1['amount'] if exch1 else 0
        amount2 = exch2['amount'] if exch2 else 0

        if exch1 and exch2 and exch1['name'] != exch2['name']:
            raise Exception("Input %s refer two different names : %s, %s" % (input, exch1['name'], exch2['name']))

        amount = interpolate(x, x1, x2, amount1 * alpha1, amount2 * alpha2)
        act = getActByCode(*input)
        res.addExchanges({act: dict(amount=amount, name=exch['name'])})
    return res




def newParamDef(name, type, **kwargs):
    """Creates a param and register it into a global registry and as a brightway parameter"""
    if type == ParamType.ENUM:
        param = EnumParam(name, **kwargs)
    elif type == ParamType.BOOL :
        param = BooleanDef(name, **kwargs)
    else:
        param = ParamDef(name, type=type, **kwargs)

    # Put it in local registry (in memory)
    if name in _param_registry():
        _eprint("Param %s was already defined : overriding" % name)
    _param_registry()[name] = param

    # Save in brightway2 project
    bwParams = [dict(name=key, amount=value) for key, value in param.expandParams().items()]
    bw.parameters.new_project_parameters(bwParams)

    return param


def newFloatParam(name, default, **kwargs):
    return newParamDef(name, ParamType.FLOAT, default=default, **kwargs)

def newBoolParam(name, default, **kwargs):
    return newParamDef(name, ParamType.BOOL, default=default, **kwargs)

def newEnumParam(name, default, **kwargs):
    return newParamDef(name, ParamType.ENUM, default=default, **kwargs)


def _amountToFormula(amount: Union[float, str, Basic], currentAmount=None):
    """Transform amount in exchange to either simple amount or formula"""
    res = dict()
    if isinstance(amount, Basic):

        if currentAmount != None:
            amount = amount.subs(old_amount, currentAmount)

        # Check the expression does not reference undefined params
        all_symbols = list([key for param in _param_registry().values() for key, val in param.expandParams().items()])
        for symbol in amount.free_symbols:
            if not str(symbol) in all_symbols:
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



def _getAmountOrFormula(ex: ExchangeDataset) -> Union[Basic, float]:
    """ Return either a fixed float value or an expression for the amount of this exchange"""
    if 'formula' in ex:
        try:
            return parse_expr(ex['formula'])
        except:
            _eprint("Error while parsing formula '%s' : backing to amount" % ex['formula'])

    return ex['amount']





def _newAct(db_name, code):

    db = _getDb(db_name)
    # Already present : delete it ?
    for act in db:
        if act['code'] == code:
            _eprint("Activity '%s' was already in '%s'. Overwriting it" % (code, db_name))
            act.delete()

    return db.new_activity(code)



def newActivity(db_name, name, unit,
        exchanges: Dict[Activity, Union[float, str]] = dict(),
        code=None,
        **argv):
    """Creates a new activity

    Parameters
    ----------
    name : Name ofthe new activity
    db_name : Destination DB : ACV DB by default
    exchanges : Dict of activity => amount. If amount is a string, is it considered as a formula with parameters
    argv : extra params passed as properties of the new activity
    """
    act = _newAct(db_name, code if code else name)
    act['name'] = name
    act['type'] = 'process'
    act['unit'] = unit
    act.update(argv)

    # Add exchanges
    act.addExchanges(exchanges)

    return act


def copyActivity(db_name, activity: ActivityExtended, code=None , withExchanges=True, **kwargs) -> ActivityExtended:
    """Copy activity into a new DB"""

    res = _newAct(db_name, code)

    for key, value in activity.items():
        if key not in ['database', 'code']:
            res[key] = value
    for k, v in kwargs.items():
        res._data[k] = v
    res._data[u'code'] = code
    res['name'] = code
    res.save()

    if withExchanges:
        for exc in activity.exchanges():
            data = deepcopy(exc._data)
            data['output'] = res.key
            # Change `input` for production exchanges
            if exc['input'] == exc['output']:
                data['input'] = res.key
            ExchangeDataset.create(**dict_as_exchangedataset(data))

    return res



def newSwitchAct(dbname, name, paramDef: ParamDef, acts_dict: Dict[str, Activity]):
    """Create a new parametrized, virtual activity, made of a map of other activities, controlled by an enum parameter.
    This enables to implement a "Switch" with brightway parameters
    Internally, this will create a linear sum of other activities controlled by <param_name>_<enum_value> : 0 or 1

    Parameters
    ----------
    paramDef : parameter definition of type enum
    acts_dict : dict of <enumValue> => activity
    """

    # Transform map of enum values to correspoding formulas <param_name>_<enum_value>
    exch = {act: paramDef.symbol(key) for key, act in acts_dict.items()}
    res = newActivity(
        dbname,
        name,
        unit=list(acts_dict.values())[0]['unit'],
        exchanges=exch)


    # Unit of switch activity is the one of the children
    for key, act in acts_dict.items():
        if 'unit' in act:
            res['unit'] = act['unit']
            res.save()
    return res


def _actName(act: Activity):
    """Generate pretty name for activity, appending location if not 'GLO' """
    res = act['name']
    if act['location'] != 'GLO':
        res += "[%s]" % act["location"]
    return res

def _actDesc(act: Activity):
    """Generate pretty name for activity + basic information """
    name = _actName(act)
    amount = 1
    for ex in act.exchanges() :
        if ex['type'] == 'production' :
            amount = ex['amount']

    return "%s (%f %s)" % (name, amount, act['unit'])

def _multiLCA(activities, methods):
    """Simple wrapper around brightway API"""
    bw.calculation_setups['process'] = {'inv': activities, 'ia': methods}
    lca = bw.MultiLCA('process')
    cols = [_actName(act) for act_amount in activities for act, amount in act_amount.items()]
    return pd.DataFrame(lca.results.T, index=[method_name(method) for method in methods], columns=cols)


def _listOfDictToDictOflist(LD):
    return {k: [dic[k] for dic in LD] for k in LD[0]}


def _completeParamValues(params):
    """Check parameters and expand enum params.

    Returns
    -------
        Dict of param_name => float value
    """

    # undef_params = param_registry.keys() - params.keys()
    # if undef_params :
    #    raise Exception("Some model parameters are not set : %s" % undef_params)

    res = dict()
    for key, val in params.items():
        if key in _param_registry():
            param = _param_registry()[key]
        else:
            raise Exception("Parameter not found : %s. Valid parameters : %s" % (key, list(_param_registry().keys())))

        if isinstance(val, list):
            newvals = [param.expandParams(val) for val in val]
            res.update(_listOfDictToDictOflist(newvals))
        else:
            res.update(param.expandParams(val))
    return res


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




def preMultiLCAAlgebric(model:ActivityExtended, methods, param_names=None, amount=1) :
    '''
        This method transforms an activity into a set of functions ready to compute LCA very fast on a set on methods.
        You may use is and pass the result to postMultiLCAAlgebric for fast computation on a model that does not change.

        This method is used by multiLCAAlgebric
    '''

    # print("computing model to expression for %s" % model)
    expr, actBySymbolName = actToExpression(model)

    dbname = model.key[0]

    param_names = _expand_param_names(param_names)

    # Check missing params
    free_names = set([str(symb) for symb in expr.free_symbols])
    act_names = set([str(symb) for symb in actBySymbolName.keys()])
    expected_names  = free_names - act_names
    missing_names = expected_names - set(param_names)

    if len(missing_names) > 0 :
        raise Exception('Missing parameter of the model : %s' % str(missing_names))

    debug(param_names)

    for name in param_names :
        registry_names = _expand_param_names(_param_registry())
        if name not in registry_names :
            raise Exception('Model refers to unknown param "%s"' % name)

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
    lambdas = []
    for imethod, method in enumerate(methods):
        # print("Generating lamba function for %s / %s" % (model, method))

        # Replace activities by their value in expression for this method
        sub = dict({symbol: lca.iloc[imethod, iact] for iact, symbol in enumerate(pureTechActBySymbol.keys())})
        method_expr = expr.xreplace(sub)

        # Tranform Sympy expression to lambda function, based on numpy to fast vectorial evaluation
        lambd = lambdify(param_names, method_expr, 'numpy')
        lambdas.append(lambd)

    return lambdas

def method_name(method) :
    return method[1] + " - " + method[2]

def postMultiLCAAlgebric(methods, lambdas, alpha=1, **params):
    '''
        Compute LCA for a given set of parameters and pre-compiled lambda functions.
        This function is used by **multiLCAAlgebric**

        Parameters
        ----------
        methodAndLambdas : Output of preMultiLCAAlgebric
        **params : Parameters of the model
    '''

    # Check and expand params
    params = _completeParamValues(params)

    # Expand parameters as list of parameters
    param_length = 1

    for key, val in params.items():
        if isinstance(val, list):
            if param_length == 1:
                param_length = len(val)
            elif param_length != len(val):
                raise Exception("Parameters should be a single value or a list of same number of values")

    # Expand params and transform lists to np.array for vector computation
    for key in params.keys():
        val = params[key]
        if not isinstance(val, list):
            val = list([val] * param_length)
        params[key] = np.array(val)


    res = np.zeros((len(methods), param_length))

    # Compute result on whole vectors of parameter samples at a time : lambdas use numpy for vector computation
    for imethod, lambd in enumerate(lambdas):
        res[imethod, :] = alpha * lambd(**params)

    return pd.DataFrame(res, index=[method_name(method) for method in methods]).transpose()

def _expand_param_names(param_names) :
    '''Expand parameters names (with enum params) '''
    return [name for key in param_names for name in _param_registry()[key].names()]

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
        if type(model) is tuple :
            model, alpha=model

        lambdas = preMultiLCAAlgebric(model, methods, params.keys())

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


def actToExpression(act: Activity):
    """Computes a symbolic expression of the model, referencing background activities and model parameters as symbols

    Returns
    -------
        (sympy_expr, dict of symbol => activity)
    """

    act_symbols = dict()  # Dict of  act = > symbol

    def act_to_symbol(db_name, code):

        act = _getDb(db_name).get(code)
        name = act['name']
        base_slug = slugify(name, separator='_')

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
            if input_db in [BIOSPHERE3_DB_NAME, ECOINVENT_DB_NAME]:
                if not (input_db, input_code) in act_symbols:
                    act_symbols[(input_db, input_code)] = act_to_symbol(input_db, input_code)
                act_expr = act_symbols[(input_db, input_code)]

            # Our model : recursively transform it to a symbolic expression
            else:

                if input_db == act['database'] and input_code == act['code']:
                    raise Exception("Recursive exchange : %s" % (act.__dict__))

                sub_act = _getDb(input_db).get(input_code)
                act_expr = rec_func(sub_act)

            res += formula * act_expr

        return  res / outputAmount

    expr = rec_func(act)

    return (expr, _reverse_dict(act_symbols))


def _reverse_dict(dic):
    return {v: k for k, v in dic.items()}


def _heatmap(df, title, vmax, ints=False):
    ''' Produce heatmap of a dataframe'''
    fig, ax = plt.subplots(figsize=(17, 17))
    sns.heatmap(df.transpose(), cmap="gist_heat_r", vmax=vmax, annot=True, fmt='.0f' if ints else 'f', square=True)
    plt.title(title, fontsize=20)
    ax.tick_params(axis="x", labelsize=18)
    ax.tick_params(axis="y", labelsize=18)


def oat_matrix(model, impacts, n=10) :
    '''Generates a heatmap of the incertitude of the model, varying input parameters one a a time.'''

    # Compile model into lambda functions for fast LCA
    lambdas = preMultiLCAAlgebric(model, impacts, _param_registry().keys())

    change = np.zeros((len(_variable_params()), len(impacts)))

    for iparam, param in enumerate(_variable_params().values()) :
        params = {param.name: param.default for param in _param_registry().values()}

        # Compute range of values for given param
        params[param.name] = param.range(n)

        # Compute LCA
        df = postMultiLCAAlgebric(impacts, lambdas, **params)

        # Compute change
        change[iparam] =  (df.max() - df.min()) / df.median() * 100

    # Build final heatmap
    change = pd.DataFrame(change, index=_variable_params().keys(), columns=[imp[2] for imp in impacts])
    _heatmap(change.transpose(), 'Change of impacts per variability of the input parameters (%)', 100, ints=True)


def _method_unit(method) :
    return bw.Method(method).metadata['unit']


def _display_tabs(titlesAndContentF) :
    '''Generate tabs'''
    tabs = []
    titles= []
    for title, content_f in titlesAndContentF :
        titles.append(title)

        tab = widgets.Output()
        with tab :
            content_f()
        tabs.append(tab)

    res = widgets.Tab(children=tabs)
    for i, title in enumerate(titles) :
        res.set_title(i, title)
    display(res)


def oat_dasboard(modelOrLambdas, impacts, param: ParamDef, n=10) :
    '''
    Analyse the evolution of impacts for a single parameter. The other parameters are set to their default values.

    Parameters
    ----------
    model : activity, or lambdas as precomputed by preMultiLCAAlgebric, for faster computation
    impacts : set of methods
    param: parameter to analyse
    n: number of samples of the parameter
    '''

    params = {param.name : param.default for param in _param_registry().values()}

    # Compute range of values for given param
    params[param.name] = param.range(n)


    #print("Params: ", params)

    if isinstance(modelOrLambdas, Activity) :
        df = multiLCAAlgebric(modelOrLambdas, impacts, **params)
    else :
        df = postMultiLCAAlgebric(impacts, modelOrLambdas, **params)

    # Add X values in the table
    pname = param.name
    if param.unit :
        pname = '%s [%s]' % (pname, param.unit)
    df.insert(0, pname, param.range(n))
    df = df.set_index(pname)


    def table() :
        display(df)

    def graph() :

        with warnings.catch_warnings():

            warnings.simplefilter("ignore")

            nb_rows = len(impacts) // 3 + 1

            fig, axes = plt.subplots(figsize=(15, 15))

            axes = df.plot(
                ax=axes, sharex=True, subplots=True,
                layout=(nb_rows, 3),
                #legend=None,
                kind = 'line' if param.type == ParamType.FLOAT else 'bar')

            axes = axes.flatten()

            for ax, impact in zip(axes, impacts) :
                ax.set_ylim(ymin=0)
                ax.set_ylabel(_method_unit(impact))

            plt.show(fig)

    def change() :

        ch = (df.max() - df.min()) / df.median() * 100
        fig, ax = plt.subplots(figsize=(9, 6))
        plt.title('Relative change for %s' % df.index.name)
        ch.plot(kind='barh', rot=30)
        ax.set_xlabel('Relative change of the median value (%)')
        plt.tight_layout()
        plt.show(fig)

    _display_tabs([
        ("Graphs", graph),
        ("Data", table),
        ("Variation", change)
    ])


def oat_dashboard_interact(model, methods):
    '''Interative dashboard, with a dropdown for selecting parameter'''

    lambdas = preMultiLCAAlgebric(model, methods, _param_registry().keys())

    def process_func(param) :
        oat_dasboard(lambdas, methods, _param_registry()[param])

    paramlist = list(_variable_params().keys())
    interact(process_func, param=paramlist)


def _stochastics(modelOrLambdas, methods, n=1000) :

    ''' Compute stochastic impacts for later analysis of incertitude '''

    # Extract variable names
    param_names = list(_variable_params().keys())
    problem = {
        'num_vars': len(param_names),
        'names': param_names,
        'bounds': [[0, 1]] * len(param_names)
    }

    print("Generating samples ...")
    X = saltelli.sample(problem, n, calc_second_order=True)

    # Map normalized 0-1 random values into real values
    print("Transforming samples ...")
    params = dict()
    for i, param_name in enumerate(param_names) :
        param = _param_registry()[param_name]
        vals = list(map(lambda v : param.rand(v), X[:, i]))
        params[param_name] = vals

    # Add static parameters
    for param in _fixed_params().values() :
        params[param.name] = param.default

    print("Processing LCA ...")
    if isinstance(modelOrLambdas, Activity):
        Y = multiLCAAlgebric(modelOrLambdas, methods, **params)
    else:
        Y = postMultiLCAAlgebric(methods, modelOrLambdas, **params)

    return problem, X, Y


def _variable_params():
    return {key : param for key, param in _param_registry().items() if param.distrib != DistributionType.FIXED}

def _fixed_params():
    return {key : param for key, param in _param_registry().items() if param.distrib == DistributionType.FIXED}

def _sobols(methods, problem, Y) :
    ''' Computes sobols indices'''
    s1 = np.zeros((len(problem['names']), len(methods)))
    st = np.zeros((len(problem['names']), len(methods)))

    for i, method in enumerate(methods) :

        try:
            y = Y[Y.columns[i]]
            res = sobol.analyze(problem, y.to_numpy(), calc_second_order=True)
            st[:, i] = res["ST"]
            s1[:, i] = res["S1"]

        except Exception as e:
            error("Sobol failed on %s" % method[2], e)
    return (s1, st)

def _incer_stochastic_matrix(methods, param_names, Y, st):

    ''' Internal method computing matrix of parameter importance '''
    def draw(mode) :

        if mode == 'sobol' :
            data = st
        else :
            # If percent, express result as percentage of standard deviation / mean
            data = np.zeros((len(param_names), len(methods)))
            for i, method in enumerate(methods):
                # Total variance
                var = np.var(Y[Y.columns[i]])
                mean = np.mean(Y[Y.columns[i]])
                if mean != 0 :
                    data[:, i] = np.sqrt((st[:, i] * var)) / mean * 100


        df = pd.DataFrame(data, index=param_names, columns=[method_name(method) for method in methods])
        _heatmap(
            df.transpose(),
            title="Relative deviation of impacts (%)" if mode == 'percent' else "Sobol indices (part of variability)",
            vmax=100 if mode == 'percent' else 1,
            ints= mode == 'percent')

    interact(draw, mode=[('Raw sobol indices (ST)', 'sobol'), ('Deviation (ST) / mean', 'percent')])


def incer_stochastic_matrix(modelOrLambdas, methods, n=1000):
    ''' Method computing matrix of parameter importance '''

    problem, X, Y = _stochastics(modelOrLambdas, methods, n)

    print("Processing Sobol indices ...")
    s1, st = _sobols(methods, problem, Y)

    _incer_stochastic_matrix(methods, problem['names'], Y, st)


def _incer_stochastic_violin(methods, Y) :
    ''' Internal method for computing violin graph of impacts '''

    nb_rows = math.ceil(len(methods) / 3)
    fig, axes = plt.subplots(nb_rows, 3, figsize=(15, 15), sharex=True)

    for imethod, method, ax in zip(range(len(methods)), methods, axes.flatten()) :
        ax.violinplot(Y[Y.columns[imethod]], showmedians=True)
        ax.title.set_text(method_name(method))
        ax.set_ylim(ymin=0)
        ax.set_ylabel(_method_unit(method))

    plt.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
    plt.show(fig)

def incer_stochastic_violin(modelOrLambdas, methods, n=1000):

    ''' Method for computing violin graph of impacts '''

    problem, X, Y = _stochastics(modelOrLambdas, methods, n)

    _incer_stochastic_violin(methods, Y)

def _incer_stochastic_variations(methods, Y, param_names, sobols1):

    ''' Method for computing violin graph of impacts '''
    method_names=[method_name(method) for method in methods]

    std = np.std(Y)
    mean = np.mean(Y)

    fig = plt.figure(num=None, figsize=(12, 6), dpi=80, facecolor='w', edgecolor='k')
    ax = plt.gca()
    tab20b = plt.get_cmap('tab20b')
    tab20c = plt.get_cmap('tab20c')
    ax.set_prop_cycle('color', [tab20b(k) if k < 1 else tab20c(k-1) for k in np.linspace(0, 2, 40)])

    relative_variance_pct =  std*std / (mean*mean) * 100
    totplt = plt.bar(np.arange(len(method_names)), relative_variance_pct, 0.8)

    sum = np.zeros(len(methods))

    plots = [totplt[0]]

    data = np.zeros((len(param_names) + 2, len(methods)))
    data[0, :] = mean
    data[1, :] = std


    for i_param, param_name in enumerate(param_names) :
        s1 = sobols1[i_param, :]
        data[i_param+2, :] = s1

        curr_bar = s1 * relative_variance_pct
        curr_plt = plt.bar(np.arange(len(method_names)), curr_bar, 0.8, bottom=sum)
        sum += curr_bar
        plots.append(curr_plt[0])


    plt.legend(plots, ['Higher order'] + param_names)
    plt.xticks(np.arange(len(method_names)), method_names, rotation=90)
    plt.title("variance / mean² (%)")
    plt.show(fig)

    # Show raw data
    rows = ["mean", "std"] + ["s1(%s)" % param for param in param_names]
    df = pd.DataFrame(data, index=rows, columns=[method_name(method) for method in methods])
    display(df)




def incer_stochastic_dasboard(model, methods, n=1000) :
    ''' Generates a dashboard with several statistics : matrix of parameter incertitude, violin diagrams, ...'''

    problem, X, Y = _stochastics(model, methods, n)
    param_names = problem['names']

    print("Processing Sobol indices ...")
    s1, st = _sobols(methods, problem, Y)


    def violin() :
        _incer_stochastic_violin(methods, Y)

    def variation():
        _incer_stochastic_variations(methods, Y, param_names, s1)

    def matrix() :
        _incer_stochastic_matrix(methods, problem['names'], Y, st)

    _display_tabs([
        ("Violin graphs", violin),
        ("Impact variations", variation),
        ("Sobol matrix", matrix)
    ])
