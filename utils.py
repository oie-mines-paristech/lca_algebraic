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
ECOINVENT_DB_NAME = 'ecoinvent 3.4 cut off'
BIOSPHERE3_DB_NAME = 'biosphere3'
ACV_DB_NAME = 'incer-acv'

DEFAULT_PARAM_GROUP = "acv"


# Global
def param_registry():
    """ """
    if not 'param_registry' in builtins.__dict__:
        builtins.param_registry = dict()
    return builtins.param_registry


# Sympy symbols
old_amount = symbols(
    "old_amount")  # Can be used in epxression of amount for updateExchanges, in order to reference the previous value


NumOrExpression = Union[float, Basic]


# Type of parameters
class ParamType:
    ENUM = "enum"
    BOOL = "bool"
    FLOAT = "float"


class ParamDef(Symbol):
    """Generic definition of a parameter, with name, bound, type, distribution
    This definition will serve both to generate brightway2 parameters and to evaluate.

    This class inherits sympy Symbol, making it possible to use in standard arithmetic python
    while keeping it as a symbolic expression (delayed evaluation).
    """

    def __new__(cls, name, *karg, **kargs):
        return Symbol.__new__(cls, name)

    def __init__(self, name, type: str, default, description=""):
        self.name = name
        self.type = type
        self.default = default
        self.description = description

    # Expand parameter (usefull for enum param)
    def expandParams(self, value=None) -> Dict[str, float]:
        if value == None:
            value = self.default
        return {self.name: value}

    def __repr__(self):
        return self.name


class EnumParam(ParamDef):
    """Enum param is a facility representing a choice / switch as many 0/1 parameters.
    It is not itself a Sympy symbol. use #symbol("value") to access it"""

    def __init__(self, name, values: List[str], **argv):
        super(EnumParam, self).__init__(name, ParamType.ENUM, **argv)
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


class ActivityExtended(Activity):
    """Improved API for activity : adding a few useful methods.
    Those methods are backported to #Activity in order to be directly available on all existing instances
    """

    def getExchange(self, name=None, input=None, single=True):
        """Get exchange by name or input

        Parameters
        ----------
        name : name of the exchange. Name can be suffixed with '#LOCATION' to distinguish several exchanges with same name. \
            It can also be suffised by '*' to match on exchange starting with this name (Default value = None)
        single :True if a single match is expected. Otherwize, a list of result is returned

        Returns
        -------
            Single exchange or list of exchanges (if _single is False or "name" contains a '*')
            raise Exception if not matching exchange found
        """

        _single=single

        def single_match(name, exch) :
            nonlocal _single
            # Name can be "Elecricity#RER"
            if "#" in name:
                name, loc = name.split("#")
                act = getActByCode(*exch['input'])
                if not 'location' in act or not act['location'] == loc :
                    return False
            if '*' in name :
                _single=False
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

        if _single and len(exchs) != 1:
            raise Exception("Expected 1 exchange with name '%s' found %d" % (name, len(exchs)))
        if _single:
            return exchs[0]
        else:
            return exchs

    def setOutputAmount(self, amount):
        self.addExchanges({self : amount})

    def updateExchanges(self, updates: Dict[str, any] = dict()):
        """Update existing exchanges, by name.

        Parameters
        ----------
        updates : Dict of exchange name => either single value (float or SympPy expression) for updating only amount, \
            or dict of attributes, for updating several at a time. The sympy expression can reference the symbol 'old_amount' \
            that will be replaced with the current value.
        """

        # Update exchanges
        for name, attrs in updates.items():

            exchs = self.getExchange(name)
            if not isinstance(exchs, list):
                exchs = [exchs]
            for exch in exchs:

                if attrs is None:
                    exch.delete()
                    exch.save()
                    continue

                # Single value ? => amount
                if not isinstance(attrs, dict):
                    attrs = dict(amount=attrs)

                if 'amount' in attrs:
                    attrs.update(amountToFormula(attrs['amount'], exch['amount']))

                exch.update(attrs)
                exch.save()

                # We have a formula now ? => register it to parametrized exchange
                if 'formula' in attrs:
                    bw.parameters.add_exchanges_to_group(DEFAULT_PARAM_GROUP, self)

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

        prev_amount = amount if amount else getAmountOrFormula(current_exch)

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
            exch.update(amountToFormula(amount))
            if 'formula' in exch:
                parametrized = True

            exch.save()
        self.save()
        if parametrized:
            bw.parameters.add_exchanges_to_group(DEFAULT_PARAM_GROUP, self)

    def getAmount(self, *args, sum=False, **kargs):
        """
        Get the amount of one or several exchanges, selected by name or input. See #get_dict_as_exchange()
        """
        exchs = self.getExchange(*args, single=not sum, **kargs)
        if sum:
            res = 0
            if len(exchs) == 0:
                raise Exception("No exchange found")
            for exch in exchs:
                res += getAmountOrFormula(exch)
            return res
        else:
            return getAmountOrFormula(exchs)

    def exchangesNp(self):
        """ """
        for exch in self.exchanges():
            if exch['input'] != exch['output']:
                yield exch


# Backport new methods to vanilla Activity class in order to benefit from it for all existing instances
for name, item in ActivityExtended.__dict__.items():
    if isinstance(item, types.FunctionType):
        setattr(Activity, name, item)


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


def isnumber(value):
    return isinstance(value, int) or isinstance(value, float)


def printAct(*activities,  **params):
    """
    Print activities and their exchanges.
    If parameter values are provided, formulas will be evaluated accordingly
    """
    tables = []
    names = []

    if 'hide_same' in params :
        hide_same = params.pop('hide_same')
    else :
        hide_same = False

    for act in activities:
        df = pd.DataFrame(index=['input', 'amount', 'unit', 'type'])
        data = dict()
        for (i, exc) in enumerate(act.exchanges()):
            input = bw.get_activity(exc.input.key)
            amount = getAmountOrFormula(exc)

            # Params provided ? Evaluate formulas
            if len(params) > 0 and isinstance(amount, Basic):
                new_params = [(name, value) for name, value in completeParamValues(params).items()]
                amount = amount.subs(new_params)

            name = exc['name']
            if 'location' in input and input['location'] != "GLO":
                name += "[%s]" % input['location']
            if exc.input.key[0] == ACV_DB_NAME :
                name += " {incer}"

            iname = name
            i=1
            while iname in data :
                iname = "%s#%d" % (name, i)
                i += 1

            data[iname] = [str(input), amount, exc.unit, exc['type']]

        for key, values in data.items() :
            df[key] = values

        tables.append(df.T)
        names.append(actDesc(act))

    full = pd.concat(tables, axis=1, keys=names, sort=True)
    if hide_same :
        full = full.loc[full[names[0]]['amount'] != full[names[1]]['amount']]

    display(full)



def resetDb(db_name=ACV_DB_NAME):
    """ Create or cleanup DB """
    if db_name in bw.databases:
        eprint("Db %s was here. Reseting it" % db_name)
        del bw.databases[db_name]
    db = bw.Database(ACV_DB_NAME)
    db.write(dict())


def importDb(dbname, path):

    if dbname in bw.databases:
        eprint("Database '%s' has already been imported " % dbname)
    else:
        ei34 = bw.SingleOutputEcospold2Importer(path, dbname)
        ei34.apply_strategies()
        ei34.statistics()
        ei34.write_database()


dbs = dict()


def getDb(dbname) -> bw.Database:
    """Pool of Database instances"""
    if not dbname in dbs:
        dbs[dbname] = bw.Database(dbname)
    return dbs[dbname]


def resetParams(db_name=ACV_DB_NAME):
    """Reset project and activity parameters"""
    param_registry().clear()
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
        db_index[db_name] = _build_index(getDb(db_name))
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
    """Get activity by code"""
    return getDb(db_name).get(code)


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


def interpolate(x, x1, x2, y1, y2):
    """Build an expression for linear interpolation between two points"""
    return y1 + (y2 - y1) * (x - x1) / (x2 - x1)


def newInterpolatedAct(name: str, act1: ActivityExtended, act2: ActivityExtended, x1, x2, x, alpha1=1, alpha2=1, **kwargs):

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
    res = copyActivity(act1, name, withExchanges=False, **kwargs)

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
    else:
        param = ParamDef(name, type=type, **kwargs)

    # Put it in local registry (in memory)
    if name in param_registry():
        eprint("Param %s was already defined : overriding" % name)
    param_registry()[name] = param

    # Save in brightway2 project
    bwParams = [dict(name=key, amount=value) for key, value in param.expandParams().items()]
    bw.parameters.new_project_parameters(bwParams)

    return param


def newFloatParam(name, default, **kwargs):
    return newParamDef(name, ParamType.FLOAT, default=default, **kwargs)


def newEnumParam(name, default, **kwargs):
    return newParamDef(name, ParamType.ENUM, default=default, **kwargs)


def amountToFormula(amount: Union[float, str, Basic], currentAmount=None):
    """Transform amount in exchange to either simple amount or formula"""
    res = dict()
    if isinstance(amount, Basic):

        if currentAmount != None:
            amount = amount.subs(old_amount, currentAmount)

        # Check the expression does not reference undefined params
        all_symbols = list([key for param in param_registry().values() for key, val in param.expandParams().items()])
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










def getAmountOrFormula(ex: ExchangeDataset) -> Union[Basic, float]:
    """ Return either a fixed float value or an expression for the amount of this exchange"""
    if 'formula' in ex:
        try:
            return parse_expr(ex['formula'])
        except:
            eprint("Error while parsing formula '%s' : backing to amount" % ex['formula'])

    return ex['amount']





def _newAct(code, db_name=ACV_DB_NAME):

    db = getDb(db_name)
    # Already present : delete it ?
    for act in db:
        if act['code'] == code:
            eprint("Activity '%s' was already in '%s'. Overwriting it" % (code, db_name))
            act.delete()

    return db.new_activity(code)



def newActivity(
        name, unit,
        exchanges: Dict[Activity, Union[float, str]] = dict(),
        db_name=ACV_DB_NAME,
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
    act = _newAct(code if code else name , db_name)
    act['name'] = name
    act['type'] = 'process'
    act['unit'] = unit
    act.update(argv)

    # Add exchanges
    act.addExchanges(exchanges)

    return act


def copyActivity(activity: ActivityExtended, code=None, db_name=ACV_DB_NAME, withExchanges=True, **kwargs) -> ActivityExtended:
    """Copy activity into a new DB"""

    res = _newAct(code, db_name)

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



def newSwitchAct(name, paramDef: ParamDef, acts_dict: Dict[str, Activity]):
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
        name,
        unit=list(acts_dict.values())[0]['unit'],
        exchanges=exch)


    # Unit of switch activity is the one of the children
    for key, act in acts_dict.items():
        if 'unit' in act:
            res['unit'] = act['unit']
            res.save()
    return res


def actName(act: Activity):
    """Generate pretty name for activity, appending location if not 'GLO' """
    res = act['name']
    if act['location'] != 'GLO':
        res += "[%s]" % act["location"]
    return res

def actDesc(act: Activity):
    """Generate pretty name for activity + basic information """
    name = actName(act)
    amount = 1
    for ex in act.exchanges() :
        if ex['type'] == 'production' :
            amount = ex['amount']

    return "%s (%f %s)" % (name, amount, act['unit'])

def _multiLCA(activities, methods):
    """Simple wrapper around brightway API"""
    bw.calculation_setups['process'] = {'inv': activities, 'ia': methods}
    lca = bw.MultiLCA('process')
    cols = [actName(act) for act_amount in activities for act, amount in act_amount.items()]
    return pd.DataFrame(lca.results.T, index=[method[2] for method in methods], columns=cols)


def listOfDictToDictOflist(LD):
    return {k: [dic[k] for dic in LD] for k in LD[0]}


def completeParamValues(params):
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
        if key in param_registry():
            param = param_registry()[key]
        else:
            raise Exception("Parameter not found : %s. Valid parameters : %s" % (key, list(param_registry().keys())))

        if isinstance(val, list):
            newvals = [param.expandParams(val) for val in val]
            res.update(listOfDictToDictOflist(newvals))
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
    params = completeParamValues(params)

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


def multiLCAAlgebric(models, methods, **params):
    """Compute LCA by expressing the foreground model as symbolic expression of background activities and parameters.
    Then, compute 'static' inventory of the referenced background activities.
    This enables a very fast recomputation of LCA with different parameters, useful for stochastic evaluation of parametrized model

    Parameters
    ----------
    models : Single model or list of models : if list of models, you cannot use param lists
    methods : List of methods / impacts to consider
    params : You should provide named values of all the parameters declared in the model. \
             Values can be single value or list of samples, all of the same size
    """
    dfs = dict()

    if not isinstance(models, list):
        models = [models]

    # Check and expand params
    params = completeParamValues(params)

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

    for model in models:

        #print("computing model to expression for %s" % model)
        expr, actBySymbolName = actToExpression(model)

        # Check missing parameters
        param_names = set([str(symb) for symb in expr.free_symbols])
        act_names = set([str(symb) for symb in actBySymbolName.keys()])
        remaining = param_names - act_names - params.keys()
        if remaining:
            raise Exception("Some model parameters are not set : %s" % remaining)

        # Create dummy reference to biosphere
        # We cannot run LCA to biosphere activities
        # We create a technosphere activity mapping exactly to 1 biosphere item
        pureTechActBySymbol = OrderedDict()
        for name, act in actBySymbolName.items():
            if act[0] == BIOSPHERE3_DB_NAME:
                act = _getOrCreateDummyBiosphereActCopy(act[1])
            else:
                act = getActByCode(*act)
            pureTechActBySymbol[name] = act

        # List of activities, ordered
        acts = pureTechActBySymbol.values()

        # Transform to [{act1:1], {act2:1}, etc] for MultiLCA
        actsWithAmount = [{act: 1} for act in acts]

        #print("Computing full LCA of background activities for '%s'" % model)

        # Compute LCA for all background activities and methods
        lca = _multiLCA(actsWithAmount, methods)

        # For each method, compute an algebric expression with activities replaced by their values
        lambdas = []
        for imethod, method in enumerate(methods):

            #print("Generating lamba function for %s / %s" % (model, method))

            # Replace activities by their value in expression for this method
            sub = [(symbol, lca.iloc[imethod, iact]) for iact, symbol in enumerate(pureTechActBySymbol.keys())]
            method_expr = expr.subs(sub)

            # Tranform Sympy expression to lambda function, based on numpy to fast vectorial evaluation
            lambd = lambdify(params.keys(), method_expr, 'numpy')
            lambdas.append(lambd)

        res = np.zeros((len(methods), param_length))

        # Compute result on whole vectors of parameter samples at a time : lambdas use numpy for vector computation
        for imethod, lambd in enumerate(lambdas):
            #print("Computing results for %s / %s" % (model, methods[imethod]))
            res[imethod, :] = lambd(**params)

        model_name = actName(model)
        df = pd.DataFrame(res, index=[method[2] for method in methods]).transpose()

        # Single params ? => give the single row the name of the model activity
        if df.shape[0] == 1:
            df = df.rename(index={0: model_name})

        dfs[model_name] = df

    if len(dfs) == 1:
        df = list(dfs.values())[0]
        display(df)
    else:
        # Concat several dataframes for several models
        display(pd.concat(list(dfs.values())))


def _getOrCreateDummyBiosphereActCopy(code):
    """
        We cannot reference directly biosphere in the model, since LCA can only be applied to products
        We create a dummy activity in our DB, with same code, and single exchange of '1'
    """

    code_to_find = code + "#asTech"
    try:
        return getDb(ACV_DB_NAME).get(code_to_find)
    except:
        bioAct = getDb(BIOSPHERE3_DB_NAME).get(code)
        name = bioAct['name'] + ' # asTech'
        res = newActivity(name, bioAct['unit'], {bioAct: 1}, code=code_to_find)
        return res


def actToExpression(act: Activity):
    """Computes a symbolic expression of the model, referencing background activities and model parameters as symbols

    Returns
    -------
        (sympy_expr, dict of symbol => activity)
    """

    act_symbols = dict()  # Dict of  act = > symbol

    def act_to_symbol(db_name, code):

        act = getDb(db_name).get(code)
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

            formula = getAmountOrFormula(exch)

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
            if input_db != ACV_DB_NAME:
                if not (input_db, input_code) in act_symbols:
                    act_symbols[(input_db, input_code)] = act_to_symbol(input_db, input_code)
                act_expr = act_symbols[(input_db, input_code)]

            # Our model : recursively transform it to a symbolic expression
            else:

                if input_db == act['database'] and input_code == act['code']:
                    raise Exception("Recursive exchange : %s" % (act.__dict__))

                sub_act = getDb(input_db).get(input_code)
                act_expr = rec_func(sub_act)

            res += formula * act_expr

        return  res / outputAmount

    expr = rec_func(act)

    return (expr, reverse_dict(act_symbols))


def reverse_dict(dic):
    return {v: k for k, v in dic.items()}



# ---   Graphviz tools -------

from graphviz import Digraph

def graph(act:ActivityExtended) :
    g = Digraph(node_attr=dict(
        shape ='plaintext',
        rankdir='LR'))

    aName = actName(act)

    lines = ['<tr><td colspan="2" > %s </td></tr>' % aName]

    for exch in act.exchanges() :
        amount = getAmountOrFormula(exch)
        name = exch["name"]
        lines.append('<tr><td>%s</td><td label="foo" port="%s">%s</td></tr>' % (name, name, str(amount)))

    g.node(aName, "<<table>%s</table>>" % "\n".join(lines))

    display(g)