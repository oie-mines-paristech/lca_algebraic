import functools
import inspect
import re
import types

from copy import deepcopy
from itertools import chain

import pandas as pd
from bw2data.backends.peewee.utils import dict_as_exchangedataset
from bw2data.meta import databases as dbmeta
from sympy import symbols

from .base_utils import *
from .base_utils import _getDb, _actDesc, _getAmountOrFormula, _actName, _isOutputExch
from .params import *
from .params import _param_registry, _completeParamValues
from typing import Tuple, Dict
import inspect


BIOSPHERE3_DB_NAME="biosphere3"

_metaCache = defaultdict(lambda : {})

def _setMeta(dbname, key, value) :
    """Set meta param on DB"""
    _metaCache[dbname][key] = value

    data = dbmeta[dbname]
    data[key] = value
    dbmeta[dbname] = data
    dbmeta.flush()


def _getMeta(db_name, key) :

    if key in _metaCache[db_name] :
        return _metaCache[db_name][key]

    val = dbmeta[db_name].get(key)
    _metaCache[db_name][key] = val
    return val

FOREGROUND_KEY = "fg"

def _isForeground(db_name) :
    """ Check is db is marked as foreground DB : which means activities may be parametrized / should be developped. """
    return _getMeta(db_name, FOREGROUND_KEY)

def setForeground(db_name) :
    """ Set a db as being a foreground database, meaning it is parametrized and lca_algebraic should develop its activities"""
    return _setMeta(db_name, FOREGROUND_KEY, True)

def setBackground(db_name) :
    """ Set a db as being a foreground database, meaning it should be considred as static"""
    return _setMeta(db_name, FOREGROUND_KEY, False)

def SET_USER_DB(db_name) :
    """Deprecated, use #setForeground() / setBackground() instead"""
    error("Deprecated, use #setForeground() / setBackground() instead")
    setForeground(db_name)

def _listTechBackgroundDbs() :
    """List all background databases technosphere (non biosphere) batabases"""
    return list(name for name in bw.databases if not _isForeground(name) and not name == BIOSPHERE3_DB_NAME)

old_amount = symbols("old_amount")  # Can be used in expression of amount for updateExchanges, in order to reference the previous value
NumOrExpression = Union[float, Basic]


def list_databases() :
    """List of databases and their status"""
    data = list(
        dict(
            name=name,
            backend=_getMeta(name, "backend"),
            type="foreground" if _isForeground(name) else "background") for name in bw.databases)
    res = pd.DataFrame(data)
    return res.set_index("name")


def with_db_context(func=None, arg="self"):
    """ Internal decorator wrapping function into DbContext, using its first parameters (either Activity, Db or Db name)"""

    if func is None:
        return functools.partial(with_db_context, arg=arg)

    param_specs = inspect.signature(func).parameters

    if not arg in param_specs :
        raise Exception("No param %s in signature of %s" % (arg, func))

    @functools.wraps(func)
    def wrapper(*args, **kwargs):

        # Transform all parameters (positionnal and named) to named ones
        all_param = {
            k: args[n] if n < len(args) else v.default
            for n, (k, v) in enumerate(param_specs.items()) if k != 'kwargs'
        }
        all_param.update(kwargs)

        val = all_param[arg]
        if hasattr(val, "key") :
            # value is an activity
            dbname = val.key[0]
        elif isinstance(val, str) :
            # Value is directly a  db_name
            dbname = val
        else:
            raise Exception("Param %s is neither an Activity or a db_name : %s" % (arg, val))

        with DbContext(dbname) :
            return func(*args, **kwargs)

    return wrapper


class ActivityExtended(Activity):
    """Improved API for activity : adding a few useful methods.
    Those methods are backported to #Activity in order to be directly available on all existing instances
    """

    @with_db_context
    def listExchanges(self) :
        """ Iterates on all exchanges (except "production") and return a list of (exch-name, target-act, amount) """
        res = []
        for exc in self.exchanges():

            # Don't show production
            if _isOutputExch(exc) :
                continue

            input = bw.get_activity(exc.input.key)
            amount = _getAmountOrFormula(exc)
            res.append((exc["name"], input, amount))
        return res

    @with_db_context
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

        def single_match(name, exch):

            # Name can be "Elecricity#RER"
            if "#" in name:
                name, loc = name.split("#")
                negative = False
                if loc.startswith("!"):
                    negative = True
                    loc = loc[1:]
                act = getActByCode(*exch['input'])

                if not 'location' in act or (negative and act['location'] == loc) or (
                        not negative and act['location'] != loc):
                    return False

            if '*' in name:
                name = name.replace('*', '')
                if "name" in exch:
                    return name in exch['name']
                else:
                    return name == exch.input
           
        def match(exch):
            if name:
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
        self.addExchanges({self: amount})

    @with_db_context
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

        parametrized = False

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
                    else:
                        attrs = dict(amount=attrs)

                if 'amount' in attrs:
                    attrs.update(_amountToFormula(attrs['amount'], exch['amount']))

                exch.update(attrs)
                exch.save()

                # We have a formula now ? => register it to parametrized exchange
                if 'formula' in attrs:
                    parametrized = True


    def deleteExchanges(self, name, single=True):
        ''' Remove matching exchanges '''
        exchs = self.getExchange(name, single=single)
        if not isinstance(exchs, list):
            exchs = [exchs]
        if len(exchs) == 0:
            raise Exception("No exchange found for '%s'" % name)
        for ex in exchs:
            ex.delete()
            ex.save()
        self.save()

    @with_db_context
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

    @with_db_context
    def addExchanges(self, exchanges: Dict[Activity, Union[NumOrExpression, dict]] = dict()):
        """Add exchanges to an existing activity, with a compact syntax :

        Parameters
        ----------
        exchanges : Dict of activity => amount or activity => attributes_dict. \
            Amount being either a fixed value or Sympy expression (arithmetic expression of Sympy symbols)
        """

        with DbContext(self.key[0]) :

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
                    type='production' if self == sub_act else 'technosphere' if sub_act.get('type') == 'process' else 'biosphere')

                exch.update(attrs)
                exch.update(_amountToFormula(amount))
                if 'formula' in exch:
                    parametrized = True

                exch.save()
            self.save()


    @with_db_context
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

    def getOutputAmount(self):

        """ Return the amount of the production : 1 if none is found """
        res = 1
        for exch in self.exchanges() :
            if exch['input'] == exch['output']:
                # Not 1 ?
                if exch['amount'] != 1:
                    res = exch['amount']
                break
        return res


    def exchangesNp(self):
        """ List of exchange, except production (output) one."""
        for exch in self.exchanges():
            if exch['input'] != exch['output']:
                yield exch


# Backport new methods to vanilla Activity class in order to benefit from it for all existing instances
for name, item in ActivityExtended.__dict__.items():
    if isinstance(item, types.FunctionType):
        setattr(Activity, name, item)


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


# Index of activities per name, for fast search dict[db_name][activity_word] => list of activitites
db_index = dict()


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
                 single=True, unit=None) -> ActivityExtended :
    """
        Find single activity by name & location
        Uses index for fast fetching
    """

    if name and '*' in name:
        in_name = name.replace("*", "")
        name = None

    def act_filter(act):
        if name and not name == act['name']:
            return False
        if in_name and not in_name in act['name']:
            return False
        if loc and not loc == act['location']:
            return False
        if unit and not unit == act['unit']:
            return False
        if category and not category in act['categories']:
            return False
        if categories and not tuple(categories) == act['categories']:
            return False
        return True

    if code:
        acts = [getActByCode(db_name, code)]
    else:
        search = name if name is not None else in_name

        search = search.lower()
        search = search.replace(',', ' ')

        # Find candidates via index
        # candidates = _find_candidates(db_name, name_key)
        candidates = _getDb(db_name).search(search, limit=200)

        if len(candidates) == 0 :
            # Try again removing strange caracters
            search = re.sub(r'\w*[^a-zA-Z ]+\w*', ' ', search)
            candidates = _getDb(db_name).search(search, limit=200)

        # print(search, candidates)

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
    """Alias for findActivity(name, ... db_name=BIOSPHERE3_DB_NAME). See doc for #findActivity """
    return findActivity(name=name, loc=loc, db_name=BIOSPHERE3_DB_NAME, **kwargs)


def findTechAct(name=None, loc=None, single=True, **kwargs):
    """
        Search activities in technosphere. This function try to guess which database is your background database.
        See also doc for #findActivity """
    dbs = _listTechBackgroundDbs()
    if len(dbs) > 1 :
        raise Exception("There is more than one technosphere background DB (%s) please use findActivity(..., db_name=YOUR_DB)" % str(dbs))

    return findActivity(name=name, loc=loc, db_name=dbs[0], single=single, **kwargs)



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


def _newAct(db_name, code):

    if not _isForeground(db_name) :
        error("WARNING: You are creating activity in background DB. You should only do it in your foreground / user DB : ", db_name)

    db = _getDb(db_name)
    # Already present : delete it ?
    for act in db:
        if act['code'] == code:
            error("Activity '%s' was already in '%s'. Overwriting it" % (code, db_name))
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


def copyActivity(db_name, activity: ActivityExtended, code=None, withExchanges=True, **kwargs) -> ActivityExtended:
    """Copy activity into a new DB"""

    res = _newAct(db_name, code)

    for key, value in activity.items():
        if key not in ['database', 'code']:
            res[key] = value
    for k, v in kwargs.items():
        res._data[k] = v
    res._data[u'code'] = code
    res['name'] = code
    res['type'] = 'process'
    res['inherited_from'] = activity.key
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

ActivityOrActivityAmount = Union[Activity, Tuple[Activity, float]]
def newSwitchAct(dbname, name, paramDef: ParamDef, acts_dict: Dict[str, ActivityOrActivityAmount]):
    """Create a new parametrized, virtual activity, made of a map of other activities, controlled by an enum parameter.
    This enables to implement a "Switch" with brightway parameters
    Internally, this will create a linear sum of other activities controlled by <param_name>_<enum_value> : 0 or 1

    By default, all activities have associated amount of 1.
    You can provide other amounts by providing a tuple of (activity, amount).

    Parameters
    ----------
    dbname: name of the target DB
    name: Name of the new activity
    paramDef : parameter definition of type enum
    acts_dict : dict of "enumValue" => activity or "enumValue" => (activity, amount)

    Examples
    --------

    >>> newSwitchAct(MYDB, "switchAct", switchParam, {
    >>>    "val1" : act1 # Amount is 1
    >>>    "val2" : (act2, 0.4) # Different amount
    >>>    "val3" : (act3, b + 6) # Amount with formula
    >>> }
    """

    # Transform map of enum values to corresponding formulas <param_name>_<enum_value>
    exch = defaultdict(lambda : 0)

    # Forward last unit as unit of the switch
    unit= None
    for key, act in acts_dict.items() :
        amount = 1
        if type(act) == list or type(act) == tuple :
            act, amount = act
        exch[act] += amount * paramDef.symbol(key)
        unit = act["unit"]

    res = newActivity(
        dbname,
        name,
        unit=unit,
        exchanges=exch)

    return res


def printAct(*args, impact=None, **params):
    """
    Print activities and their exchanges.
    If parameter values are provided, formulas will be evaluated accordingly.
    If impact is provided it will be computed.
    """
    tables = []
    names = []

    activities = args


    for iact, act in enumerate(activities):
        with DbContext(act.key[0]) :
            inputs_by_ex_name = dict()
            df = pd.DataFrame(index=['input', 'amount', 'unit'])
            data = dict()
            for (i, exc) in enumerate(act.exchanges()):

                # Don't show production
                if _isOutputExch(exc) :
                    continue

                input = bw.get_activity(exc.input.key)
                amount = _getAmountOrFormula(exc)

                # Params provided ? Evaluate formulas
                if len(params) > 0 and isinstance(amount, Basic):
                    new_params = [(name, value) for name, value in _completeParamValues(params).items()]
                    amount = amount.subs(new_params)
                if "name" in exc :
                    ex_name = exc['name']
                else: 
                    ex_name=str(exc.input)
                    
                #if 'location' in input and input['location'] != "GLO":
                #    name += "#%s" % input['location']
                #if exc.input.key[0] not in [BIOSPHERE3_DB_NAME, ECOINVENT_DB_NAME()]:
                #    name += " {user-db}"

                # Unique name : some exchanges may havve same names
                _name = ex_name
                i = 1
                while ex_name in data:
                    ex_name = "%s#%d" % (_name, i)
                    i += 1

                inputs_by_ex_name[ex_name] = input

                input_name = _actName(input)
                if _isForeground(input.key[0]):
                    input_name += "{FG}"

                data[ex_name] = [input_name, amount, exc.unit]

            # Provide impact calculation if impact provided


            for key, values in data.items():
                df[key] = values

            tables.append(df.T)

            # ensure the name is unique
            name = _actDesc(act) + ("#%d" % iact if iact > 0 else "")
            names.append(name)

    full = pd.concat(tables, axis=1, keys=names, sort=True)

    # Highlight differences in case two activites are provided
    if len(activities) == 2:
        yellow = "background-color:yellow"
        iamount1 = full.columns.get_loc((names[0], "amount"))
        iamount2 = full.columns.get_loc((names[1], "amount"))
        iact1 = full.columns.get_loc((names[0], "input"))
        iact2 = full.columns.get_loc((names[1], "input"))

        def same_amount(row):
            res = [""] * len(row)

            if row[iamount1] != row[iamount2]:
                res[iamount1] = yellow
                res[iamount2] = yellow
            if row[iact1] != row[iact2]:
                res[iact1] = yellow
                res[iact2] = yellow
            return res

        full = full.style.apply(same_amount, axis=1)

    display(full)


def newInterpolatedAct(dbname: str, name: str, act1: ActivityExtended, act2: ActivityExtended, x1, x2, x, alpha1=1,
                       alpha2=1, **kwargs):
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



def findMethods(search=None, mainCat=None) :

    """
    Find impact method. Search in all methods against a list of match strings.
    Each parameter can be either an exact match match, or case insenstive search, if suffixed by '*'

    Parameters
    ----------
    search : String to search
    mainCat : if specified, limits the research for method[0] == mainCat.
    """
    res = []
    search = search.lower()
    for method in bw.methods:
        text = str(method).lower()
        match  = search in text
        if mainCat :
            match = match and (mainCat == method[0])
        if match :
            res.append(method)
    return res
