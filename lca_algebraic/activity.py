import re
from collections import defaultdict
from copy import deepcopy
from types import FunctionType
from typing import Dict, Tuple, Union

import brightway2 as bw
import pandas as pd
from bw2data.backends.peewee import Activity, ExchangeDataset
from bw2data.backends.peewee.utils import dict_as_exchangedataset
from pint import DimensionalityError, Quantity
from sympy import Basic, simplify, symbols

from lca_algebraic.base_utils import _actName, _getDb, _isOutputExch
from lca_algebraic.database import (
    _find_biosphere_db,
    _isForeground,
    _listTechBackgroundDbs,
    with_db_context,
)
from lca_algebraic.params import (
    DbContext,
    ParamDef,
    _complete_and_expand_params,
    _getAmountOrFormula,
    _param_registry,
)

from .base_utils import ValueOrExpression, getActByCode
from .log import logger, warn
from .settings import Settings
from .units import is_dimensionless, is_equivalent, parse_db_unit
from .units import unit_registry as u

# Can be used in expression of amount for updateExchanges, in order to reference the previous value
old_amount = symbols("old_amount")
old_amount_with_unit = u.Quantity(old_amount, u.old_unit)


def _exch_name(exch):
    return exch["name"] if "name" in exch else str(exch.input)


class ActivityExtended(Activity):
    """Improved API for activity : adding a few useful methods.
    Those methods are backported to #Activity in order to be directly available on all existing instances
    """

    @with_db_context
    def listExchanges(self):
        """Iterates on all exchanges (except "production") and return a list of (exch-name, target-act, amount)"""
        res = []
        for exc in self.exchanges():
            # Don't show production
            if _isOutputExch(exc):
                continue

            input = bw.get_activity(exc.input.key)
            amount = _getAmountOrFormula(exc)
            res.append((exc["name"], input, amount))
        return res

    @with_db_context
    def getExchange(self, name: str = None, input: Activity = None, single=True):
        """Get exchange by name or input

        Parameters
        ----------
        name :
            name of the exchange. Name can be suffixed with '#LOCATION' to distinguish several exchanges with same name. \
            It can also be suffised by '*' to match an exchange starting with this name. Location can be a negative match '!'
            Example  "Wood*#!RoW" matches any exchange with name  containing Wood, and location not "RoW"

        input :
            input activity

        single :
            True if a single match is expected. Otherwize, a list of result is returned

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
                act = getActByCode(*exch["input"])

                if "location" not in act or (negative and act["location"] == loc) or (not negative and act["location"] != loc):
                    return False

            if "*" in name:
                name = name.replace("*", "")
                return name in _exch_name(exch)
            else:
                return name == _exch_name(exch)

        def match(exch):
            if name:
                if isinstance(name, list):
                    return any(single_match(iname, exch) for iname in name)
                else:
                    return single_match(name, exch)

            if input:
                return input == exch["input"]

        exchs = list(exch for exch in self.non_production_exchanges() if match(exch))
        if len(exchs) == 0:
            raise Exception("Found no exchange matching name : %s" % name)

        if single and len(exchs) != 1:
            raise Exception("Expected 1 exchange with name '%s' found %d" % (name, len(exchs)))
        if single:
            return exchs[0]
        else:
            return exchs

    def setOutputAmount(self, amount):
        """Set the amount for the single output exchange (1 by default)"""

        output_exchange = self.getOutputExchange()
        output_exchange["amount"] = amount
        output_exchange.save()
        output_exchange.save()

    @with_db_context
    def updateExchanges(self, updates: Dict[str, any] = dict()):
        """Update existing exchanges, by name.

        Parameters
        ----------
        updates : Dict of "<exchange name>" => <new value>

            <exchange name> can be suffixed with '#LOCATION' to distinguish several exchanges with same name. \
            It can also be suffixed by '*' to match an exchange starting with this name. Location can be a negative match '!'
            Example : "Wood*#!RoW" matches any exchange with name  containing Wood, and location not "RoW"

            <New Value>  : either single value (float or SympPy expression) for updating only amount, \
                or activity for updating only input,
            or dict of attributes, for updating both at once, or any other attribute.
            The amount can reference the symbol 'old_amount' that will be replaced with the current amount of the exchange.
        """

        # Update exchanges
        for ex_name, updates in updates.items():
            # Build input & amount
            if updates is not None and not isinstance(updates, dict):
                if isinstance(updates, Activity):
                    updates = dict(input=updates)
                else:
                    updates = dict(amount=updates)

            # Find echzanges matching name
            matching_exchanges = self.getExchange(ex_name, single="*" not in ex_name)
            if not isinstance(matching_exchanges, list):
                matching_exchanges = [matching_exchanges]

            # Loop on matching echanges to update
            for exch in matching_exchanges:
                # Delete exchange
                if updates is None:
                    exch.delete()
                    exch.save()
                    continue
                else:
                    self._update_exchange(exch, updates)

    def deleteExchanges(self, name, single=True):
        """Remove matching exchanges

        Parameters
        ----------
        name:
            Name of the exchange to delete. Can contain wildcards. See #getExchange for more details.

        single:
            If true (default) expect to only delete a single exchange
        """
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
    def addExchanges(self, exchanges: Dict[Activity, Union[ValueOrExpression, dict]] = dict()):
        """Add exchanges to an existing activity, with a compact syntax :

        Parameters
        ----------
        exchanges :
            Dict of activity => amount or activity => attributes_dict. \
            Amount being either a fixed value or Sympy expression (arithmetic expression of Sympy symbols)
        """

        with DbContext(self.key[0]):
            for sub_act, updates in exchanges.items():
                if not isinstance(updates, dict):
                    updates = dict(amount=updates)

                exch = self.new_exchange(
                    input=sub_act.key,
                    name=sub_act["name"],
                    unit=sub_act["unit"] if "unit" in sub_act else None,
                    type="technosphere" if sub_act.get("type") == "process" else "biosphere",
                )

                self._update_exchange(exch, updates)

            self.save()

    def _transform_unit(self, amount: ValueOrExpression, exchange_unit: str):
        if not Settings.units_enabled:
            return amount

        # Parse unit
        if exchange_unit is None:
            logger.warn("Missing unit for target activity; Assuming dimensionless")
            exchange_unit = u.dimensionless
        else:
            exchange_unit = parse_db_unit(exchange_unit)

        if not isinstance(amount, Quantity):
            if not is_dimensionless(exchange_unit) and not self.isSwitch() and amount != 0:
                raise Exception(
                    f"Unit '{exchange_unit}' expected, and dimensionless amount provided in a non-switch Activity: {amount}"
                )
            else:
                return amount

        # Help the compiler
        amount: Quantity

        act_unit = parse_db_unit(self["unit"])

        if self.isSwitch() and (exchange_unit != act_unit):
            raise Exception(f"Units should be the same in a switch activity {exchange_unit} != {act_unit}")

        # Using 'old_amount' ? => replace with requested unit
        if amount.units == u.old_unit:
            amount = u.Quantity(amount.magnitude, exchange_unit)

        # We try to transform either to the exchange unit, or ex_unit/act_unit
        for target_unit in [exchange_unit, exchange_unit / act_unit]:
            # We'll catch error at the end if none of the target units worked
            if not is_equivalent(target_unit, amount.units):
                continue

            # Try to convert
            new_amount = amount.to(target_unit).magnitude

            # Auto scale disabld ?
            if not _equals(amount.magnitude, new_amount) and not u.auto_scale:
                raise Exception(f"auto_scale is disabled. '{amount}' should be explicity transformed to {target_unit}")

            return new_amount

        # At this point, no convertion worked :
        raise DimensionalityError(
            amount.units,
            exchange_unit,
            f"Unit of amount '{amount}' is not compatible with physical unit of exchange '{exchange_unit}' or the "
            f"unit of exchange divided by the unit of the activity : '{exchange_unit / act_unit}'",
        )

    def _amount_to_formula(self, amount: ValueOrExpression, exchange: ExchangeDataset):
        res = dict()
        if isinstance(amount, Basic):
            current_amount = exchange.get("amount", None)
            if current_amount is not None:
                amount = amount.subs(old_amount, current_amount)

            # Check the expression does not reference undefined params
            all_symbols = list([key for param in _param_registry().values() for key, val in param.expandParams().items()])
            for symbol in amount.free_symbols:
                if not str(symbol) in all_symbols:
                    raise Exception("Symbol '%s' not found in params : %s" % (symbol, all_symbols))

            res["formula"] = str(amount)
            res["amount"] = 0
        elif isinstance(amount, float) or isinstance(amount, int):
            res["amount"] = amount
        else:
            raise Exception(
                "Amount should be either a constant number or a Sympy expression (expression of ParamDef). Was : %s"
                % type(amount)
            )
        return res

    def _update_exchange(self, exchange: ExchangeDataset, updates):
        """Update a single exchange. Take care of setting amount / formula accordingly"""
        amount = updates.pop("amount") if "amount" in updates else None

        if amount is not None:
            # Update units
            amount = self._transform_unit(amount, exchange["unit"])

            # Extract formula if two separate field "amount" and "formula"
            # Update the list of updates
            updates.update(self._amount_to_formula(amount, exchange))

        exchange.update(updates)
        exchange.save()

    def isSwitch(self):
        return self.get("switch", False)

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

    def getOutputExchange(self):
        for exch in self.exchanges():
            if (exch["input"] == exch["output"]) and (exch["type"] == "production"):
                return exch

    def getOutputAmount(self):
        """Return the amount of the production : 1 if none is found"""
        output_exchange = self.getOutputExchange()
        return 1.0 if output_exchange is None else output_exchange["amount"]

    def non_production_exchanges(self):
        """List of exchange, except production (output) one."""
        for exch in self.exchanges():
            if exch["input"] != exch["output"]:
                yield exch

    def updateMeta(self, **kwargs):
        """Update any property. Useful to update axes"""
        for key, val in kwargs.items():
            self._data[key] = val
        self.save()


def findActivity(
    name=None,
    loc=None,
    code=None,
    categories=None,
    category=None,
    db_name=None,
    single=True,
    case_sensitive=False,
    unit=None,
    limit=1500,
) -> ActivityExtended:
    """
        Find activity by name & location
        Uses index for fast fetching

    :param name: Name of the activity. Can contain '*' for searching partial chain
    :param loc: optional location
    :param code: Unique code. If provided alone, returns the activity for this code
    :param categories: Optional : exact list of catagories
    :param category: Optional : single category that should be part of the list of categories of the selected activities
    :param db_name: Name of the database
    :param single: If False, returns a list of matching activities. If True (default) fails if more than one activity fits.
    :param case_sensitive: If True (default) ignore the case
    :param unit: If provided, only match activities with provided unit
    :return: Either a single activity (if single is True) or a list of activities, possibly empty.
    """

    in_name = None

    if name and "*" in name:
        in_name = name.replace("*", "")
        name = None

    if not case_sensitive:
        if name:
            name = name.lower()
        if in_name:
            in_name = in_name.lower()

    def act_filter(act):
        act_name = act["name"]
        if not case_sensitive:
            act_name = act_name.lower()

        if name and not name == act_name:
            return False
        if in_name and in_name not in act_name:
            return False
        if loc and not loc == act["location"]:
            return False
        if unit and not unit == act["unit"]:
            return False
        if category and category not in act["categories"]:
            return False
        if categories and not tuple(categories) == tuple(act["categories"]):
            return False
        return True

    if code:
        acts = [getActByCode(db_name, code)]
    else:
        search = name if name is not None else in_name

        search = search.lower()
        search = search.replace(",", " ")

        # Find candidates via index
        # candidates = _find_candidates(db_name, name_key)
        candidates = _getDb(db_name).search(search, limit=limit)

        if len(candidates) == 0:
            # Try again removing strange caracters
            search = re.sub(r"\w*[^a-zA-Z ]+\w*", " ", search)
            candidates = _getDb(db_name).search(search, limit=limit)

        # Exact match
        acts = list(filter(act_filter, candidates))

    if single and len(acts) == 0:
        any_name = name if name else in_name
        raise Exception("No activity found in '%s' with name '%s' and location '%s'" % (db_name, any_name, loc))
    if single and len(acts) > 1:
        raise Exception(
            "Several activity found in '%s' with name '%s' and location '%s':\n%s"
            % (db_name, name, loc, "\n".join(str(act) for act in acts))
        )
    if len(acts) == 1:
        return acts[0]
    else:
        return acts


def findBioAct(name=None, loc=None, **kwargs):
    """Alias for findActivity(name, ... db_name=BIOSPHERE3_DB_NAME). See doc for #findActivity"""
    return findActivity(name=name, loc=loc, db_name=_find_biosphere_db(), **kwargs)


def findTechAct(name=None, loc=None, single=True, **kwargs):
    """
    Search activities in technosphere. This function tries to guess which database is your background database.
    If you have more than one background technosphere, you should use findActivity() and specify the **db_name** directly.
    See also doc for #findActivity"""
    dbs = _listTechBackgroundDbs()
    if len(dbs) > 1:
        raise Exception(
            "There is more than one technosphere background DB (%s) please use findActivity(..., db_name=YOUR_DB)" % str(dbs)
        )

    return findActivity(name=name, loc=loc, db_name=dbs[0], single=single, **kwargs)


def _equals(val1: ValueOrExpression, val2: ValueOrExpression):
    """Compare float of Sympy values"""
    if val1 == val2:
        return True
    if isinstance(val1, Basic) != isinstance(val2, Basic):
        return False
    return simplify(val1 / val2) == 1.0


def _newAct(db_name, code):
    if not _isForeground(db_name):
        warn(
            "WARNING: You are creating activity in background DB. You should only do it in your foreground / user DB : ",
            db_name,
        )

    db = _getDb(db_name)
    # Already present : delete it ?
    for act in db:
        if act["code"] == code:
            warn("Activity '%s' was already in '%s'. Overwriting it" % (code, db_name))
            act.delete()

    return db.new_activity(code)


def newActivity(
    db_name,
    name,
    unit,
    exchanges: Dict[Activity, Union[float, str]] = dict(),
    amount=1,
    code=None,
    type="process",
    switchActivity=False,
    **argv,
) -> ActivityExtended:
    """Creates a new activity

    Parameters
    ----------
    name :
        Name of the new activity
    db_name :
        Destination DB : ACV DB by default
    unit:
        Unit of the process

    code:
        Unique code in the Db. Optional. If not provided, the name is used

    exchanges :
        Dict of activity => amount. See the doc for @addExchanges()

    argv :
        Any extra params passed as properties of the new activity

    switch:
        Activities marked as *switch* are expected to be linear combination of activities of same unit.
        This option changes how physical units are checked.

    amount:
        Production amount. 1 by default
    """

    code = code if code else name

    act = _newAct(db_name, code)
    act["name"] = name
    act["type"] = type
    act["unit"] = unit
    if switchActivity:
        act["switch"] = True

    act.update(argv)

    # Add single production exchange
    if type == "process":
        ex = act.new_exchange(
            input=act.key,
            name=act["name"],
            unit=act["unit"],
            type="production",
            amount=amount,
        )
        ex.save()

        act["reference product"] = act["name"]
        act.save()

    # Add exchanges
    act.addExchanges(exchanges)

    return act


def copyActivity(db_name, activity: ActivityExtended, code=None, withExchanges=True, **kwargs) -> ActivityExtended:
    """Copy an activity and its exchanges into another database. You usually want to copy activities from your background to
    your foreground DB to update them, keeping your background DB clean.

    Parameters
    ----------
    db_name:
        Name of the target database
    activity:
        Source activity
    code:
        Code of the target activity. Also used as its name


    Returns
    -------
        The new activity. note that is is flagged with the custom property **inherited_from**, providing the full key of the
        initial activity.
    """

    res = _newAct(db_name, code)

    # Same code if not provided
    if code is None:
        code = activity.key[1]

    for key, value in activity.items():
        if key not in ["database", "code"]:
            res[key] = value
    for k, v in kwargs.items():
        res._data[k] = v
    res._data["code"] = code
    res["name"] = code
    res["type"] = "process"
    res["inherited_from"] = activity.key
    res.save()

    if withExchanges:
        for exc in activity.exchanges():
            data = deepcopy(exc._data)
            data["output"] = res.key
            # Change `input` for production exchanges
            if exc["input"] == exc["output"]:
                data["input"] = res.key
            ExchangeDataset.create(**dict_as_exchangedataset(data))

    return res


ActivityOrActivityAmount = Union[Activity, Tuple[Activity, float]]


def newSwitchAct(dbname, name, paramDef: ParamDef, acts_dict: Dict[str, ActivityOrActivityAmount]):
    """Creates a new parametrized, virtual activity, made of a map of other activities, controlled by an enum parameter.
    This enables to implement a "Switch" with brightway parameters
    Internally, this will create a linear sum of other activities controlled by <param_name>_<enum_value> : 0 or 1

    By default, all activities have associated amount of 1.
    You can provide other amounts by providing a tuple of (activity, amount).

    Parameters
    ----------
    dbname:
        name of the target DB
    name:
        Name of the new activity
    paramDef :
        parameter definition of type enum
    acts_dict :
        dict of "enumValue" => activity or "enumValue" => (activity, amount)

    Examples
    --------

    >>> newSwitchAct(MYDB, "switchAct", switchParam, {
    >>>    "val1" : act1 # Amount is 1
    >>>    "val2" : (act2, 0.4) # Different amount
    >>>    "val3" : (act3, b + 6) # Amount with formula
    >>> }
    """

    # Transform map of enum values to corresponding formulas <param_name>_<enum_value>
    exch = defaultdict(lambda: 0)

    # Forward last unit as unit of the switch
    unit = None
    for key, act in acts_dict.items():
        amount = 1
        if isinstance(act, (list, tuple)):
            act, amount = act
        exch[act] += amount * paramDef.symbol(key)
        unit = act["unit"]

    res = newActivity(dbname, name, unit=unit, exchanges=exch, switch=True)

    return res


def _actDesc(act: ActivityExtended):
    """Generate pretty name for activity + basic information"""
    name = _actName(act)
    amount = act.getOutputAmount()

    return "%s (%f %s)" % (name, amount, act["unit"])


def printAct(*activities, **params):
    """
    Print activities and their exchanges.
    If parameter values are provided, formulas will be evaluated accordingly.

    Parameters
    ----------
    activities:
        One or two activities. If two activities are provided, differences are highlighted.

    params:
        If provided, the formulas are evaluated accordingly and the result amount is shown instead of formula

    Returns
    -------
        A Dataframe is returned, containing all information, exchange by exchange
    """
    tables = []
    names = []

    for act in activities:
        with DbContext(act.key[0]):
            inputs_by_ex_name = dict()
            df = pd.DataFrame(index=["input", "amount", "unit"])
            data = dict()
            for i, exc in enumerate(act.exchanges()):
                # Don't show production
                if _isOutputExch(exc):
                    continue

                input = bw.get_activity(exc.input.key)
                amount = _getAmountOrFormula(exc)

                # Params provided ? Evaluate formulas
                if len(params) > 0 and isinstance(amount, Basic):
                    new_params = [(name, value) for name, value in _complete_and_expand_params(params).items()]
                    amount = amount.subs(new_params)

                ex_name = _exch_name(exc)
                # if 'location' in input and input['location'] != "GLO":
                #    name += "#%s" % input['location']
                # if exc.input.key[0] not in [BIOSPHERE3_DB_NAME, ECOINVENT_DB_NAME()]:
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
            names.append(_actDesc(act))

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

            if row.iloc[iamount1] != row.iloc[iamount2]:
                res[iamount1] = yellow
                res[iamount2] = yellow
            if row.iloc[iact1] != row.iloc[iact2]:
                res[iact1] = yellow
                res[iact2] = yellow
            return res

        full = full.style.apply(same_amount, axis=1)

    return full


# Backport new methods to vanilla Activity class in order to benefit from it for all existing instances
for name, item in ActivityExtended.__dict__.items():
    if isinstance(item, FunctionType):
        setattr(Activity, name, item)
