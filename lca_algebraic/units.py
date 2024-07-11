import operator
from contextlib import contextmanager

import brightway2 as bw
from pint import DimensionalityError, OffsetUnitCalculusError, Unit, UnitRegistry
from pint.compat import _to_magnitude, zero_or_nan
from pint.facets.plain import PlainQuantity

from lca_algebraic.base_utils import getActByCode
from lca_algebraic.settings import Settings

NEW_UNITS = {"person", "old_unit"}  # Used with 'old_amount'

CORE_ALIASES = {
    "ratio": "count",
    "fraction": "count",
    "unit": "count",
    "ton": "metric_ton",  # Override silly US ton to metric ton (f*ck imperial unit system ...)
}

ALIASES = {"square_meter": "m²"}


# The global Unit registry
unit_registry = UnitRegistry()


def check_unit_consistency(db_name: str):
    """
    Check units of exchanges VS units of target activities in a single database.
    This check is done statically. The purpose is to run this on a background, non parametric, database."""
    db = bw.Database(db_name)

    errors = list()
    for act in db:
        for ex in act.exchanges():
            sub_act = getActByCode(*ex["input"])

            if sub_act["unit"] != ex["unit"]:
                errors.append(f"Unit of exxchange {ex} ({ex['unit']}) does not match unit of {sub_act} ({sub_act['unit']})")

    if errors:
        raise Exception("Error : unit inconsistancy\n" + "\n".join(errors))

    return True


def define_separate_unit(unit_name):
    """Define completely new units / dimension that cannot be added to anything else : like 'kwp' 'panel'"""
    unit_registry.define(f"{unit_name} = [{unit_name}]")


def define_alias_unit(unit_name, expression):
    """Define alias / shortcut for exiting units
    Examples:
         >>> define_alias_unit("square_meter", u.meter * u.meter)
         >>> define_alias_unit("tkm", u.ton * u.kilometer)
    """
    unit_registry.define(f"{unit_name} = {expression}")


def parse_db_unit(unit_str):
    """Convert unit found in Ecoinvent datanbase into Pint unit"""
    unit_str = unit_str.replace("-", " ").replace("standard", "")

    return unit_registry.parse_units(unit_str)


def is_dimensionless(unit):
    try:
        unit_registry.convert(1, unit, unit_registry.dimensionless)
        return True
    except DimensionalityError:
        return False


def is_equivalent(first: Unit, second: Unit):
    try:
        unit_registry.convert(1, first, second)
        return True
    except DimensionalityError:
        return False


# Add some additional units

for unit in NEW_UNITS:
    define_separate_unit(unit)

for key, val in CORE_ALIASES.items():
    unit_registry.define(f"@alias {val} = {key}")

for key, val in ALIASES.items():
    define_alias_unit(key, val)

# Hack of Pint unit to fail when autoscale is false
unit_registry.auto_scale = False


def _add_sub_modified(self: PlainQuantity, other, op):
    """Perform addition or subtraction operation and return the result.

    Parameters
    ----------
    other : pint.PlainQuantity or any type accepted by :func:`_to_magnitude`
        object to be added to / subtracted from self
    op : function
        operator function (e.g. operator.add, operator.isub)
    """

    def _safe_other_magnitude(units):
        if not self._REGISTRY.auto_scale:
            raise (Exception(f"Auto scale disabled : explicit convertio of '{other}' to {units} required"))
        return other.to(units).magnitude

    if not self._check(other):
        # other not from same Registry or not a PlainQuantity
        if zero_or_nan(other, True):
            # If the other value is 0 or NaN (but not a PlainQuantity)
            # do the operation without checking units.
            # We do the calculation instead of just returning the same
            # value to enforce any shape checking and type casting due to
            # the operation.
            units = self._units
            magnitude = op(
                self._magnitude,
                _to_magnitude(other, self.force_ndarray, self.force_ndarray_like),
            )
        elif self.dimensionless:
            units = self.UnitsContainer()

            other_magnitude = self.to(units)._magnitude

            if other_magnitude != self.magnitude and not self._REGISTRY.auto_scale:
                raise (Exception(f"Auto scale disabled : explicit convertion of '{self}' to {units} required"))

            magnitude = op(
                self.to(units)._magnitude,
                _to_magnitude(other, self.force_ndarray, self.force_ndarray_like),
            )
        else:
            raise DimensionalityError(self._units, "dimensionless")
        return self.__class__(magnitude, units)

    if not self.dimensionality == other.dimensionality:
        raise DimensionalityError(self._units, other._units, self.dimensionality, other.dimensionality)

    # Next we define some variables to make if-clauses more readable.
    self_non_mul_units = self._get_non_multiplicative_units()
    is_self_multiplicative = len(self_non_mul_units) == 0
    if len(self_non_mul_units) == 1:
        self_non_mul_unit = self_non_mul_units[0]
    other_non_mul_units = other._get_non_multiplicative_units()
    is_other_multiplicative = len(other_non_mul_units) == 0
    if len(other_non_mul_units) == 1:
        other_non_mul_unit = other_non_mul_units[0]

    # Presence of non-multiplicative units gives rise to several cases.
    if is_self_multiplicative and is_other_multiplicative:
        if self._units == other._units:
            magnitude = op(self._magnitude, other._magnitude)
            units = self._units
        # If only self has a delta unit, other determines unit of result.
        elif self._get_delta_units() and not other._get_delta_units():
            magnitude = op(self._convert_magnitude_not_inplace(other._units), other._magnitude)
            units = other._units
        else:
            units = self._units
            magnitude = op(self._magnitude, _safe_other_magnitude(self._units))

    elif (
        op == operator.sub
        and len(self_non_mul_units) == 1
        and self._units[self_non_mul_unit] == 1
        and not other._has_compatible_delta(self_non_mul_unit)
    ):
        if self._units == other._units:
            magnitude = op(self._magnitude, other._magnitude)
        else:
            magnitude = op(self._magnitude, _safe_other_magnitude(self._units))
        units = self._units.rename(self_non_mul_unit, "delta_" + self_non_mul_unit)

    elif (
        op == operator.sub
        and len(other_non_mul_units) == 1
        and other._units[other_non_mul_unit] == 1
        and not self._has_compatible_delta(other_non_mul_unit)
    ):
        # we convert to self directly since it is multiplicative
        magnitude = op(self._magnitude, _safe_other_magnitude(self._units))
        units = self._units

    elif (
        len(self_non_mul_units) == 1
        # order of the dimension of offset unit == 1 ?
        and self._units[self_non_mul_unit] == 1
        and other._has_compatible_delta(self_non_mul_unit)
    ):
        # Replace offset unit in self by the corresponding delta unit.
        # This is done to prevent a shift by offset in the to()-call.
        tu = self._units.rename(self_non_mul_unit, "delta_" + self_non_mul_unit)
        magnitude = op(self._magnitude, _safe_other_magnitude(tu))
        units = self._units
    elif (
        len(other_non_mul_units) == 1
        # order of the dimension of offset unit == 1 ?
        and other._units[other_non_mul_unit] == 1
        and self._has_compatible_delta(other_non_mul_unit)
    ):
        # Replace offset unit in other by the corresponding delta unit.
        # This is done to prevent a shift by offset in the to()-call.
        tu = other._units.rename(other_non_mul_unit, "delta_" + other_non_mul_unit)
        magnitude = op(self._convert_magnitude_not_inplace(tu), other._magnitude)
        units = other._units
    else:
        raise OffsetUnitCalculusError(self._units, other._units)

    return self.__class__(magnitude, units)


# Override the _add_sub method
PlainQuantity._add_sub = _add_sub_modified


def __quantity__or__(self: PlainQuantity, unit: Unit):
    return self.to(unit)


PlainQuantity.__or__ = __quantity__or__


def __unit__ror__(self: Unit, value):
    if isinstance(value, PlainQuantity):
        return value.to(self)
    else:
        return self._REGISTRY.Quantity(value, self)


Unit.__ror__ = __unit__ror__


@contextmanager
def switch_units(value: bool):
    """Temporary switch support of units off or on and then revert it back to preivous value"""
    old_value = Settings.units_enabled

    try:
        Settings.units_enabled = value
        yield None
    finally:
        Settings.units_enabled = old_value
