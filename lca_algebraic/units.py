import brightway2 as bw
from pint import DimensionalityError, UnitRegistry

from .base_utils import getActByCode

NEW_UNITS = {"person"}

UNIT_ALIASES = {
    "ratio": "count",
    "fraction": "count",
    "unit": "count",
    "ton": "metric_ton",  # Override silly US ton to metric ton (f*ck imperial unit system ...)
}


# The global Unit registry
unit_registry = UnitRegistry()


def check_unit_consistency(db_name: str):
    """Check units of exchanges VS units of target activities in a single database"""
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


# Add some additional units

for unit in NEW_UNITS:
    define_separate_unit(unit)

for key, val in UNIT_ALIASES.items():
    unit_registry.define(f"@alias {val} = {key}")
