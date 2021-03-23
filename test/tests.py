import os
import sys
import pytest

sys.path.insert(0, os.getcwd())
sys.path.insert(0, os.path.join(os.getcwd(), "test"))

from lca_algebraic import *
from fixtures import *

USER_DB = "fg"
BG_DB= "bg"
METHOD_PREFIX='tests'

SET_USER_DB(USER_DB)
SET_BG_DB(BG_DB)

# Reset func project, empty DB
initDb('tests')
init_acts(BG_DB)
init_methods(BG_DB, METHOD_PREFIX)

# Import Ecoinvent DB (if not already done)
# Update the name and path to the location of the ecoinvent database
# importDb("ecoinvent 3.4", './ecoinvent 3.4_cutoff_ecoSpold02/datasets')

# We use a separate DB for defining our foreground model / activities
# Choose any name


def setup_function() :
    """Before each test"""

    print("resetting fg DB")

    resetDb(USER_DB)
    resetParams(USER_DB)


def test_switch_activity_support_sevral_times_same_target() :
    """ Test that switch activity can target the same activity several times """

    # Enum param
    p1 = newEnumParam(
        'p1',
        values=["v1", "v2", "v3"],
        default="v1")

    bg_act1 = findTechAct("bg_act1")

    act = newSwitchAct(USER_DB, "switchAct", p1, {
        "v1" : bg_act1,
        "v2" : bg_act1,
        "v3" : bg_act1
    })

    impact = (METHOD_PREFIX, 'all', 'total')

    res = multiLCAAlgebric(act, [impact], p1=["v1", "v2", "v3"])
    vals = res.values

    assert vals[0] == vals[1] and vals[1] == vals[2]

def test_new_switch_act_with_tuples() :

    p1 = newEnumParam(
        'p1',
        values=["v1", "v2", "v3"],
        default="v1")

    bg_act1 = findTechAct("bg_act1")
    bg_act2 = findTechAct("bg_act2")

    newSwitchAct(USER_DB, "switchAct", p1, {
        "v1": (bg_act1, 2.0),
        "v2": (bg_act2, 3.0),
    })



def test_enum_values_are_enforced():

    # Enum param
    p1 = newEnumParam(
        'p1',
        values=["v1", "v2", "v3"], default="v1")

    act = newActivity(USER_DB, "Foo", "unit")

    climate = [m for m in bw.methods if 'ILCD 1.0.8 2016' in str(m) and 'no LT' in str(m)][1]

    with pytest.raises(Exception) as exc:
        multiLCAAlgebric(act, climate, p1="bar")

    assert 'Invalid value' in str(exc)

if __name__ == '__main__':
    pytest.main(sys.argv)





