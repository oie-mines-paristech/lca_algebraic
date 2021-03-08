import os
import sys
sys.path.insert(0, os.getcwd())
from lca_algebraic import *
import pytest

initDb('MyProject')

# Import Ecoinvent DB (if not already done)
# Update the name and path to the location of the ecoinvent database
# importDb("ecoinvent 3.4", './ecoinvent 3.4_cutoff_ecoSpold02/datasets')

# We use a separate DB for defining our foreground model / activities
# Choose any name
USER_DB = 'test-db'

def setup_function() :
    """Before each test"""

    print("resetting DB")
    SET_USER_DB(USER_DB)
    resetDb(USER_DB)
    resetParams(USER_DB)


def test_switch_activity_support_sevral_times_same_target() :
    """ Test that switch activity can target the same activity several times """

    # Enum param
    p1 = newEnumParam(
        'p1',
        values=["v1", "v2", "v3"],
        default="v1")

    alu = findTechAct("aluminium alloy production, AlMg3", loc="RER")

    act = newSwitchAct(USER_DB, "switchAct", p1, {
        "v1" : alu,
        "v2" : alu,
        "v3" : alu
    })

    climate = [m for m in bw.methods if 'ILCD 1.0.8 2016' in str(m) and 'no LT' in str(m)][1]

    res = multiLCAAlgebric(act, [climate], p1=["v1", "v2", "v3"])
    vals = res.values

    assert vals[0] == vals[1] and vals[1] == vals[2]

def test_enum_vales_are_enforced():

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





