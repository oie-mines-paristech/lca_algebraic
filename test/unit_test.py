import os
import sys
from tempfile import mkstemp

sys.path.insert(0, os.getcwd())
sys.path.insert(0, os.path.join(os.getcwd(), "test"))

import pytest
from conftest import BG_DB, METHOD_PREFIX, USER_DB
from fixtures import *
from numpy.testing import assert_array_equal

from lca_algebraic.database import _isForeground, setForeground, setBackground
from lca_algebraic.params import _param_registry
from lca_algebraic.lca import _cachedActToExpression
from pandas.testing import assert_frame_equal


def test_load_params():
    _p1 = newEnumParam("p1", values={"v1": 0.6, "v2": 0.3}, default="v1")
    _p2 = newFloatParam("p2", min=1, max=3, default=2, distrib=DistributionType.TRIANGLE)
    _p3 = newBoolParam("p3", default=1, formula=_p2 + 4)

    _param_registry().clear()

    # Params are loaded as global variable with their real names
    loadParams()

    # Get params from in memory DB
    loaded_params = {(param.name, param.dbname): param for param in _param_registry().all()}

    assert _p1.__dict__ == loaded_params[("p1", None)].__dict__
    assert _p2.__dict__ == loaded_params[("p2", None)].__dict__
    assert _p3.__dict__ == loaded_params[("p3", None)].__dict__


def test_export(data):
    p1 = newFloatParam("p1", default=0.5, distrib=DistributionType.FIXED)
    p3_fg = newBoolParam("p3", default=1)  # Param with same name linked to a user DB

    act1 = newActivity(
        USER_DB,
        "act1",
        "unit",
        {
            data.bio1: 2 * p1,
            data.bio2: 3 * p3_fg,
        },
    )

    f, filename = mkstemp()

    outfile = export_db(USER_DB, filename)

    # Clear all
    resetDb(USER_DB)
    resetParams()
    _param_registry().clear()

    import_db(filename)

    print(filename)

    # Check params are the same
    p1_ = _param_registry()["p1"]
    p3_ = _param_registry()["p3"]

    assert p1.__dict__ == p1_.__dict__
    assert p3_fg.__dict__ == p3_.__dict__

    # Test multLCA
    act1 = findActivity("act1", db_name=USER_DB)
    res = compute_impacts(act1, [data.imulti])

    assert res.values[0] == 7.0


def test_switch_activity_support_sevral_times_same_target():
    """Test that switch activity can target the same activity several times"""

    # Enum param
    p1 = newEnumParam("p1", values=["v1", "v2", "v3"], default="v1")

    bg_act1 = findTechAct("bg_act1")

    act = newSwitchAct(USER_DB, "switchAct", p1, {"v1": bg_act1, "v2": bg_act1, "v3": bg_act1})

    impact = (METHOD_PREFIX, "all", "total")

    res = compute_impacts(act, [impact], p1=["v1", "v2", "v3"])
    vals = res.values

    assert vals[0] == vals[1] and vals[1] == vals[2]


def test_new_switch_act_with_tuples():
    p1 = newEnumParam("p1", values=["v1", "v2", "v3"], default="v1")

    bg_act1 = findTechAct("bg_act1")
    bg_act2 = findTechAct("bg_act2")

    newSwitchAct(
        USER_DB,
        "switchAct",
        p1,
        {
            "v1": (bg_act1, 2.0),
            "v2": (bg_act2, 3.0),
        },
    )


def test_list_params_should_support_missing_groups():
    p1 = newFloatParam("p1", default=1.0, distrib=DistributionType.FIXED)
    p2 = newFloatParam("p2", default=2.0, distrib=DistributionType.FIXED, group="mygroup")

    list_parameters()


def test_freeze(data):
    p1 = newFloatParam("p1", default=1.0, distrib=DistributionType.FIXED)
    p2 = newFloatParam("p2", default=1.0, distrib=DistributionType.FIXED)

    newActivity(
        USER_DB,
        "act1",
        "unit",
        {
            data.bio1: 2 * p1,
            data.bio2: 3 * p2,
        },
    )

    # p1 should be set as default value 1
    freezeParams(USER_DB, p2=2)

    # Load back activity
    act1 = findActivity("act1", db_name=USER_DB)

    for exc in act1.exchanges():
        # Don't show production
        if exc["type"] == "production":
            continue

        name = exc["name"]
        amount = exc["amount"]

        # Brightway2 does not like ints ...
        assert isinstance(amount, float)

        if name == "bio1":
            # p1=1 (default) * 2
            assert amount == 2.0
        elif name == "bio2":
            # p2=2 * 3
            assert amount == 6.0


def test_find_activities():
    act1 = newActivity(USER_DB, "Activity 1", "unit", categories=["cat1", "cat2"])

    res = findActivity("Activity 1", db_name=USER_DB, case_sensitive=True)
    assert res == act1

    # Case insensitive
    assert findActivity("activity*", db_name=USER_DB) == act1

    # Test categories
    assert findActivity("activity*", categories=["cat1", "cat2"], db_name=USER_DB) == act1
    assert findActivity("activity*", categories=["cat1"], db_name=USER_DB, single=False) == []
    assert findActivity("activity*", category="cat1", db_name=USER_DB) == act1


def test_enum_values_are_enforced():
    # Enum param
    p1 = newEnumParam("p1", values=["v1", "v2", "v3"], default="v1")

    act = newActivity(USER_DB, "Foo", "unit")

    climate = [m for m in bw.methods if "ILCD 1.0.8 2016" in str(m) and "no LT" in str(m)][1]

    with pytest.raises(Exception) as exc:
        compute_impacts(act, climate, p1="bar")

    assert "Invalid value" in str(exc)


def test_setforeground():
    setForeground(USER_DB)

    assert _isForeground(USER_DB)

    setBackground(USER_DB)

    assert _isForeground(USER_DB) == False


def test_reset_params():
    # Define 3 variables with same name, attached to project or db (user or bg)
    newBoolParam("p1", False)
    newBoolParam("p2", False)
    newBoolParam("p3", False)  # Project param

    # Should delete all
    resetParams()
    assert len(_param_registry().all()) == 0

    # Check they are deleted in Db as well
    loadParams()
    assert len(_param_registry().all()) == 0


def test_simplify_model(data):
    # key param, varying from 1 to 2
    p1 = newFloatParam("p1", 1, min=1, max=2)

    # Minor param with constant value 0.001
    p2 = newFloatParam("p2", 1, min=0.001, max=0.001)

    # Model p1 + 0.001*p1 + p2
    m1 = newActivity(USER_DB, "m1", "kg", {data.bio1: p1 * (p1 + 0.001 * p1 + p2)})

    # Simplified model, removing minor sum terms (default)
    res = sobol_simplify_model(m1, [data.ibio1], simple_products=False)[0]
    assert res.expr.__repr__() == "1.0*p1**2"

    # Simplified model, without removing minor sum term
    res = sobol_simplify_model(m1, [data.ibio1], simple_sums=False, simple_products=False)[0]
    assert res.expr.__repr__() == "p1*(1.0*p1 + 0.001)"

    # -- Simplify products

    p3 = newFloatParam("p3", 1, min=1, max=1.001)
    p4 = newFloatParam("p4", 1, min=-0.999, max=-0.998)

    # Boolean should not be removed
    p5 = newBoolParam("p5", 1)
    m2 = newActivity(USER_DB, "m2", "kg", {data.bio1: 4.0 + 5 * p3 + 3 * p4 + 3 * p5})

    res = sobol_simplify_model(m2, [data.ibio1], simple_products=True)[0]
    assert res.expr.__repr__() == "3.0*p5 + 6.01"


def test_db_params_lca(data):
    """Test multiLCAAlgebraic with parameters with similar names from similar DBs"""
    USER_DB2 = "fg2"
    resetDb(USER_DB2)

    # Define 3 variables with same name, attached to project or db (user or bg)
    p1_project = newFloatParam(
        "p1",
        0,
        min=0,
        max=2,
    )
    p1_user = newFloatParam("p1", 1, min=0, max=2, dbname=USER_DB)
    p1_user2 = newFloatParam("p1", 2, min=0, max=2, dbname=USER_DB2)

    # Create 2 models : one for each user db, using different params with same name
    m1 = newActivity(USER_DB, "m1", "kg", {data.bio1: 2.0 * p1_user})
    m2 = newActivity(USER_DB2, "m2", "kg", {data.bio1: 2.0 * p1_user2})

    # p1 as default value of 1 for user db 1
    res = compute_impacts(m1, [data.ibio1])
    assert res.values[0] == 2.0

    # p1 as default value of 2 for user db 2
    res = compute_impacts(m2, [data.ibio1])
    assert res.values[0] == 4.0

    # Overriding p1 for m1
    res = compute_impacts(m1, [data.ibio1], p1=3)
    assert res.values[0] == 6.0

    # Overriding p1 for m2
    res = compute_impacts(m2, [data.ibio1], p1=4)
    assert res.values[0] == 8.0


def test_params_with_formulas(data):
    p1 = newFloatParam("p1", default=1, min=0, max=2)
    p2 = newFloatParam("p2", default=0, min=0, max=2, formula=2 * p1)
    p3 = newFloatParam("p3", default=0, min=0, max=3, formula=2 * p2)

    # act1 = bio1 : p2
    act1 = newActivity(USER_DB, "act1", "kg", {data.bio1: p2})

    # act2 = bio1 : p2
    act2 = newActivity(USER_DB, "act2", "kg", {data.bio1: p3})

    # By default, p2 should be (p1=1)*2
    res = compute_impacts(act1, [data.ibio1])
    assert res.values[0] == 2.0

    # p2 is computed automatically from p1
    res = compute_impacts(act1, [data.ibio1], p1=2)
    assert res.values[0] == 4.0

    # Overriding p2 should disable computing of formula
    res = compute_impacts(act1, [data.ibio1], p1=2, p2=1)

    assert res.values[0] == 1.0

    # p1 > p2 > p3
    res = compute_impacts(act2, [data.ibio1], p1=4)

    assert res.values[0] == 16.0

    # Should also work with lists of values
    res = compute_impacts(act1, [data.ibio1], p1=[1.0, 2.0])

    res.iloc[:, 0] == [2.0, 4.0]


def test_compute_impacts_return_params(data):
    p1 = newFloatParam("p1", default=1, min=0, max=2)
    p_default = newFloatParam("p_default", default=1, min=0, max=2)
    p_computed = newFloatParam("p_computed", default=0, min=0, max=3, formula=2 * p1)
    p_unused = newFloatParam("p_unused", default=0, min=0, max=2)

    act1 = newActivity(
        USER_DB,
        "act1",
        "kg",
        {data.bio1: p1, data.bio2: p_default, data.bio3: p_computed},
    )

    res: TabbedDataframe = compute_impacts(
        act1,
        [data.ibio1, data.ibio2, data.ibio3],
        functional_unit=p_default,
        return_params=True,
        description="Something here",
        # Params
        p1=[1, 2, 3],
    )

    assert res.dataframes["Results"].to_dict() == {
        "bio1 - total[MJ-Eq]": {1: 1.0, 2: 2.0, 3: 3.0},
        "bio2 - total[MJ-Eq]": {1: 1.0, 2: 1.0, 3: 1.0},
        "bio3 - total[MJ-Eq]": {1: 2.0, 2: 4.0, 3: 6.0},
    }

    assert res.dataframes["Parameters"].to_dict() == {
        "min": {("", "p1"): 0, ("", "p_computed"): 0, ("", "p_default"): 0},
        "max": {("", "p1"): 2, ("", "p_computed"): 3, ("", "p_default"): 2},
        "default": {("", "p1"): 1, ("", "p_computed"): 0, ("", "p_default"): 1},
        "value_1": {("", "p1"): 1.0, ("", "p_computed"): 2.0, ("", "p_default"): 1.0},
        "value_2": {("", "p1"): 2.0, ("", "p_computed"): 4.0, ("", "p_default"): 1.0},
        "value_3": {("", "p1"): 3.0, ("", "p_computed"): 6.0, ("", "p_default"): 1.0},
    }

    # res.to_excel("tst.xlsx")


def test_switch_value(data):
    p1 = newFloatParam("p1", default=0, min=0, max=2)
    p2 = newFloatParam("p2", default=0, min=0, max=2)
    switch_param = newEnumParam("switch_param", default="p1", values=["from_p1", "from_p2"])
    computed = newFloatParam("computed", default=0, min=0, max=3)

    computed.formula = switchValue(switch_param, from_p1=p1 * 2, from_p2=p2 * 3)

    act1 = newActivity(USER_DB, "act1", "kg", {data.bio1: computed})

    res = compute_impacts(act1, [data.ibio1], switch_param="from_p1", p1=1)

    assert res.values[0] == 2.0

    res = compute_impacts(act1, [data.ibio1], switch_param="from_p2", p2=1)

    assert res.values[0] == 3.0


def test_interpolation(data):
    # Common helper to check results
    def check_impacts(model, p_values, expected_results):
        # Compute impacts for several values of p
        impacts = compute_impacts(model, [data.ibio1], p=p_values)

        values = impacts[impacts.columns[0]]
        assert_array_equal(values, expected_results)

    # Define param
    p = newFloatParam("p", 1.0, min=1, max=3)

    # Create act1 act2 and act4 having respectively 1.0, 2.0 units of bio1
    act1, act2 = [newActivity(USER_DB, "act%d" % v, "unit", {data.bio1: v}) for v in [1.0, 2.0]]

    # Interpolate between 1 : act1 (1 bio1) and 3 : act2 (2 bio1)
    interp1 = interpolate_activities(USER_DB, "interp1", p, {1.0: act1, 3.0: act2})

    check_impacts(interp1, [0.0, 1.0, 2.0, 3.0, 5.0], [0.5, 1.0, 1.5, 2.0, 3.0])

    # Interpolate including zero
    interp_with_zero = interpolate_activities(USER_DB, "interp_w_zero", p, {1.0: act1, 3.0: act2}, add_zero=True)

    check_impacts(interp_with_zero, [0.0, 0.5, 1.0, 3.0], [0.0, 0.5, 1.0, 2.0])


def test_axis(data):
    p1 = newFloatParam("p1", 2, min=1, max=3)

    act1_phase_a = newActivity(USER_DB, "act1", "unit", {data.bio1: 1.0}, phase="phase a")

    act2_phase_b = newActivity(USER_DB, "act2", "unit", {data.bio1: 2.0}, phase="phase b")

    act3_no_phase = newActivity(USER_DB, "act3", "unit", {data.bio1: 3.0})

    model = newActivity(
        USER_DB,
        "model",
        "unit",
        {
            act1_phase_a: 1,
            act2_phase_b: 1,
            act3_no_phase: 1,
        },
    )

    res = compute_impacts(model, [data.ibio1], functional_unit=p1, axis="phase", p1=0.5)

    # Compute twice to warm the cache :
    # Creating dummy activities clears it and we cannot check it works properly otherwize
    res = compute_impacts(model, [data.ibio1], functional_unit=p1, axis="phase", p1=0.5)

    res = {key: val for key, val in zip(res.index.values, res[res.columns[0]].values)}

    expected = dict(phase_a=2.0, phase_b=4.0, _other_=6.0)
    expected["*sum*"] = 12.0

    assert res == expected

    res = compute_impacts(model, [data.ibio1], functional_unit=p1, p1=0.5)

    print(res)

    assert res.values[0] == 12.0


def test_compute_impacts_with_parametrized_fu(data):
    p1 = newFloatParam("p1", 2, min=1, max=3)

    m1 = newActivity(USER_DB, "m1", "kg", {data.bio1: 4 * p1})

    # Functional unit is parametrized
    fu_value = 2 * p1

    # Compute single value
    res = compute_impacts(m1, [data.ibio1], functional_unit=fu_value, p1=[1.0])

    assert res.iloc[0, 0] == 2.0

    # Compute list of values
    res = compute_impacts(m1, [data.ibio1], functional_unit=fu_value, p1=[1.0, 2.0])

    print(res)

    assert res.iloc[0, 0] == 2.0
    assert res.iloc[1, 0] == 2.0


def test_compute_inventory(data):
    p1 = newFloatParam("p1", 1, min=1, max=3)
    p2 = newFloatParam("p2", 1, min=1, max=3)

    # Two nested activities
    fg_act1 = newActivity(USER_DB, name="act1", unit="kg", exchanges={data.bg_act1: p1})

    root_act = newActivity(
        USER_DB, name="root_act", unit="kg", exchanges={fg_act1: 1, data.bg_act1: 1, data.bio1: 1, data.bg_act2: p2}
    )

    df: DataFrame = compute_inventory(root_act, functional_unit=10, p2=3)

    df_expected = DataFrame(
        [
            {
                "database": "bg",
                "name": "bg_act1",
                "location": "GLO",
                "unit": "kg",
                "value": 0.2,
            },
            {
                "database": "bg",
                "name": "bg_act2",
                "location": "GLO",
                "unit": "kg",
                "value": 0.3,
            },
            {"database": "bg", "name": "bio1", "location": "GLO", "unit": "kg", "value": 0.1},
        ]
    )

    assert_frame_equal(df_expected, df, rtol=1e-03)


def test_inventory_loops_should_work(data):
    act1 = newActivity(USER_DB, "act1", "kg")
    act1.addExchanges({act1: 0.5, data.bio1: 1})  # Loop on itself  # Bg act

    res = compute_impacts(act1, data.ibio1)

    assert res.values[0] == 12.0


def test_should_list_params_with_mixed_groups(data):
    """Test for bug #13 : https://github.com/oie-mines-paristech/lca_algebraic/issues/13"""
    p1 = newFloatParam("foo", 2, min=1, max=3, group="foo")
    bar = newFloatParam("bar", 2, min=1, max=3)

    m1 = newActivity(USER_DB, "m1", "kg", {data.bio1: 2.0 * p1 + bar})

    oat_matrix(m1, [data.ibio1, data.ibio2])


def test_oat_should_work_with_named_params(data):
    """Test for bug #12 : https://github.com/oie-mines-paristech/lca_algebraic/issues/12"""
    p1 = newFloatParam("foo", 2, min=1, max=3, group="foo")
    bar = newFloatParam("bar", 2, min=1, max=3)

    m1 = newActivity(USER_DB, "m1", "kg", {data.bio1: 2.0 * p1 + bar})

    oat_matrix(model=m1, impacts=[data.ibio1, data.ibio2])


def test_multiLCAAlgebric_with_dict(data):
    """Tests parameters can be used in 'power'"""

    m1 = newActivity(USER_DB, "m1", "kg", {data.bio1: 1})

    m2 = newActivity(USER_DB, "m2", "kg", {data.bio2: 1})

    res = compute_impacts({m1: 1, m2: 2}, [data.ibio1, data.ibio2])

    assert res.iloc[0, 0] == 1.0
    assert res.iloc[1, 1] == 2.0
    assert res.iloc[0, 1] == 0.0
    assert res.iloc[1, 0] == 0.0


def test_params_as_power(data):
    """Tests parameters can be used in 'power'"""

    p1 = newFloatParam("p1", 2, min=0, max=2)

    m1 = newActivity(USER_DB, "m1", "kg", {data.bio1: 2.0**p1})

    res = compute_impacts(m1, [data.ibio1], p1=2)
    assert res.values[0] == 4.0


def test_brightway_lca(data):
    """Tests parameters can be used in 'power'"""

    p1 = newFloatParam("p1", 2, min=0, max=2)

    act1 = newActivity(USER_DB, "act1", "kg", {data.bio1: 2.0 * p1})
    act2 = newActivity(USER_DB, "act2", "kg", {act1: 1})

    res = multiLCA(act2, [data.ibio1], p1=2)
    assert res.values[0] == 4.0


def test_named_parameters_for_with_db_context(data):
    """Tests functions annotated with with_context_db, still support named db .
    See: https://github.com/oie-mines-paristech/lca_algebraic/issues/12
    """
    m1 = newActivity(USER_DB, "m1", "kg", {data.bio1: 1})


if __name__ == "__main__":
    pytest.main(sys.argv)
