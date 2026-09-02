import pytest
from numpy.ma.testutils import assert_array_equal

from lca_algebraic.params import _getAmountOrFormula, STORE_FORMULA_TAG
from lca_algebraic.units import unit_registry as u
from lca_algebraic.units import *
from lca_algebraic import (
    Settings,
    Min,
    Max,
    newFloatParam,
    newActivity,
    copyActivity,
    compute_impacts,
    interpolate_activities)
from test.conftest import USER_DB
from sympy import Float
import sympy

@pytest.fixture(scope="module", autouse=True)
def enable_units():
    # Enable units
    Settings.units_enabled = True

    yield

    # Disable units
    Settings.units_enabled = False


@pytest.fixture(scope="function", autouse=True)
def disable_autoscale():
    u.auto_scale = False
    yield


def test_dimensionless_units():
    """Added several alias to dimensionless units"""
    assert 2 * u.fraction + 1 * u.ratio + 1 * u.unit == 4 * u.dimensionless


def test_new_separate_units():
    define_separate_unit("kwp")

    assert 2 * u.kwp + 3 * u.kwp == 5 * u.kwp


def test_new_alias_unit():
    define_alias_unit("tkm", u.ton * u.km)
    define_alias_unit("biton", 1000 * u.ton)

    assert 1 * u.tkm == 1 * u.ton * u.km
    assert (1000 * u.ton + (1 * u.kton).to(u.ton)) == 2000.0 * u.ton


def test_auto_scale():
    u.auto_scale = True

    # Should scale autoamtically quantities of compatible units
    assert 1 * u.m + 1 * u.km == 1001 * u.meter


def test_newp_param():
    # When units are activated newXXXValue returns a parameter together with its unit (Quantity)
    p_with_unit = newFloatParam("p", 0, min=0, max=2, unit="kWh")

    p = p_with_unit.magnitude

    assert isinstance(p_with_unit, u.Quantity)
    assert p.with_unit() == p_with_unit
    assert p_with_unit.units == u.kWh


def test_new_activity_with_units(data):
    p1_meter = newFloatParam("p1", default=0, min=0, max=1, unit="m")
    p2_kg = newFloatParam("p2", default=0, min=0, max=1, unit="kg")
    p3_ton = newFloatParam("p3", default=0, min=0, max=1, unit="ton")

    unit_registry.auto_scale = True

    # Should fail : BG activities are all in kg
    with pytest.raises(DimensionalityError):
        act1 = newActivity(USER_DB, "act1", "kg", exchanges={data.bg_act1: 2.0 * p1_meter})

    # Should pass
    act1 = newActivity(USER_DB, "act1", "kg", exchanges={data.bg_act1: 2.0 * p2_kg})

    ex = act1.getExchange(input=data.bg_act1)
    assert _getAmountOrFormula(ex) == 2.0 * p2_kg.magnitude

    # Should convert ton to kg
    act1 = newActivity(USER_DB, "act1", "kg", exchanges={data.bg_act1: 2.0 * p3_ton})

    ex = act1.getExchange(input=data.bg_act1)
    assert _getAmountOrFormula(ex) == 2000.0 * p3_ton.magnitude

    # Should pass
    act1 = newActivity(USER_DB, "act1", "kg", exchanges={data.bg_act1: 2.0 | u.kg})

    ex = act1.getExchange(input=data.bg_act1)
    assert _getAmountOrFormula(ex) == 2.0

    # Should convert ton to kg
    act1 = newActivity(USER_DB, "act1", "kg", exchanges={data.bg_act1: 2.0 | u.ton})

    ex = act1.getExchange(input=data.bg_act1)
    assert _getAmountOrFormula(ex) == 2000.0

    unit_registry.auto_scale = False

    # Should fail (autoscale disabled)
    with pytest.raises(Exception) as e:
        act1 = newActivity(USER_DB, "act1", "kg", exchanges={data.bg_act1: 2 * p3_ton})
    assert "auto_scale" in str(e.value)


def test_update_activities_with_units(data):
    unit_registry.auto_scale = True

    p1_meter = newFloatParam("p1", default=0, min=0, max=1, unit="m")
    p2_ton = newFloatParam("p2", default=0, min=0, max=1, unit="ton")
    p2 = p2_ton.magnitude

    copyActivity(USER_DB, data.bg_act1)

    # Should fail : BG activities are all in kg
    with pytest.raises(DimensionalityError):
        act1 = newActivity(USER_DB, "act1", "kg", exchanges={data.bg_act1: 2 * p1_meter})

    # Should convert ton to kg
    act1 = newActivity(USER_DB, "act1", "kg", exchanges={data.bg_act1: 2 * p2_ton})

    ex = act1.getExchange(name="bg_act1")
    assert STORE_FORMULA_TAG in ex
    assert _getAmountOrFormula(ex) == 2000.0 * p2

    # Should convert ton to kg
    act1.updateExchanges({data.bg_act1: 1.0 | u.ton})
    act1.save()

    ex = act1.getExchange(name="bg_act1")
    assert STORE_FORMULA_TAG not in ex
    assert _getAmountOrFormula(ex) == 1000.0

    # Check if sympy formula are simplified to scala
    act1.updateExchanges({data.bg_act1: Float(1.0)*Float(10.0) | u.kg})
    act1.save()

    ex = act1.getExchange(name="bg_act1")
    assert STORE_FORMULA_TAG not in ex
    assert _getAmountOrFormula(ex) == 10.0


def test_interpolation_with_units(data):
    # Common helper to check results
    def check_impacts(model, p_values, expected_results):
        # Compute impacts for several values of p
        methods = list(expected_results)

        impacts = compute_impacts(model, methods, p=p_values)
        print(impacts)
        for i, k in enumerate(methods):
            values = impacts.iloc[:, i]
            assert_array_equal(values, expected_results[k])

    # The param is in meter
    p = newFloatParam("p", 1.0, min=1, max=3, unit="m")

    # Create act1 act2  having respectively 1.0, 2.0 kg of bio1, output unit : kg too
    act1 = newActivity(USER_DB, "act1", "kg", {data.bio1: 1.0 | u.kg})
    act2 = newActivity(USER_DB, "act2", "kg", {data.bio2: 1.0 | u.kg})

    # Interpolate including zero
    interp_with_zero = interpolate_activities(USER_DB, "interp_with_zero", p, {1.0 | u.m: act1, 3.0 | u.m: act2}, add_zero=True)

    check_impacts(
        interp_with_zero,
        p_values=[0.0, 1.0, 2.0, 3.0, 4.0],
        expected_results={
            data.ibio1: [0.0, 1.0, 0.5, 0.0, 0.0],
            data.ibio2: [0.0, 0.0, 0.5, 1.0, 1.0],
        },
    )


def test_interpolation_with_units(data):
    # Common helper to check results
    def check_impacts(model, p_values, expected_results):
        # Compute impacts for several values of p
        methods = list(expected_results)

        impacts = compute_impacts(model, methods, p=p_values)
        print(impacts)
        for i, k in enumerate(methods):
            values = impacts.iloc[:, i]
            assert_array_equal(values, expected_results[k])

    # The param is in meter
    p = newFloatParam("p", 1.0, min=1, max=3, unit="m")

    # Create act1 act2  having respectively 1.0, 2.0 kg of bio1, output unit : kg too
    act1 = newActivity(USER_DB, "act1", "kg", {data.bio1: 1.0 | u.kg})
    act2 = newActivity(USER_DB, "act2", "kg", {data.bio2: 1.0 | u.kg})

    # Interpolate including zero
    interp_with_zero = interpolate_activities(USER_DB, "interp_with_zero", p, {1.0 | u.m: act1, 3.0 | u.m: act2}, add_zero=True)

    check_impacts(
        interp_with_zero,
        p_values=[0.0, 1.0, 2.0, 3.0, 4.0],
        expected_results={
            data.ibio1: [0.0, 1.0, 0.5, 0.0, 0.0],
            data.ibio2: [0.0, 0.0, 0.5, 1.0, 1.0],
        },
    )


def test_compute_impact_with_functional_unit(data):
    # P1 in meter
    p1_m = newFloatParam("p1", default=1, min=0, max=2, unit="m")

    # P2 in Kg
    p2_kg = newFloatParam("p2", default=2, min=0, max=2, unit="kg")

    # Create activity with units
    act1 = newActivity(USER_DB, "act1", "kg", exchanges={data.bg_act1: 2 * p2_kg})

    functional_unit = 2 * p1_m

    # Ask with fonctional units of unit "meter
    res = compute_impacts(act1, data.ibio1, functional_unit=functional_unit)

    # Result should contain physical units in method names
    assert res.to_dict() == {"bio1 - total[MJ-Eq / meter]": {"act1": 2.0}}


def test_persist_load_params():
    """Custom units should be persisted to db and loaded correctly"""


def test_parse_db_unit():
    assert parse_db_unit("km-person") == u.km * u.person

def test_min_max_with_units():
    p0 = newFloatParam("p0", default=1, min=0, max=2, unit="m")
    p1 = newFloatParam("p1", default=1, min=0, max=2, unit="m")
    p2 = newFloatParam("p2", default=1, min=0, max=2, unit="m")
    p3 = newFloatParam("p3", default=1, min=0, max=2, unit="kg")

    emin = Min(p0, p1, p2)
    assert isinstance(emin, u.Quantity)
    assert emin.magnitude.func == sympy.Min
    assert emin.units == u.m

    emax = Max(p0, p1, p2)
    assert isinstance(emax, u.Quantity)
    assert emax.magnitude.func == sympy.Max
    assert emax.units == u.m

    with pytest.raises(Exception):
        emin = Min(p0, p1, p3)

    with pytest.raises(Exception):
        emax = Max(p0, p1, p3)
