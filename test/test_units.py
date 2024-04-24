import pytest

from lca_algebraic.params import _getAmountOrFormula
from lca_algebraic.units import unit_registry as u
from lca_algebraic.units import *
from lca_algebraic import Settings, newFloatParam, newActivity, copyActivity, compute_impacts
from test.conftest import USER_DB


@pytest.fixture(scope="module", autouse=True)
def enable_units():
    # Enable units
    Settings.units_enabled = True

    yield

    # Enable units
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


def test_add_exchanges(data):
    p1_meter = newFloatParam("p1", default=0, min=0, max=1, unit="m")
    p2_kg = newFloatParam("p2", default=0, min=0, max=1, unit="kg")
    p3_ton = newFloatParam("p3", default=0, min=0, max=1, unit="ton")

    unit_registry.auto_scale = True

    # Should fail : BG activities are all in kg
    with pytest.raises(DimensionalityError):
        act1 = newActivity(USER_DB, "act1", "kg", exchanges={data.bg_act1: 2 * p1_meter})

    # Should pass
    act1 = newActivity(USER_DB, "act1", "kg", exchanges={data.bg_act1: 2 * p2_kg})

    unit_registry.auto_scale = False

    # Should fail (autoscale disabled)
    with pytest.raises(Exception) as e:
        act1 = newActivity(USER_DB, "act1", "kg", exchanges={data.bg_act1: 2 * p3_ton})
    assert "auto_scale" in str(e.value)


def test_update_exchanges(data):
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

    assert _getAmountOrFormula(act1.getExchange(name="bg_act1")) == 2000.0 * p2


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
