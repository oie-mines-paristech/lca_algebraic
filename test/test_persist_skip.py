"""TDD: smart parameter persistence.

Re-creating a parameter that is *identical* to the one already persisted in
Brightway must be a no-op: rewriting it triggers a full parameter
re-computation (`recalculate`) and defeats the purpose of the parameter cache.
This is what makes repeated `new*Param` calls (e.g. rebuilding descriptors on
every run) cheap.
"""
import pytest

import lca_algebraic.params as P
from lca_algebraic import DistributionType, loadParams, newBoolParam, newEnumParam, newFloatParam, resetParams


@pytest.fixture
def writes(monkeypatch):
    """Records every call to `new_project_parameters` (a DB write)."""
    calls = []
    real = P.new_project_parameters

    def spy(out, *args, **kwargs):
        calls.append(out)
        return real(out, *args, **kwargs)

    monkeypatch.setattr(P, "new_project_parameters", spy)
    return calls


def test_first_creation_persists(writes):
    newFloatParam("p", 3.0, min=0.0, max=5.0, unit="kg", group="g")
    assert len(writes) == 1


def test_repeated_identical_creation_is_not_repersisted(writes):
    newFloatParam("p", 3.0, min=0.0, max=5.0, unit="kg", group="g")
    newFloatParam("p", 3.0, min=0.0, max=5.0, unit="kg", group="g")
    assert len(writes) == 1


def test_identical_recreation_after_loadparams_is_not_repersisted(writes):
    newFloatParam("p", 3.0, min=0.0, max=5.0, unit="kg", group="g")
    loadParams()
    newFloatParam("p", 3.0, min=0.0, max=5.0, unit="kg", group="g")
    assert len(writes) == 1


def test_int_default_matches_float_amount(writes):
    """A default given as int (0) must match a stored amount (0) too."""
    newFloatParam("p", 0, min=0, max=5, unit="kg", group="g")
    loadParams()
    newFloatParam("p", 0, min=0, max=5, unit="kg", group="g")
    assert len(writes) == 1


def test_changed_default_is_repersisted(writes):
    newFloatParam("p", 3.0, min=0.0, max=5.0, unit="kg", group="g")
    loadParams()
    newFloatParam("p", 4.0, min=0.0, max=5.0, unit="kg", group="g")
    assert len(writes) == 2


def test_changed_bound_is_repersisted(writes):
    newFloatParam("p", 3.0, min=0.0, max=5.0, unit="kg", group="g")
    loadParams()
    newFloatParam("p", 3.0, min=1.0, max=5.0, unit="kg", group="g")
    assert len(writes) == 2


def test_changed_unit_is_repersisted(writes):
    newFloatParam("p", 3.0, min=0.0, max=5.0, unit="kg", group="g")
    loadParams()
    newFloatParam("p", 3.0, min=0.0, max=5.0, unit="t", group="g")
    assert len(writes) == 2


def test_changed_group_is_repersisted(writes):
    newFloatParam("p", 3.0, min=0.0, max=5.0, unit="kg", group="g")
    loadParams()
    newFloatParam("p", 3.0, min=0.0, max=5.0, unit="kg", group="other")
    assert len(writes) == 2


def test_bool_identical_roundtrip(writes):
    newBoolParam("b", 1, group="g")
    loadParams()
    newBoolParam("b", 1, group="g")
    assert len(writes) == 1


def test_bool_changed_is_repersisted(writes):
    newBoolParam("b", 1, group="g")
    loadParams()
    newBoolParam("b", 0, group="g")
    assert len(writes) == 2


def test_enum_identical_roundtrip(writes):
    newEnumParam("e", default="a", values=["a", "b"], group="g")
    loadParams()
    newEnumParam("e", default="a", values=["a", "b"], group="g")
    assert len(writes) == 1


def test_enum_changed_default_is_repersisted(writes):
    newEnumParam("e", default="a", values=["a", "b"], group="g")
    loadParams()
    newEnumParam("e", default="b", values=["a", "b"], group="g")
    assert len(writes) == 2


def test_enum_changed_weights_is_repersisted(writes):
    newEnumParam("e", default="a", values={"a": 1, "b": 1}, group="g")
    loadParams()
    newEnumParam("e", default="a", values={"a": 1, "b": 2}, group="g")
    assert len(writes) == 2


def test_reset_params_forgets_signatures(writes):
    """After resetParams the DB is empty, so a re-creation must persist again."""
    newFloatParam("p", 3.0, min=0.0, max=5.0, unit="kg", group="g")
    resetParams()
    newFloatParam("p", 3.0, min=0.0, max=5.0, unit="kg", group="g")
    assert len(writes) == 2
