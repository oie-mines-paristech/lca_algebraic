"""Cache pitfalls for axis MatricesResult variants (AxisMatricesSet)."""

from unittest.mock import patch

import pytest

from lca_algebraic import compute_impacts, newActivity, newFloatParam
from lca_algebraic.cache import ExprCache, clear_caches
from lca_algebraic.lca import (
    AxisMatricesSet,
    _computeAxisMatricesSet,
    _expr_base_key,
    _matrices_set_cache_key,
    _model_cache_key,
)
from lca_algebraic.settings import temp_settings
from test.conftest import USER_DB


def _minimal_axis_model(data):
    """Small 3-branch model used in test_axis."""
    act_internal0 = newActivity(USER_DB, "act_internal0", "unit", {data.bio1: 1.0})
    act1 = newActivity(USER_DB, "act1", "unit", {act_internal0: 1.0}, phase="phase a")
    act2 = newActivity(USER_DB, "act2", "unit", {act_internal0: 2.0}, phase="phase b")
    act3 = newActivity(USER_DB, "act3", "unit", {act_internal0: 3.0})
    return newActivity(USER_DB, "model", "unit", {act1: 1, act2: 1, act3: 1})


def _axis_expected():
    return {
        "phase_a": 2.0,
        "phase_b": 4.0,
        "*other*": 6.0,
        "*all*": 12.0,
    }


def test_axis_matrices_set_isolated_from_no_axis_cache(data):
    """Pitfall: axis and no-axis must not share the same matrices_set entry."""
    model = _minimal_axis_model(data)

    compute_impacts(model, [data.ibio1], axis="phase")
    compute_impacts(model, [data.ibio1])

    with ExprCache(USER_DB) as cache:
        axis_key = _matrices_set_cache_key(_expr_base_key("phase"))
        no_axis_key = _matrices_set_cache_key(_expr_base_key(None))

        assert axis_key in cache.data
        assert no_axis_key in cache.data
        assert cache.data[axis_key] is not cache.data[no_axis_key]

        axis_set: AxisMatricesSet = cache.data[axis_key]
        no_axis_set: AxisMatricesSet = cache.data[no_axis_key]

        assert set(axis_set.variants.keys()) >= {"_all_", "phase_a", "phase_b"}
        assert set(no_axis_set.variants.keys()) == {"_all_"}


def test_axis_matrices_set_reused_across_models(data):
    """Pitfall: matrices_set is per DB/axis, model cache is per model."""
    model_a = _minimal_axis_model(data)
    model_b = newActivity(
        USER_DB,
        "model_b",
        "unit",
        exchanges={data.bg_act1: 1.0},
        phase="phase a",
    )

    calls = []

    def counting_compute(db_name, axis_attr):
        calls.append((db_name, axis_attr))
        return _computeAxisMatricesSet(db_name, axis_attr)

    with patch("lca_algebraic.lca._computeAxisMatricesSet", side_effect=counting_compute):
        compute_impacts(model_a, [data.ibio1], axis="phase")
        compute_impacts(model_b, [data.ibio1], axis="phase")

    assert len(calls) == 1

    with ExprCache(USER_DB) as cache:
        mset_key = _matrices_set_cache_key(_expr_base_key("phase"))
        model_keys = [k for k in cache.data if isinstance(k, tuple) and len(k) == 4 and k[1] == "model"]

        assert mset_key in cache.data
        assert len(model_keys) == 2
        assert {k[2] for k in model_keys} == {str(model_a), str(model_b)}


def test_model_cache_key_includes_alpha(data):
    """Pitfall: model expressions must not be reused across different alphas."""
    model = newActivity(USER_DB, "m_alpha", "unit", {data.bio1: 1.0})

    r1 = compute_impacts({model: 1.0}, [data.ibio1])
    r2 = compute_impacts({model: 2.0}, [data.ibio1])

    assert r1.values[0] == 1.0
    assert r2.values[0] == 2.0

    with ExprCache(USER_DB) as cache:
        base_key = _expr_base_key(None)
        key_alpha_1 = _model_cache_key(base_key, model, 1.0)
        key_alpha_2 = _model_cache_key(base_key, model, 2.0)

        assert key_alpha_1 in cache.data
        assert key_alpha_2 in cache.data
        assert cache.data[key_alpha_1] is not cache.data[key_alpha_2]


def test_factorized_cache_namespace_is_separate(data):
    """Pitfall: factorize_static_bg toggles a different cache namespace."""
    model = _minimal_axis_model(data)

    with temp_settings(factorize_static_bg=True):
        compute_impacts(model, [data.ibio1], axis="phase")
        with ExprCache(USER_DB) as cache:
            assert (("phase", "factorized"), "matrices_set") in cache.data

    clear_caches()

    with temp_settings(factorize_static_bg=False):
        compute_impacts(model, [data.ibio1], axis="phase")
        with ExprCache(USER_DB) as cache:
            assert (("phase",), "matrices_set") in cache.data
            assert (("phase", "factorized"), "matrices_set") not in cache.data


def test_axis_variants_share_background_activities(data):
    """Pitfall: all variants must expose the same bg_acts for impacts[] indexing."""
    matrices_set = _computeAxisMatricesSet(USER_DB, axis_attr="phase")
    matrices_set.ensure_bg_acts_consistency()

    ref = matrices_set.variants["_all_"].bg_acts
    for tag, matrices in matrices_set.variants.items():
        assert matrices.bg_acts == ref


def test_warm_axis_cache_still_returns_correct_results(data):
    """Regression: cache hits must preserve axis arithmetic (_all_ - tag)."""
    p1 = newFloatParam("p1", 2, min=1, max=3)
    model = _minimal_axis_model(data)

    # Warm caches (same setup as test_axis)
    compute_impacts(model, [data.ibio1], functional_unit=p1, axis="phase", p1=0.5)
    res = compute_impacts(model, [data.ibio1], functional_unit=p1, axis="phase", p1=0.5)
    got = {key: val for key, val in zip(res.index.values, res[res.columns[0]].values)}

    assert got == _axis_expected()
