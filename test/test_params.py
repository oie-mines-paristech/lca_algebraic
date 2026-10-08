import math

import numpy as np
import pytest
from bw2data.parameters import ProjectParameter
from stats_arrays import LognormalUncertainty, UncertaintyBase

from lca_algebraic import DistributionType, loadParams, newFloatParam
from lca_algebraic.params import _param_registry

# ppf(0) and ppf(1) give the bounds of the support
BOUNDS = np.array([0.0, 1.0])


def normal(min=None, max=None):
    return newFloatParam("p", 10, min=min, max=max, std=2, distrib=DistributionType.NORMAL, save=False)


def test_normal_truncated_at_zero():
    p = newFloatParam("p", 1, min=0, max=5, std=2, distrib=DistributionType.NORMAL, save=False)
    assert p.rand(0.0) == 0.0


def test_normal_respects_max():
    assert list(normal(min=8, max=11).rand(BOUNDS)) == pytest.approx([8, 11])


def test_normal_min_only():
    assert list(normal(min=8).rand(BOUNDS)) == pytest.approx([8, math.inf])


def test_normal_max_only():
    assert list(normal(max=11).rand(BOUNDS)) == pytest.approx([-math.inf, 11])


def test_normal_without_bounds():
    p = normal()
    assert list(p.rand(BOUNDS)) == [-math.inf, math.inf]
    assert p.rand(0.5) == pytest.approx(10)


def lognormal():
    return newFloatParam("p", 10, std=0.5, distrib=DistributionType.LOGNORMAL)


def test_lognormal_default_is_median():
    p = lognormal()
    assert p.rand(0.5) == pytest.approx(10)
    # std is the standard deviation of log(x)
    assert np.log(p.rand(norm_cdf(1.0))) - np.log(10) == pytest.approx(0.5)


def norm_cdf(x):
    return 0.5 * (1 + math.erf(x / math.sqrt(2)))


def test_lognormal_matches_brightway():
    """What we sample matches what Brightway (stats_arrays) samples from the exported parameter"""
    p = lognormal()
    data = ProjectParameter.get(ProjectParameter.name == "p").data

    assert data["uncertainty type"] == LognormalUncertainty.id
    params = UncertaintyBase.from_dicts({"loc": data["loc"], "scale": data["scale"]})
    alpha = np.array([0.1, 0.5, 0.9])
    assert list(LognormalUncertainty.ppf(params, alpha.reshape(1, -1))[0]) == pytest.approx(list(p.rand(alpha)))


def test_lognormal_load_roundtrip():
    p = lognormal()
    _param_registry().clear()
    loadParams()
    loaded = _param_registry()["p"]
    assert (loaded.default, loaded.std) == (p.default, p.std)
    assert loaded.rand(0.5) == pytest.approx(10)
