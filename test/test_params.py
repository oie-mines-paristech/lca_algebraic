import math

import numpy as np
import pytest

from lca_algebraic import DistributionType, newFloatParam

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
