import pytest
import sympy

from sympy import Symbol, Min, Max

from lca_algebraic import newFloatParam
from numpy import array


def test_simplify_sums():
    a = newFloatParam("a", min=1e-3, max=2e-3, default=1.5e-3)
    b = newFloatParam("b", min=-2e-3, max=-1e-3, default=-1.5e-3)
    c = newFloatParam("c", min=3.0, max=4.0, default=3.5)
    d = newFloatParam("d", min=-4.0, max=-3.0, default=-3.5)

    e0 = a + b + c + d

    param_values = {
        "a": [a.min, a.max],
        "b": [b.min, b.max],
        "c": [c.min, c.max],
        "d": [d.min, d.max],
    }

    import lca_algebraic.stats

    e0sum = lca_algebraic.stats._simplify_sums(e0, param_values)
    e0prd = lca_algebraic.stats._simplify_products(e0, param_values)

    assert e0sum == (d + c)
    assert e0prd == e0


def test_simplify_products():
    a = newFloatParam("a", min=0.999, max=1.001, default=1.0)
    b = newFloatParam("b", min=-1.001, max=-0.999, default=-1.0)
    c = newFloatParam("c", min=3.0, max=4.0, default=3.5)
    d = newFloatParam("d", min=-4.0, max=-3.0, default=-3.5)

    e0 = a * c + b * d

    print(e0)

    param_values = {
        "a": array([a.min, a.max]),
        "b": array([b.min, b.max]),
        "c": array([c.min, c.max]),
        "d": array([d.min, d.max]),
    }

    import lca_algebraic.stats

    e0sum = lca_algebraic.stats._simplify_sums(e0, param_values)
    e0prd = lca_algebraic.stats._simplify_products(e0, param_values)

    assert e0sum == e0
    assert e0prd == (c - d)


def test_simplify_min_max():
    a = newFloatParam("a", min=1.0, max=2.0, default=1.5)
    b = newFloatParam("b", min=-2.0, max=-1.0, default=-1.5)
    c = newFloatParam("c", min=3.0, max=4.0, default=3.5)
    d = newFloatParam("d", min=-4.0, max=-3.0, default=-3.5)

    e0 = Min(a, b, c, d)
    e1 = Max(a, b, c, d)

    param_values = {
        "a": array([a.min, a.max]),
        "b": array([b.min, b.max]),
        "c": array([c.min, c.max]),
        "d": array([d.min, d.max]),
    }

    import lca_algebraic.stats

    e0min = lca_algebraic.stats._simplify_min(e0, param_values)
    e0max = lca_algebraic.stats._simplify_max(e0, param_values)
    e1min = lca_algebraic.stats._simplify_min(e1, param_values)
    e1max = lca_algebraic.stats._simplify_max(e1, param_values)

    assert e0min == d
    assert e0max == e0
    assert e1min == e1
    assert e1max == c


def test_sobols_failure_is_reported(monkeypatch):
    """A failed Sobol analysis is reported with the method name, not hidden by an error in the handler"""
    import lca_algebraic.stats as stats

    warnings = []
    monkeypatch.setattr(stats, "_parallel_map", lambda f, items: [(0, {})])
    monkeypatch.setattr(stats, "warn", lambda *args: warnings.append(args))

    stats._sobols([("EF v3.1", "climate change", "GWP100")], {"names": ["a"]}, None)

    assert warnings[0][0] == "Sobol failed on climate change - GWP100"


def test_sobol_without_second_order(data):
    """second_order=False skips the second-order samples, leaves s2 empty and still simplifies the model"""
    import lca_algebraic.stats as stats
    from lca_algebraic import newActivity, sobol_simplify_model
    from test.conftest import USER_DB

    p1 = newFloatParam("p1", 1, min=1, max=2)
    p2 = newFloatParam("p2", 1, min=0.001, max=0.001)
    m1 = newActivity(USER_DB, "m1", "kg", {data.bio1: p1 * (p1 + 0.001 * p1 + p2)})

    problem, params, Y = stats._stochastics(m1, [data.ibio1], 64, var_params=[p1, p2], second_order=False)
    assert len(params["p1"]) == 64 * (2 + 2)

    sob = stats._sobols([data.ibio1], problem, Y, second_order=False)
    assert sob.s2 is None and sob.s2_conf is None
    assert sob.s1.shape == (2, 1)

    res = sobol_simplify_model(m1, [data.ibio1], n=64, simple_products=False, second_order=False)[0]
    assert res.expr.__repr__() == "1.0*p1**2"
