from sympy import symbols, lambdify, simplify

from lca_algebraic.axis_dict import AxisDict, NO_AXIS


def test_sum():
    a, b = symbols("a b")

    a1 = AxisDict({a: 1})
    a2 = AxisDict({a: 2})
    b2 = AxisDict({b: 2})

    assert a1 + b2 == AxisDict({a: 1, b: 2})
    assert a1 + a2 == AxisDict({a: 3})
    assert a1 + 1 == AxisDict({a: 1, NO_AXIS: 1})


def test_mul():
    a = symbols("a")

    a1 = AxisDict({a: 2})
    assert a1 * 2 == AxisDict({a: 4})
    assert simplify(a1 / 2) == AxisDict({a: 1})


def test_free_symbols():
    dic = AxisDict({"a": "b"})
    assert dic.free_symbols == set([symbols("b")])


def test_lambdify():
    a, b = symbols("a b")

    a1 = AxisDict({a: b * 2})

    lambd = lambdify([b], a1)

    res = lambd(2)

    assert res == {a: 4}
