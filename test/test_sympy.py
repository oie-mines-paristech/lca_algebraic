from sympy import symbols

from lca_algebraic import settings, temp_settings, newFloatParam
from lca_algebraic.settings import CSE
from lca_algebraic.lambda_expression import _lambdify


def test_lambdify_cse_modes():
    p1 = newFloatParam("p1", 0, 0, 1)
    p2 = newFloatParam("p2", 0, 0, 1)

    expr1 = p1 + p2 + 1
    expr = expr1 * expr1

    with temp_settings(cse_mode=CSE.DAG_REFS):
        lambd = _lambdify(expr, [p1, p2])
        res = lambd(p1=1, p2=2)
        assert res == 16

    with temp_settings(cse_mode=CSE.NONE):
        lambd = _lambdify(expr, [p1, p2])
        res = lambd(p1=1, p2=2)
        assert res == 16

    with temp_settings(cse_mode=CSE.DEFAULT):
        lambd = _lambdify(expr, [p1, p2])
        res = lambd(p1=1, p2=2)
        assert res == 16
