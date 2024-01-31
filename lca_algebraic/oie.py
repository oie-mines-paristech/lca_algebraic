from pandas import DataFrame
from sympy import Symbol, Basic, Expr
from .params import ParamDef, _completeParamValues, _param_registry


def compute_intermediate_value(formula):
    """Compute intermediate formula with default value for parameters"""
    replace = {param: param.default for param in _param_registry().values()}
    return formula.evalf(subs=replace)








