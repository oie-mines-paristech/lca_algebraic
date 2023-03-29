from pandas import DataFrame
from sympy import Symbol, Basic, Expr
from .params import ParamDef, _completeParamValues, _param_registry


def compute_intermediate_value(formula):
    """Compute intermediate formula with default value for parameters"""
    replace = {param: param.default for param in _param_registry().values()}
    return formula.evalf(subs=replace)



def intermediate_values_table(globals) :
    """Display a table of all intermediate formula : their name, formula and default value."""

    res = []
    for name, formula in list(globals.items()) :
        if isinstance(formula, Expr) and not isinstance(formula, ParamDef) :
            has_params = any(symbol for symbol in formula.free_symbols if isinstance(symbol, ParamDef))
            if not has_params :
                continue
            res.append(dict(
                name=name,
                formula=formula,
                value=compute_intermediate_value(formula)))

    res = DataFrame.from_dict(res)

    if len(res) > 0 :
        res = res.set_index("name")

    return res




