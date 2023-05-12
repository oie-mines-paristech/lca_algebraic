from _ast import expr
from typing import Dict

from sympy import parse_expr, Expr, Float

from lca_algebraic import SymDict
from lca_algebraic.lca import _preMultiLCAAlgebric, _lambdify
from lca_algebraic.params import _param_registry
from lca_algebraic.stats import _round_expr
from enum import Enum

def export_lca(
        system,
        functional_units,
        methods_dict,
        axis=None,
        num_digits=3):
    """
    :param system: Root inventory
    :param functional_units : Dict of formula to divide
    :param methods_dict: dict of method_name => method tuple
    :param axis: axis
    :param num_digits: Number of digits
    :return: an instance of "Model"
    """
    lambdas = _preMultiLCAAlgebric(system, list(methods_dict.values()), axis=axis)

    # Simplify
    for lambd in lambdas:
        if isinstance(lambd.expr, SymDict):
            lambd.expr = lambd.expr.dict
        lambd.expr = round_expr(lambd.expr, num_digits=num_digits)

    # Gather all params
    param_names = list(set().union(*[list(lambd.params) for lambd in lambdas]))
    all_params = {param.name: Param.from_ParamDef(param) for param in _param_registry().all() if param.name in param_names}

    # Dict of lambdas
    impacts = {method : Lambda(lambd.expr, all_params) for method, lambd in zip(methods_dict.keys(), lambdas)}
    functional_units = {key: Lambda(expr, all_params) for key, expr in functional_units.items()}

    return Model(
        params=all_params,
        functional_units=functional_units,
        impacts=impacts)

class ParamType(str, Enum) :
    BOOLEAN = "bool"
    ENUM = "enum"
    FLOAT = "float"


def round_expr(exp_or_dict, num_digits):
    if isinstance(exp_or_dict, dict) :
        return dict({key: _round_expr(val, num_digits) for key, val in exp_or_dict.items()})
    else:
        return _round_expr(exp_or_dict, num_digits)

class Lambda :
    """
    This class represents a compiled (lambdified) expression together with the list of requirement parameters and the source expression
    """
    def __init__(self, expr, all_params):

        if isinstance(expr, dict):
            self.lambd = dict()
            all_expanded_params = set()
            for key, expr in expr.items():
                expanded_params = list(str(symbol) for symbol in expr.free_symbols)
                all_expanded_params.update(expanded_params)
                self.lambd[key] = _lambdify(expr, expanded_params)
            self.params = unexpand_param_names(all_params, all_expanded_params)

        else:
            if not isinstance(expr, Expr):
                expr = Float(expr)

            expanded_params = list(str(symbol) for symbol in expr.free_symbols)
            self.lambd = _lambdify(expr, expanded_params)
            self.params = unexpand_param_names(all_params, expanded_params)

        self.expr = expr

    def evaluate(self, all_params, param_values):

        # First, set default values
        values = {key: all_params[key].default for key in self.params}

        # Override with actual values
        values.update({key: val for key, val in param_values.items() if key in self.params})

        # Expand
        expanded_values = dict()
        for param_name, val in values.items():
            param = all_params[param_name]
            expanded_values.update(param.expand_values(val))


        if isinstance(self.lambd, dict) :
            return {key: lambd(**expanded_values) for key, lambd in self.lambd.items()}
        else:
            return self.lambd(**expanded_values)


    def __json__(self):

        if isinstance(self.expr, dict) :
            expr = {key: str(expr) for key, expr in self.expr}
        else:
            expr = str(self.expr)

        return dict(
            params=self.params,
            expr=expr)

    @classmethod
    def from_json(cls, js, all_params):
        expr = js["expr"]
        if isinstance(expr, dict):
            expr = {key:parse_expr(expr) for key, expr in expr.items()}
        else:
            expr = parse_expr(expr)

        return cls(expr=expr, all_params=all_params)

class Param :

    def __init__(self, name, type, unit, default, values=None):
        self.name = name
        self.type = type
        self.default = default
        self.unit = unit
        if values:
            self.values = values

    @classmethod
    def from_json(cls, js):
        return cls(**js)

    @classmethod
    def from_ParamDef(cls, paramDef):
        return cls(
            name=paramDef.name,
            type=paramDef.type,
            unit=paramDef.unit,
            default=paramDef.default,
            values=getattr(paramDef, "values", None))

    def expand_values(self, value):

        # Simple case
        if self.type != ParamType.ENUM :
            return {self.name:value}

        # Enum ? generate individual boolean param values
        return {"%s_%s" % (self.name, enum): 1 if value == enum else 0 for enum in self.values}

    def expand_names(self):

        if self.type != ParamType.ENUM :
            return [self.name]

        # Enum ? generate individual boolean param values
        return ["%s_%s" % (self.name, enum) for enum in self.values]

    def __json__(self):
        return self.__dict__

def expand_param_names(all_params, param_names):
    res = []
    for param_name in param_names :
        param = all_params[param_name]
        res.extend(param.expand_names())
    return res

def unexpand_param_names(all_params, expanded_param_names):
    """Build a dict of expended_param => param"""
    expanded_params_to_params = {name:param.name for param in all_params.values() for name in param.expand_names() }
    return list(set(expanded_params_to_params[name] for name in expanded_param_names))

class Model :

    def __init__(
            self,
            params:Dict,
            impacts:Dict[str, Lambda],
            functional_units:Dict[str, Lambda]) :

        self.params = params
        self.impacts = impacts
        self.functional_units = functional_units

    def __json__(self):
        return self.__dict__

    def evaluate(self, impact, functional_unit, **param_values):

        impact_lambd = self.impacts[impact]
        functional_unit_lambd = self.functional_units[functional_unit]

        fu_val = functional_unit_lambd.evaluate(self.params, param_values)
        impacts = impact_lambd.evaluate(self.params, param_values)
        if isinstance(impacts, dict) :
            return {key: val * fu_val for key, val in impacts.items()}
        else:
            return impacts * fu_val


    @classmethod
    def from_json(cls, js) :
        params = {key: Param.from_json(val) for key, val in js["params"].items()}
        impacts = {key: Lambda.from_json(val, params) for key, val in js["impacts"].items()}
        functional_units = {key: Lambda.from_json(val, params) for key, val in js["functional_units"].items()}
        return cls(params, impacts, functional_units)

def serialize_model(obj) :
    if isinstance(obj, dict) :
        return {key: serialize_model(val) for key, val in obj.items()}

    if hasattr(obj, "__json__") :
        return serialize_model(obj.__json__())

    if hasattr(obj, "__dict__") :
        return serialize_model(obj.__dict__)

    return obj






