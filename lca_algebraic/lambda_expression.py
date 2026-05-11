from copy import copy
from dataclasses import dataclass
from typing import Dict, List

from bw2data.backends.peewee import Activity
from sympy import Basic, Expr, lambdify, parse_expr
from sympy.printing.numpy import NumPyPrinter

from lca_algebraic.axis_dict import AxisDict
from lca_algebraic.base_utils import _user_functions
from lca_algebraic.params import (
    _complete_params,
    _expand_param_names,
    _expand_params,
    _expanded_names_to_names,
)
from lca_algebraic.settings import Settings


@dataclass
class ValueContext:
    """Represents a result value, with all parameters values used in context"""

    value: float
    context: Dict[str, float]


class LambdaExpr:
    """
    This class represents a compiled (lambdified) expression together
    with the list of requirement parameters and the source expression
    """

    def __init__(self, expr: Expr, expanded_params=None, params=None, sobols=None):
        """Computes a lamdda function from expression and list of expected parameters.
        you can provide either the list pf expanded parameters (full vars for enums) for the 'user' param names
        """

        # List of background activities, indexed in the same order as impacts[i]
        self.background_activities: List[Activity] = list()
        self.impacts: List[float] = None

        if isinstance(expr, dict):
            # Come from JSON serialization
            obj = expr
            # LIst of required params for this lambda
            self.params: List[str] = obj["params"]

            # Full names
            self.expanded_params = _expand_param_names(self.params)
            local_dict = {x[0].name: x[0] for x in _user_functions.values()}
            self.expr = parse_expr(obj["expr"], local_dict=local_dict)
            self.lambd = _lambdify(self.expr, self.expanded_params + ["impacts"])
            self.sobols = obj["sobols"]

        else:
            self.expr = expr
            self.params = params

            if expanded_params is None:
                if params is None:
                    expanded_params = _free_symbols(expr)
                    params = _expanded_names_to_names(expanded_params)
                    self.params = params

                # We expand again the parameters
                # If we expect an enum param name, we also expect the other ones :
                # enumparam_val1 => enumparam_val1, enumparam_val2, ...
                expanded_params = _expand_param_names(params)

            elif self.params is None:
                self.params = _expanded_names_to_names(expanded_params)

            self.lambd = _lambdify(expr, expanded_params + ["impacts"])
            self.expanded_params = expanded_params
            self.sobols = sobols

    @property
    def has_axis(self):
        return isinstance(self.expr, AxisDict)

    @property
    def axis_keys(self):
        if self.has_axis:
            return self.expr.str_keys()
        else:
            return None

    def compute(self, **params) -> ValueContext:
        """Compute result value based of input parameters"""

        if self.impacts is None:
            raise Exception("Trying to compute impacts on a generic LambdaExpr object, without impact values attached.")

        # Add default or computed values
        completed_params = _complete_params(params, self.params)

        # Expand enums
        expanded_params = _expand_params(completed_params)

        # Remove parameters that are not required
        expanded_params = _filter_param_values(expanded_params, self.expanded_params)

        value = self.lambd(impacts=self.impacts, **expanded_params)

        return ValueContext(value=value, context=completed_params)

    def serialize(self):
        expr = str(self.expr)
        return dict(params=self.params, expr=expr, sobols=self.sobols)

    @staticmethod
    def use_sympy_cse(b=True):
        LambdaExpr._use_sympy_cse = b

    def __repr__(self):
        return repr(self.expr)

    def _repr_latex_(self):
        return self.expr._repr_latex_()

    def with_impacts(self, impacts: Dict[Activity, float]):
        """Transforms a generic LambaExpr toa one with custom impact values"""
        if self.impacts is not None:
            raise Exception("You are specializing a LambaExpression with impacts already set.")
        res = copy(self)
        res.impacts = [impacts[bg_act] for bg_act in self.background_activities]

        return res


def _filter_param_values(params, expanded_param_names):
    return {key: val for key, val in params.items() if key in expanded_param_names}


def _free_symbols(expr: Basic):
    if isinstance(expr, Basic):
        return set([str(symb) for symb in expr.free_symbols if "impacts" not in str(symb)])
    else:
        # Static value
        return set()


class LambdWrapper:
    """Wrapper of lambda function. required for pickling in cache"""

    def __init__(self, lambd):
        self.lambd = lambd

    def __call__(self, *args, **kwargs):
        res = self.lambd(*args, **kwargs)
        if isinstance(res, dict):
            # Transform key symbols into Str
            return {str(k): v for k, v in res.items()}
        else:
            return res


def _lambdify(expr: Basic, expanded_params):
    """Lambdify, handling manually the case of SymDict (for impacts by axis)"""

    printer = NumPyPrinter(
        {"fully_qualified_modules": False, "inline": True, "allow_unknown_functions": True, "user_functions": dict()}
    )

    modules = [{x[0].name: x[1] for x in _user_functions.values()}, "numpy"]

    if isinstance(expr, Basic):
        lambd = lambdify(expanded_params, expr, modules, printer=printer, cse=Settings.lambdify_cse)
        return LambdWrapper(lambd=lambd)

    else:
        # Not an expression : return static func
        def static_func(*args, **kargs):
            return expr

        return static_func
