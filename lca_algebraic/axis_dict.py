from sympy import symbols
from sympy.core.containers import Dict as SympyDict

NO_AXIS_NAME = "_other_"
NO_AXIS = symbols(NO_AXIS_NAME)


class AxisDict(SympyDict):
    """This class acts like a dict with arithmetic operations. It is useful to process 'axes' LCA computations"""

    def _apply_op(self, other, fop, null_val):
        # None is the key for non flagged values
        if not isinstance(other, AxisDict):
            other = AxisDict({NO_AXIS: other})

        all_keys = set(other._dict.keys()) | set(self._dict.keys())
        return AxisDict({key: fop(self._dict.get(key, null_val), other._dict.get(key, null_val)) for key in all_keys})

    def __repr__(self):
        """Custom representation that returns string as key instead of symbols"""
        return "{%s}" % ",".join("'%s': %s" % (k.__repr__(), v.__repr__()) for k, v in self._dict.items())

    def __str__(self):
        return self.__repr__()

    def _apply_self(self, fop):
        return AxisDict({key: fop(val) for key, val in self._dict.items()})

    def __add__(self, other):
        return self._apply_op(other, lambda a, b: a + b, 0)

    def __sub__(self, other):
        return self._apply_op(other, lambda a, b: a - b, 0)

    def __rsub__(self, other):
        return self._apply_op(other, lambda a, b: b - a, 0)

    def __radd__(self, other):
        return self._apply_op(other, lambda a, b: b + a, 0)

    def __mul__(self, other):
        return self._apply_self(lambda a: a * other)

    def __rmul__(self, other):
        return self._apply_self(lambda a: other * a)

    def __truediv__(self, other):
        return self._apply_self(lambda a: a / other)

    def __rtruediv__(self, other):
        return self._apply_self(lambda a: other / a)

    def as_coeff_Mul(self, rational=False):
        """Efficiently extract the coefficient of a product."""
        return 1, self

    def _defer(self, funcname, args, kwargs):
        return AxisDict(
            {
                key: val if not hasattr(val, funcname) else getattr(val, funcname)(*args, **kwargs)
                for key, val in self._dict.items()
            }
        )

    def str_keys(self):
        # Return a list to ensure the order is kept
        return list(str(key) for key in self._dict.keys())

    def equals(self, other):
        return isinstance(other, AxisDict) and self._dict == other._dict

    @property
    def free_symbols(self):
        """Only return free symbol for values (not keys)"""
        res = set()
        for key, val in self._dict.items():
            if hasattr(val, "free_symbols"):
                res |= val.free_symbols
        return res
