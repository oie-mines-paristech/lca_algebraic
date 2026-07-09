from sympy import Expr


class AxisDict(dict):
    """This class acts like a dict with arithmetic operations. It is useful to process 'axes' LCA computations"""

    def __init__(self, src: dict = None, **kwargs):
        # Create an empty dict
        super().__init__()

        # Add src
        if src is not None:
            super().update(src)

        # add kwargs if any
        super().update(kwargs)

    def _apply_op(self, other, fop, null_val):
        # None is the key for non flagged values
        if isinstance(other, int | float | Expr):
            return AxisDict({k: fop(v, other) for k, v in self.items()})
        elif isinstance(other, AxisDict):
            all_keys = set(other.keys()) | set(self.keys())
            return AxisDict({key: fop(self.get(key, null_val), other.get(key, null_val)) for key in all_keys})
        raise Exception("AxisDict: unsuported operation")

    def __repr__(self):
        """Custom representation that returns string as key instead of symbols"""
        return "{%s}" % ",".join("'%s': %s" % (k.__repr__(), v.__repr__()) for k, v in self.items())

    def __str__(self):
        return self.__repr__()

    def __add__(self, other):
        return self._apply_op(other, lambda a, b: a + b, 0)

    def __sub__(self, other):
        return self._apply_op(other, lambda a, b: a - b, 0)

    def __rsub__(self, other):
        return self._apply_op(other, lambda a, b: b - a, 0)

    def __radd__(self, other):
        return self._apply_op(other, lambda a, b: b + a, 0)

    def __mul__(self, other):
        return self._apply_op(other, lambda a, b: a * b, 1)

    def __rmul__(self, other):
        return self._apply_op(other, lambda a, b: a * b, 1)

    def __truediv__(self, other):
        return self._apply_op(other, lambda a, b: a / b, 1)

    def __rtruediv__(self, other):
        return self._apply_op(other, lambda a, b: b / a, 1)

    def _deprecated_equals(self, other):
        return isinstance(other, AxisDict) and self._dict == other._dict

    @property
    def free_symbols(self):
        """Only return free symbol for values (not keys)"""
        res = set()
        for key, val in self.items():
            if hasattr(val, "free_symbols"):
                res |= val.free_symbols
        return res
