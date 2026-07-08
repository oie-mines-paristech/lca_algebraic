from lca_algebraic.axis_dict import AxisDict
from sympy import symbols

def test_axisdict_op_axisdict():

    def _axisdict_op_axisdict(op, identity_value):
        a1 = AxisDict({"a": 11, "b": 21, "c": 31})
        b1 = AxisDict({"a": 12, "b": 22})
        c1 = op(a1, b1)
        d1 = op(b1, a1)

        # check that a1 and b1 aren't modified
        assert a1 == AxisDict({"a": 11, "b": 21, "c": 31})
        assert b1 == AxisDict({"a": 12, "b": 22})
        assert c1 == AxisDict({"a": op(11,12), "b": op(21,22), "c": op(31,identity_value)})
        assert d1 == AxisDict({"a": op(12,11), "b": op(22,21), "c": op(identity_value,31)})

    _axisdict_op_axisdict(lambda a, b: a * b, 1)
    _axisdict_op_axisdict(lambda a, b: a + b, 0)
    _axisdict_op_axisdict(lambda a, b: a - b, 0)
    _axisdict_op_axisdict(lambda a, b: a / b, 1)

def test_free_symbols():
    dic = AxisDict({"a": symbols("b")})
    assert dic.free_symbols == set([symbols("b")])

def test_scalar_op_axisdict():

    def _scalar_op_axisdict(op):
        a1 = AxisDict({"a": 11, "b": 21})
        b1 = op(3, a1)
        c1 = op(a1, 3)

        # check that a1 is not modified
        assert a1 == AxisDict({"a": 11, "b": 21})
        assert b1 == AxisDict({"a": op(3,11), "b": op(3,21)})
        assert c1 == AxisDict({"a": op(11,3), "b": op(21,3)})

    _scalar_op_axisdict(lambda a, b: a * b)
    _scalar_op_axisdict(lambda a, b: a + b)
    _scalar_op_axisdict(lambda a, b: a - b)
    _scalar_op_axisdict(lambda a, b: a / b)

