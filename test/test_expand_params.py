import numpy as np

from lca_algebraic import newBoolParam, newEnumParam, newFloatParam
from lca_algebraic.params import _expand_params, _listOfDictToDictOflist, _param_registry


def test_array_expansion_matches_per_sample():
    newFloatParam("xf", default=1, min=0, max=2, save=False)
    newBoolParam("xb", default=0, save=False)
    newEnumParam("xe", default="a", values=["a", "b", "c"], save=False)
    values = {"xf": [0.5, 1.5, 2.0, 0.1], "xb": [0, 1, 1, 0], "xe": ["a", "c", None, "c"]}

    expected = {}
    for key, vals in values.items():
        per_sample = _listOfDictToDictOflist([_param_registry()[key].expandParams(v) for v in vals])
        expected.update({k: np.array(v, float) for k, v in per_sample.items()})

    res = _expand_params(values)
    assert list(res) == list(expected)
    for key in expected:
        np.testing.assert_array_equal(res[key], expected[key])


def test_none_sample_takes_the_default():
    newFloatParam("nf", default=7, min=0, max=10, save=False)
    np.testing.assert_array_equal(_expand_params({"nf": [1.0, None]})["nf"], [1.0, 7.0])
