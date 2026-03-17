import pytest
import os.path

from lca_algebraic import newActivity, newFloatParam, compute_impacts
from lca_algebraic.edges_utils import register_custom_edge_method
from test.conftest import TEST_FOLDER, DataFixture, USER_DB, assert_impacts

EDGE_METHOD_KEY = ("edge", "total", "bio1")
EDGE_METHOD_FILE = os.path.join(TEST_FOLDER, "res", "custom_edge_method.json")


@pytest.fixture(autouse=True, scope="module")
def setup_edge_method():
    register_custom_edge_method(EDGE_METHOD_KEY, EDGE_METHOD_FILE)


"""
def test_parametrized_lcia(data: DataFixture):
    p1 = newFloatParam("p1", default=1, min=0, max=10)

    fg = newActivity(db_name=USER_DB, name="test_act", unit="kg", exchanges={data.bg_act1: p1})

    res = compute_impacts(models=fg, methods=EDGE_METHOD_KEY, p1=1)

    assert_impacts(res, 1.0)

    # Shoulduse cache of edge impacts
    res = compute_impacts(models=fg, methods=EDGE_METHOD_KEY, p1=2.0)

    assert_impacts(res, 2.0)
"""
