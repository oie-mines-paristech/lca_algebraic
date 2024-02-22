from test.fixtures import init_acts, init_methods
from types import SimpleNamespace

import brightway2 as bw
import pytest

from lca_algebraic import initProject, resetDb, resetParams
from lca_algebraic.cache import clear_caches

USER_DB = "fg"
BG_DB = "bg"
METHOD_PREFIX = "tests"


@pytest.fixture(autouse=True, scope="module")
def data():
    """Setup background data"""

    # Reset func project, empty DB
    bw.projects.set_current("tests")
    bw.bw2setup()

    # Create 3 bio activities
    bio1, bio2, bio3 = init_acts(BG_DB)

    # Create one method per bio activity, plus 1 method with several
    ibio1, ibio2, ibio3, imulti = init_methods(BG_DB, METHOD_PREFIX)

    resetDb(USER_DB, True)

    return SimpleNamespace(**locals())


@pytest.fixture(autouse=True, scope="function")
def reset_db():
    """Before each test"""

    for db_name in list(bw.databases):
        if db_name != BG_DB and not db_name.startswith("bio"):
            del bw.databases[db_name]

    resetDb(USER_DB, foreground=True)
    resetParams()
    clear_caches()
