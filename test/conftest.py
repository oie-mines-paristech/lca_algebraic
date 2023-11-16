from types import SimpleNamespace

import brightway2
import pytest

from lca_algebraic import initProject, resetDb, resetParams
from lca_algebraic.globs import _clearLCACache
from test.fixtures import init_acts, init_methods
import brightway2 as bw


USER_DB = "fg"
BG_DB= "bg"
METHOD_PREFIX='tests'

@pytest.fixture(autouse=True, scope="module")
def data() :
    """Setup background data """

    # Reset func project, empty DB
    initProject('tests')

    # Create 3 bio activities
    bio1, bio2, bio3 = init_acts(BG_DB)

    # Create one method per bio activity, plus 1 method with several
    ibio1, ibio2, ibio3, imulti = init_methods(BG_DB, METHOD_PREFIX)

    resetDb(USER_DB, True)

    # Pull local vars at global scale
    #for key, val in locals().items():
    #    globals()[key]=val

    return SimpleNamespace(**locals())


@pytest.fixture(autouse=True, scope="function")
def reset_db() :
    """Before each test"""

    for db_name in list(bw.databases) :
        if db_name != BG_DB :
            del bw.databases[db_name]

    resetDb(USER_DB, foreground=True)
    resetParams()
    _clearLCACache()
