from test.fixtures import init_acts, init_methods
from types import SimpleNamespace

import bw2data
import bw2io
import pytest

from lca_algebraic import newActivity, resetDb, resetParams
from lca_algebraic.cache import clear_caches
from lca_algebraic.settings import Settings

USER_DB = "fg"
BG_DB = "bg"
METHOD_PREFIX = "tests"


@pytest.fixture(autouse=True, scope="module")
def data():
    """Setup background data"""

    # Reset func project, empty DB
    if "tests" not in bw2data.projects:
        bw2data.projects.set_current("Glop")
        bw2data.projects.copy_project("tests", switch=False)
    bw2data.projects.set_current("tests")

    # Clear DB
    resetDb(BG_DB, False)

    # Biosphere activities
    bio1 = newActivity(BG_DB, "bio1", type="emission", unit="kg")
    bio2 = newActivity(BG_DB, "bio2", type="emission", unit="kg")
    bio3 = newActivity(BG_DB, "bio3", type="emission", unit="kg")

    # Process activities
    bg_act1 = newActivity(BG_DB, "bg_act1", "kg", {bio1: 1}, location="GLO")
    bg_act2 = newActivity(BG_DB, "bg_act2", "kg", {bio2: 1}, location="GLO")
    bg_act3 = newActivity(BG_DB, "bg_act3", "kg", {bio3: 1}, location="GLO")

    # Create one method per bio activity, plus 1 method with several
    ibio1, ibio2, ibio3, imulti = init_methods(BG_DB, METHOD_PREFIX)

    # Cleanup foreground DB
    resetDb(USER_DB, True)

    return SimpleNamespace(**locals())


@pytest.fixture(autouse=True, scope="function")
def reset_db():
    """Before each test"""

    for db_name in list(bw2data.databases):
        if db_name != BG_DB and not db_name.startswith("bio"):
            del bw2data.databases[db_name]

    resetDb(USER_DB, foreground=True)
    resetParams()
    clear_caches()
