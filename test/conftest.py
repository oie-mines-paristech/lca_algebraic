import warnings
from dataclasses import dataclass
from logging import info

from bw2data.backends import Activity

from test.fixtures import init_methods

import bw2data
import pytest

from lca_algebraic import newActivity, resetDb, resetParams
from lca_algebraic.cache import clear_caches

USER_DB = "fg"
BG_DB = "bg"
METHOD_PREFIX = "tests"
TST_PROJECT_NAME = "tests-bw25"

MethodKey = tuple[str, str, str]


@dataclass
class DataFixture:
    bio1: Activity
    bio2: Activity
    bio3: Activity

    bg_act1: Activity
    bg_act2: Activity
    bg_act3: Activity

    ibio1: MethodKey
    ibio2: MethodKey
    ibio3: MethodKey
    imulti: MethodKey


@pytest.fixture(autouse=True, scope="module")
def data() -> DataFixture:
    """Setup background data"""

    # Reset func project, empty DB
    if TST_PROJECT_NAME in bw2data.projects:
        info("Deleting old tests project")
        try:
            bw2data.projects.delete_project(TST_PROJECT_NAME, delete_dir=True)
        except Exception:
            warnings.warn("Couldn't remove previous test project.", e)

    bw2data.projects.set_current(TST_PROJECT_NAME)
    bw2data.projects

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

    return DataFixture(bio1, bio2, bio3, bg_act1, bg_act2, bg_act3, ibio1, ibio2, ibio3, imulti)


@pytest.fixture(autouse=True, scope="function")
def reset_db():
    """Before each test"""

    for db_name in list(bw2data.databases):
        if db_name != BG_DB and "biosphere" not in db_name:
            del bw2data.databases[db_name]

    resetDb(USER_DB, foreground=True)
    resetParams()
    clear_caches()
