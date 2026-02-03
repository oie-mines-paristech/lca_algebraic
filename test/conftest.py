from dataclasses import dataclass
from logging import info
import pandas as pd
import os

import brightway2 as bw
import pytest

from lca_algebraic import resetDb, resetParams, newActivity, ActivityExtended, getActByCode
from lca_algebraic.cache import clear_caches

USER_DB = "fg"
BG_DB = "bg"
METHOD_PREFIX = "tests"

MethodKey = tuple[str, str, str]

TEST_FOLDER = os.path.dirname(__file__)


@dataclass
class DataFixture:
    bio1: ActivityExtended
    bio2: ActivityExtended
    bio3: ActivityExtended

    bg_act1: ActivityExtended
    bg_act2: ActivityExtended
    bg_act3: ActivityExtended

    ibio1: MethodKey
    ibio2: MethodKey
    ibio3: MethodKey
    imulti: MethodKey


@pytest.fixture(autouse=True, scope="module")
def data() -> DataFixture:
    """Setup background data"""

    # Reset func project, empty DB
    if "tests" in bw.projects:
        info("Deleting old tests project")
        bw.projects.delete_project("tests", delete_dir=True)

    bw.projects.set_current("tests")

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

    for db_name in list(bw.databases):
        if db_name != BG_DB and not db_name.startswith("bio"):
            del bw.databases[db_name]

    resetDb(USER_DB, foreground=True)
    resetParams()
    clear_caches()


def assert_impacts(res: pd.DataFrame, value: float):
    assert res.values[0][0] == value


def init_methods(db, prefix):
    "Create impact methods for bio activities"
    res = []

    # One for each bio act
    for nbio in range(1, 4):
        bioname = "bio" + str(nbio)

        act = getActByCode(db, bioname)

        method = bw.Method((prefix, bioname, "total"))
        method.register(unit="MJ-Eq", description="quantity of " + bioname)
        method.write([(act.key, 1)])

        res.append((prefix, bioname, "total"))

    # Digital : one digit per bio activity
    method = bw.Method((prefix, "all", "total"))
    method.register(unit="1", description="quantity of " + bioname)
    method.write(
        [
            ((db, "bio1"), 1),
            ((db, "bio2"), 2),
            ((db, "bio3"), 4),
        ]
    )
    res.append((prefix, "all", "total"))

    return res
