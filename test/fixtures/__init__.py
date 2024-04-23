# Defines basic activities and methods for test
from lca_algebraic import *
from lca_algebraic import getActByCode, newActivity


def init_acts(db):
    """Init bg acts (bio and techno)"""

    return bio1, bio2, bio3


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
