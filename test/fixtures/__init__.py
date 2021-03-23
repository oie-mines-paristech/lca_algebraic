# Defines basic activities and methods for test
from lca_algebraic import *

def init_acts(db) :
    """Init bg acts (bio and techno) """

    # Clear DB
    resetDb(db)

    # Biosphere activities
    bio1 = newActivity(db, "bio1", "unit", type="biosphere")
    bio2 = newActivity(db, "bio2", "unit", type="biosphere")
    bio3 = newActivity(db, "bio3", "unit", type="biosphere")

    # Process activities
    bg_act1 = newActivity(db, "bg_act1", "kg", {
        bio1 : 1,
        bio2 : 2
    })
    bg_act2 = newActivity(db, "bg_act2", "kg", {
        bio1 : 2,
        bio2 : 1
    })
    bg_act3 = newActivity(db, "bg_act3", "m3" , {
        bio3: 1,
    })


def init_methods(db, prefix) :
    "Create impact methods for bio activities"
    res = []

    # One for each bio act
    for nbio in range(1, 3) :
        bioname = "bio" + str(nbio)

        act = getActByCode(db, bioname)

        method = bw.Method((prefix, bioname, 'total'))
        method.register(
            unit='MJ-Eq',
            description='quantity of ' + bioname)
        method.write([(act.key, 1)])

        res.append(method)

    # Digital : one digit per bio activity
    method = bw.Method((prefix, "all", "total"))
    method.register(
        unit='1',
        description='quantity of ' + bioname)
    method.write([
        ((db, "bio1"), 1),
        ((db, "bio2"), 2),
        ((db, "bio3"), 4),
    ])
    res.append(method)

    return res
