import functools
import inspect
from collections import defaultdict
from typing import Union

import brightway2 as bw
import pandas as pd
from bw2data import databases as dbmeta
from bw2data.backends import LCIBackend
from bw2data.proxies import ActivityProxyBase
from typing_extensions import deprecated

from .base_utils import one
from .log import logger

BIOSPHERE_PREFIX = "biosphere"
FOREGROUND_KEY = "fg"


class DbContext:
    """
    Context class specifying the current foreground DB in use. in internal
    Used internally to distinguish database parameters with same names

    usage :
    with DbContext("db") :
        <some code>

    """

    stack = []

    @staticmethod
    def current_db():
        if len(DbContext.stack) == 0:
            return None
        return DbContext.stack[-1]

    def __init__(self, db: Union[str, ActivityProxyBase, LCIBackend]):
        if db is None:
            self.db = None
        elif isinstance(db, ActivityProxyBase):
            self.db = db.key[0]
        elif isinstance(db, str):
            self.db = db
        else:
            self.db = db.name

    def __enter__(self):
        if self.db is not None:
            DbContext.stack.append(self.db)

    def __exit__(self, exc_type, exc_value, exc_traceback):
        if self.db is not None:
            DbContext.stack.pop()


def deleteDb(db_name):
    """Delete a database"""
    del bw.databases[db_name]


def resetDb(db_name, foreground=True):
    """
    Creates and/or cleanup a database.

    Parameters
    ----------

    db_name : str
        Name of the database

    foreground:
        If true (default), the database is set as foreground.
    """
    if db_name in bw.databases:
        logger.warning("Db %s was here. Reseting it" % db_name)
        del bw.databases[db_name]

    db = bw.Database(db_name)
    db.write(dict())
    if foreground:
        setForeground(db_name)
    else:
        setBackground(db_name)


@deprecated("DEPRECATED : Use bw2io.import_ecoinvent_release() instead")
def initProject(project_name):
    """Setup the project if not already done."""
    bw.projects.set_current(project_name)
    bw.bw2setup()


@deprecated("DEPRECATED : Use the new bw2io.import_ecoinvent_release instead")
def importDb(dbname, path, parallel=False):
    """Import eco invent DB

    DEPRECATED : Use the new bw2io.import_ecoinvent_release instead
    """
    if dbname in bw.databases:
        logger.warning("Database '%s' has already been imported " % dbname)
    else:
        ei34 = bw.SingleOutputEcospold2Importer(path, dbname, use_mp=parallel)
        ei34.apply_strategies()
        ei34.statistics()
        ei34.write_database()


_metaCache = defaultdict(lambda: {})


def _setMeta(dbname, key, value):
    """Set meta param on DB"""
    _metaCache[dbname][key] = value

    data = dbmeta[dbname]
    data[key] = value
    dbmeta[dbname] = data
    dbmeta.flush()


def _getMeta(db_name, key):
    if key in _metaCache[db_name]:
        return _metaCache[db_name][key]

    val = dbmeta[db_name].get(key)
    _metaCache[db_name][key] = val
    return val


def _isForeground(db_name):
    """Check is db is marked as foreground DB : which means activities may be parametrized / should be developped."""
    return _getMeta(db_name, FOREGROUND_KEY)


def setForeground(db_name):
    """Set a db as being a foreground database. Foreground databases are considered to be parametric by *lca_algebraic*.
    Internally, their activities are developped as Sympy formulas."""
    return _setMeta(db_name, FOREGROUND_KEY, True)


def setBackground(db_name):
    """
    Set a db as being a background database.

    *lca_algebraic* considers background databases as bring static / non parametric.
    It does not perform the algebraic expansion meccanism on them and just call **Brightway** to compute the impacts for
    their activities.

    """
    return _setMeta(db_name, FOREGROUND_KEY, False)


def _listTechBackgroundDbs():
    """List all background databases technosphere (non biosphere) batabases"""
    return list(name for name in bw.databases if not _isForeground(name) and BIOSPHERE_PREFIX not in name)


def _find_biosphere_db():
    """List all background databases technosphere (non biosphere) batabases"""
    return one(name for name in bw.databases if BIOSPHERE_PREFIX in name)


def list_databases():
    """Returns a pandas dataframe listing all database, their status (foreground/background) and numer of activities"""
    data = list(
        dict(
            name=name,
            backend=_getMeta(name, "backend"),
            nb_activities=len(bw.Database(name)),
            type="biosphere" if BIOSPHERE_PREFIX in name else "foreground" if _isForeground(name) else "background",
        )
        for name in bw.databases
    )

    res = pd.DataFrame(data)
    return res.set_index("name")


def with_db_context(func=None, arg="self"):
    """Internal decorator wrapping function into DbContext, using its first parameters (either Activity, Db or Db name)"""

    if func is None:
        return functools.partial(with_db_context, arg=arg)

    param_specs = inspect.signature(func).parameters

    if arg not in param_specs:
        raise Exception("No param %s in signature of %s" % (arg, func))

    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        # Transform all parameters (positionnal and named) to named ones
        all_param = {k: args[n] if n < len(args) else v.default for n, (k, v) in enumerate(param_specs.items()) if k != "kwargs"}
        all_param.update(kwargs)

        val = all_param[arg]
        if hasattr(val, "key"):
            # value is an activity
            dbname = val.key[0]
        elif isinstance(val, str):
            # Value is directly a  db_name
            dbname = val
        else:
            raise Exception("Param %s is neither an Activity or a db_name : %s" % (arg, val))

        with DbContext(dbname):
            return func(*args, **kwargs)

    return wrapper
