import os
import pickle
from datetime import datetime
from os import path

import brightway2 as bw
import brightway2 as bw2

from .database import _getMeta
from .log import logger
from .settings import PROXY_DB_FLAG, Settings

LCIA_CACHE = "lcia"
EXPR_CACHE = "expr"


# Overide the behaviour for pickling sympy.UndefineFunction
class Pickler(pickle.Pickler):
    from sympy.core.function import UndefinedFunction

    def reducer_override(self, obj):
        # FIXME: maybe too gready, we may check if obj is an instance of
        #        registered functions instead.
        if obj.__class__ is Pickler.UndefinedFunction:
            return type, (obj.__name__, obj.__bases__, dict(obj.__dict__))
        return NotImplemented


def last_db_update():
    """Get the last update of current database project"""
    filename = path.join(bw.projects.dir, "lci", "databases.db")

    return path.getmtime(filename)


def disable_cache():
    Settings.cache_enabled = False


class _Caches:
    "Singleton instance holding caches"
    caches = dict()


def get_dependant_dbs(db_name):
    """Recursively get list of dependant db names, including the current one"""

    # Skip proxy databases
    if _getMeta(db_name, PROXY_DB_FLAG):
        res = set()
    else:
        res = set([db_name])

    for dep in bw2.databases[db_name]["depends"]:
        res.update(get_dependant_dbs(dep))
    return res


def get_last_update(db_name):
    def last_update(db_name):
        return datetime.fromisoformat(bw2.databases[db_name]["modified"])

    res = max(last_update(db) for db in get_dependant_dbs(db_name))
    return res.timestamp()


class _CacheDict:
    """A smart cache that get cleared whenever database changes, and dumped to file whenever we exit from it"""

    def __init__(self, name, db_name):
        self.name = name
        self.db_name = db_name

    def __enter__(self):
        # No cache ? => LOCAL DICT
        if not Settings.cache_enabled:
            self.data = dict()
            return

        filename = self.filename()
        if path.exists(filename):
            if get_last_update(self.db_name) > path.getmtime(filename):
                logger.info(f"Db {self.db_name} changed recently, clearing cache {self.name}")

                # Reset cache on disk and locally
                os.remove(filename)
                _Caches.caches[(self.name, self.db_name)] = dict()

            else:
                # Cache not already loaded in memory ?
                if self.name not in _Caches.caches:
                    # Load cache from disk
                    with open(filename, "rb") as pickleFile:
                        _Caches.caches[(self.name, self.db_name)] = pickle.load(pickleFile)
        else:
            # No file yet, init local cache
            _Caches.caches[(self.name, self.db_name)] = dict()

        # Point to local cache
        self.data = _Caches.caches[(self.name, self.db_name)]

        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        # Save data on exit
        if Settings.cache_enabled and self.data:
            with open(self.filename(), "wb") as pickleFile:
                self.data = Pickler(pickleFile).dump(self.data)

    def filename(self):
        return path.join(bw.projects.dir, f"lca_algebraic_cache-{self.name}-{self.db_name}.pickle")


class LCIACache(_CacheDict):
    def __init__(self, db_name):
        super().__init__(LCIA_CACHE, db_name)


class ExprCache(_CacheDict):
    def __init__(self, db_name):
        super().__init__(EXPR_CACHE, db_name)


def clear_caches(local=True, disk=True):
    if local:
        _Caches.caches = dict()

    if disk:
        for db_name in bw2.databases:
            for cache_name in [LCIA_CACHE, EXPR_CACHE]:
                filename = _CacheDict(cache_name, db_name).filename()
                if path.exists(filename):
                    os.remove(filename)
