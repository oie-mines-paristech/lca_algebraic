import os
import pickle
from os import path

import brightway2 as bw

from .log import logger
from .settings import Settings

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


class _CacheDict:
    """A smart cache that get cleared whenever database changes, and dumped to file whenever we exit from it"""

    def __init__(self, name):
        self.name = name

        # No cache ? => LOCAL DICT
        if not Settings.cache_enabled:
            self.data = dict()
            return

        filename = _CacheDict.filename(self.name)
        if path.exists(filename):
            if last_db_update() > path.getmtime(filename):
                logger.info(f"Db changed recently, clearing cache {self.name}")

                # Reset cache on disk and locally
                os.remove(filename)
                _Caches.caches[name] = dict()

            else:
                # Cache not already loaded in memory ?
                if name not in _Caches.caches:
                    # Load cache from disk
                    with open(filename, "rb") as pickleFile:
                        _Caches.caches[name] = pickle.load(pickleFile)
        else:
            # No file yet, init local cache
            _Caches.caches[name] = dict()

        # Point to local cache
        self.data = _Caches.caches[name]

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        # Save data on exit
        if Settings.cache_enabled and self.data:
            with open(_CacheDict.filename(self.name), "wb") as pickleFile:
                self.data = Pickler(pickleFile).dump(self.data)

    @classmethod
    def filename(cls, name):
        return path.join(bw.projects.dir, f"lca_algebraic_cache-{name}.pickle")


class LCIACache(_CacheDict):
    def __init__(self):
        _CacheDict.__init__(self, LCIA_CACHE)


class ExprCache(_CacheDict):
    def __init__(self):
        _CacheDict.__init__(self, EXPR_CACHE)


def clear_caches(local=True, disk=True):
    if local:
        _Caches.caches = dict()

    if disk:
        for cache_name in [LCIA_CACHE, EXPR_CACHE]:
            filename = _CacheDict.filename(cache_name)
            if path.exists(filename):
                os.remove(filename)
