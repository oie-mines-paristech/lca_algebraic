import brightway2 as bw
from os import path
import os
import pickle

from .log import logger

LCIA_CACHE = "lcia"
EXPR_CACHE = "expr"

def last_db_update() :
    """Get the last update of current database project"""
    filename = path.join(
        bw.projects.dir, "lci", "databases.db")

    return path.getmtime(filename)

class _CacheDict:
    """ A smart cache that get cleared whenever database changes, and dumped to file whenever we exit from it """
    def __init__(self, name):
        self.name = name
        self.data = None

        filename = _CacheDict.filename(self.name)
        if path.exists(filename) :
            if last_db_update() > path.getmtime(filename) :
                logger.info(f"Db changed recently, clearing cache {self.name}")
                os.remove(filename)
            else:
                # Load cache from disk
                with open(filename, 'rb') as pickleFile:
                    self.data = pickle.load(pickleFile)

        if self.data is None :
            self.data = dict()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):

        # Save data on exit
        if self.data :
            with open(_CacheDict.filename(self.name), 'wb') as pickleFile:
                self.data = pickle.dump(self.data, pickleFile)
    @classmethod
    def filename(cls, name):
        return path.join(bw.projects.dir, f"lca_algebraic_cache-{name}.pickle")


class LCIACache(_CacheDict) :
    def __init__(self):
        _CacheDict.__init__(self, LCIA_CACHE)

class ExprCache(_CacheDict) :
    def __init__(self):
        _CacheDict.__init__(self, EXPR_CACHE)

def clear_caches() :
    for cache_name in [LCIA_CACHE, EXPR_CACHE] :
        filename = _CacheDict.filename(cache_name)
        if path.exists(filename):
            os.remove(filename)