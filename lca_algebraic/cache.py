import atexit
import os
import signal
import sys
import time
from base64 import b32encode
from collections.abc import MutableMapping
from datetime import datetime
from os import path
from typing import Dict, Tuple

from dill import Pickler, load
from sympy.core.function import UndefinedFunction

from lca_algebraic.bw_wrapper import Database, databases, projects

from .database import _getMeta
from .log import logger
from .settings import PROXY_DB_FLAG, Settings

LCIA_CACHE = "lcia"
EXPR_CACHE = "expr"
MAPPING_CACHE = "mapping"


def _filename_encode(s: str):
    return b32encode(s.encode("utf-8")).decode("utf-8")


# Overide the behaviour for pickling sympy.UndefineFunction
class MyPickler(Pickler):
    def reducer_override(self, obj):
        # FIXME: maybe too gready, we may check if obj is an instance of
        #        registered functions instead.
        if obj.__class__ is UndefinedFunction:
            return type, (obj.__name__, obj.__bases__, dict(obj.__dict__))
        return NotImplemented


def last_db_update():
    """Get the last update of current database project"""
    filename = path.join(projects.dir, "lci", "databases.db")
    return path.getmtime(filename)


def disable_cache():
    Settings.cache_enabled = False


def flush_caches():
    for cache in _Caches.caches.values():
        cache.sync()


def flush_caches_on_exit():
    """Best-effort flush of the caches, safe to call several times.

    Does nothing when ``Settings.auto_flush_on_exit`` is disabled. Errors are
    swallowed (there is nothing useful to do while the interpreter is exiting).
    """
    if not Settings.auto_flush_on_exit:
        return
    try:
        flush_caches()
    except Exception as e:
        logger.warning(f"Could not flush caches on exit : {e}")


_exit_hooks_installed = False


def install_exit_hooks():
    """Flush the caches when the interpreter stops : normal exit, uncaught
    exception, or SIGINT / SIGTERM (see ``Settings.auto_flush_on_exit``).

    Idempotent. Useful when ``Settings.auto_flush`` is off, so that a run (or a
    crash, or a Ctrl-C) does not lose the LCIA/expression caches computed so far.
    SIGKILL cannot be caught.
    """
    global _exit_hooks_installed
    if _exit_hooks_installed:
        return
    _exit_hooks_installed = True

    # Normal exit
    atexit.register(flush_caches_on_exit)

    # Uncaught exception : flush, then keep the default reporting
    _previous_excepthook = sys.excepthook

    def _excepthook(exc_type, exc, tb):
        flush_caches_on_exit()
        _previous_excepthook(exc_type, exc, tb)

    sys.excepthook = _excepthook

    # Termination signals : flush, then die with the default disposition
    def _signal_handler(signum, _frame):
        flush_caches_on_exit()
        signal.signal(signum, signal.SIG_DFL)
        os.kill(os.getpid(), signum)

    for _sig in (signal.SIGINT, signal.SIGTERM):
        try:
            signal.signal(_sig, _signal_handler)
        except (ValueError, OSError):
            # Not in the main thread (some notebook kernels) : signals are skipped
            pass


def get_dependant_dbs(db_name, skip_proxy=None):
    """Recursively get list of dependant db names, including the current one

    Proxy databases are skipped unless background factorization is enabled : with
    ``factorize_static_bg`` the developed expressions / LCIA of a database refer to
    the acts of its ``<db>-proxy`` databases, so a change of those proxies must be
    visible (see the staleness check in ``_CacheDict``).
    """
    if skip_proxy is None:
        skip_proxy = not Settings.factorize_static_bg

    # Only work if database is up to date (processed)
    if databases[db_name].get("dirty"):
        Database(db_name).process()

    # Skip proxy databases
    if skip_proxy and _getMeta(db_name, PROXY_DB_FLAG):
        res = set()
    else:
        res = set([db_name])

    for dep in databases[db_name]["depends"]:
        res.update(get_dependant_dbs(dep))
    return res


def get_last_update(db_name):
    def last_update(db_name):
        return datetime.fromisoformat(databases[db_name]["modified"])

    dependant_dbs = get_dependant_dbs(db_name)

    if len(dependant_dbs) == 0:
        logger.warning(f"Fonud no dependant dbs found for {db_name}. Returning empty update time")
        return 0

    res = max(last_update(db) for db in dependant_dbs)
    return res.timestamp()


class SyncDict(MutableMapping):
    """
    A dict tat loads its values from a file, track the latest updates, and sync its content to a file
    """

    def __init__(self, name, db_name):
        self._data = {}
        self.name = name
        self.db_name = db_name
        self.last_update = 0.0
        self.load()

    def _filename(self):
        return path.join(
            projects.dir, f"lca_algebraic_cache-{_filename_encode(self.name)}-{_filename_encode(self.db_name)}.pickle"
        )

    # ---------- core MutableMapping ----------
    def __getitem__(self, key):
        return self._data[key]

    def __setitem__(self, key, value):
        self._data[key] = value
        self._touch()

    def __delitem__(self, key):
        del self._data[key]
        self._touch()

    def __iter__(self):
        return iter(self._data)

    def __len__(self):
        return len(self._data)

    # ---------- mutation helpers ----------
    def _touch(self):
        self.last_update = time.time()

    def clear(self, disk=False):
        self._data.clear()
        if disk:
            if os.path.exists(self._filename()):
                os.remove(self._filename())
        self._touch()

    def update(self, *args, **kwargs):
        self._data.update(*args, **kwargs)
        self._touch()

    # ---------- persistence ----------
    def load(self):
        if not Settings.cache_enabled:
            return

        if not os.path.exists(self._filename()):
            return

        with open(self._filename(), "rb") as f:
            try:
                self._data = load(f)
            except Exception as e:
                logger.error(f"Error while loading cache {self._filename()}: {e}. Ignoring and overriding it")

        self.last_update = os.path.getmtime(self._filename())

    def sync(self):
        if not Settings.cache_enabled:
            return

        if len(self._data) == 0:
            return

        if os.path.exists(self._filename()):
            file_mtime = os.path.getmtime(self._filename())
        else:
            file_mtime = 0.0

        if self.last_update <= file_mtime:
            return

        # info(f"Flushing cache {self.name} / {self.db_name}")

        tmp = self._filename() + ".tmp"
        with open(tmp, "wb") as f:
            pickler = MyPickler(f)
            pickler.dump(self._data)

        os.replace(tmp, self._filename())
        return True


class _Caches:
    """Singleton instance holding caches"""

    caches: Dict[Tuple[str, str], SyncDict] = dict()


class _CacheDict:
    """A smart cache that get cleared whenever database changes, and dumped to file whenever we exit from it"""

    def __init__(self, name, db_name):
        self.name = name
        self.db_name = db_name

        key = (name, db_name)

        # Not initialized yet ?

        if key not in _Caches.caches:
            _Caches.caches[key] = SyncDict(name, db_name)

        # Data points to cache
        self.data = _Caches.caches[key]

        # Db more recent ? clean it
        if os.path.exists(self.data._filename()) and (get_last_update(db_name) > self.data.last_update):
            logger.info(f"Db changed recently, clearing cache {self.name}")
            self.data.clear(disk=True)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        # Save data when exiting the context manager (or only once at exit)
        if Settings.auto_flush:
            self.data.sync()


class LCIACache(_CacheDict):
    def __init__(self, db_name):
        super().__init__(LCIA_CACHE, db_name)


class ExprCache(_CacheDict):
    def __init__(self, db_name):
        super().__init__(EXPR_CACHE, db_name)


class MappingCache(_CacheDict):
    def __init__(self, db_name):
        super().__init__(MAPPING_CACHE, db_name)


def clear_caches(local=True, disk=True):
    for cache_name in [LCIA_CACHE, EXPR_CACHE]:
        for db_name in databases:
            cache = SyncDict(cache_name, db_name)
            if disk:
                cache.clear(disk=True)
            key = (cache_name, db_name)
            if local and key in _Caches.caches:
                del _Caches.caches[key]
