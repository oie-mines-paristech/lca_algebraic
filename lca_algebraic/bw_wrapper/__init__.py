from typing import Callable, Iterable, List, Tuple, Type

import pandas as pd

from .protocol import BrightwayBackend

_INSTALL_MSG = (
    "Brightway is required. Install one of:\n" "  pip install 'lca_algebraic[bw2]'\n" "  pip install 'lca_algebraic[bw25]'"
)


def _bw2data_major_version(bw2data) -> int:
    version = bw2data.__version__
    if isinstance(version, tuple):
        return int(version[0])
    return int(str(version).split(".")[0])


def _load_backend() -> BrightwayBackend:
    try:
        import bw2data
    except ImportError as exc:
        raise ImportError(_INSTALL_MSG) from exc

    if _bw2data_major_version(bw2data) >= 4:
        from .bw25 import BW25Backend

        return BW25Backend()

    from .bw2 import BW2Backend

    return BW2Backend()


_backend = _load_backend()

Activity = _backend.Activity
Exchange = _backend.Exchange
ExchangeDataset = _backend.ExchangeDataset
ActivityProxyBase = _backend.ActivityProxyBase
LCIBackend = _backend.LCIBackend
sqlite3_lci_db = _backend.sqlite3_lci_db
dict_as_exchangedataset = _backend.dict_as_exchangedataset

name = _backend.name

_LAZY_ATTRS = frozenset({"databases", "projects", "methods"})


def __getattr__(attr: str):
    if attr in _LAZY_ATTRS:
        return getattr(_backend, attr)
    raise AttributeError(f"module {__name__!r} has no attribute {attr!r}")


def Database(db_name: str):
    return _backend.Database(db_name)


def Method(key: Tuple):
    return _backend.Method(key)


def get_activity(key):
    return _backend.get_activity(key)


def iter_database(db_name: str) -> Iterable:
    return _backend.iter_database(db_name)


def multi_lca(
    activities: List[dict],
    methods: List[Tuple],
    act_name_fn: Callable,
    method_name_fn: Callable,
) -> pd.DataFrame:
    return _backend.multi_lca(activities, methods, act_name_fn, method_name_fn)


def is_output_exchange(exc) -> bool:
    return _backend.is_output_exchange(exc)


def exchange_type_for_sub_act(sub_act) -> str:
    return _backend.exchange_type_for_sub_act(sub_act)


def process_node_default() -> str:
    return _backend.process_node_default()


def production_edge_default() -> str:
    return _backend.production_edge_default()


def activity_not_found() -> Type[Exception]:
    return _backend.activity_not_found()


def new_database_parameters(params, db_name: str) -> None:
    _backend.new_database_parameters(params, db_name)


def new_project_parameters(params) -> None:
    _backend.new_project_parameters(params)


def copy_activity_skip_keys() -> set:
    return _backend.copy_activity_skip_keys()


def bw2setup() -> None:
    _backend.bw2setup()


def SingleOutputEcospold2Importer(path, dbname, use_mp=False):
    return _backend.SingleOutputEcospold2Importer(path, dbname, use_mp=use_mp)


__all__ = [
    "Activity",
    "ActivityProxyBase",
    "BrightwayBackend",
    "Database",
    "Exchange",
    "ExchangeDataset",
    "LCIBackend",
    "Method",
    "SingleOutputEcospold2Importer",
    "activity_not_found",
    "bw2setup",
    "copy_activity_skip_keys",
    "databases",
    "dict_as_exchangedataset",
    "exchange_type_for_sub_act",
    "get_activity",
    "is_output_exchange",
    "iter_database",
    "methods",
    "multi_lca",
    "name",
    "new_database_parameters",
    "new_project_parameters",
    "process_node_default",
    "production_edge_default",
    "projects",
    "sqlite3_lci_db",
]
