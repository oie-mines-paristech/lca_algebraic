import json
from functools import cache
from os import path
from typing import Dict

from bw2calc import LCA
from edges import EdgeLCIA, get_available_methods
from edges.filesystem_constants import DATA_DIR
from edges.georesolver import GeoResolver
from pypardiso import factorized

from lca_algebraic import ActivityExtended
from lca_algebraic.cache import _CacheDict
from lca_algebraic.log import info, logger

_CUSTOM_META: Dict[tuple, dict] = {}


def register_custom_edge_method(key, filename):
    _CUSTOM_META[key] = _load_metadata(filename)


def get_edge_methods_metadata() -> Dict[tuple, dict]:
    """REturn dict of edge metatada. both builtin and custom registered ones"""
    res = _load_builtin_metadata().copy()
    res.update(_CUSTOM_META)
    return res


@cache
def _load_metadata(filename):
    """Load metadata for single method"""
    with open(filename, "rt") as f:
        res = json.load(f)
        res["filename"] = filename
        del res["exchanges"]
        return res


@cache
def _load_builtin_metadata():
    """Load metadata for builtin edge methods"""
    res = dict()
    for method_key in get_available_methods():
        key = "_".join(method_key)
        filename = path.join(DATA_DIR, f"{key}.json")
        res[method_key] = _load_metadata(filename)
    return res


def setup_edges_serialization():
    """Function for settings custom methods to enable pickling / unpickling of EdgeLCIA"""

    def _getstate_remove_logger(obj):
        state = obj.__dict__.copy()
        if "logger" in state:
            del state["logger"]
        return state

    def _setstate_reset_logger(obj, state: dict):
        obj.__dict__.update(state)
        obj.logger = logger

    EdgeLCIA.__getstate__ = _getstate_remove_logger
    EdgeLCIA.__setstate__ = _setstate_reset_logger
    GeoResolver.__getstate__ = _getstate_remove_logger
    GeoResolver.__setstate__ = _setstate_reset_logger

    def _getstate_lca(lca: LCA):
        state = _getstate_remove_logger(lca)
        del state["solver"]
        return state

    def _setstate_lca(lca: LCA, state: dict):
        _setstate_reset_logger(lca, state)
        lca.solver = factorized(lca.technosphere_matrix.tocsc())

    LCA.__getstate__ = _getstate_lca
    LCA.__setstate__ = _setstate_lca


class EdgeCache(_CacheDict):
    """Custom cache for EdgeLCIA. We use one separate cache per method, because each pickled file is big (100Mb)"""

    def __init__(self, db_name, method: tuple):
        key = "_".join(method)
        super().__init__(f"edge_lcia_{key}", db_name)


def setup_edge_lcia(method_key, act: ActivityExtended):
    """Done once then cached"""

    # Pass either a tuple or file path
    method = method_key

    if method_key in _CUSTOM_META:
        method = _CUSTOM_META[method_key]["filename"]

    lcia = EdgeLCIA(demand={act: 1}, method=method)

    lcia.lci()
    lcia.map_exchanges()
    lcia.map_aggregate_locations()
    lcia.map_dynamic_locations()
    lcia.map_contained_locations()
    lcia.map_remaining_locations_to_global()
    lcia.evaluate_cfs()

    return lcia


MAIN_KEY = "lcia"


def compute_edge_impacts(db_name: str, method: tuple, acts: list[ActivityExtended]) -> dict[ActivityExtended, float]:
    with EdgeCache(db_name, method) as cache:
        if MAIN_KEY not in cache.data:
            # Miss
            info(f"Edge LCIA not found for {db_name} / {method}. Building it once.")

            cache.data[MAIN_KEY] = setup_edge_lcia(method, acts[0])

        lcia: EdgeLCIA = cache.data[MAIN_KEY]
        res = dict()
        for act in acts:
            lcia.redo_lcia(demand={act: 1.0})
            res[act] = lcia.score

        return res


# Call once to setup custom serialization
setup_edges_serialization()
