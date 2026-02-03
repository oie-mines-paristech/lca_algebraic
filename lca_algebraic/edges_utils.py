import json
from functools import cache
from os import path
from typing import Dict

from edges import EdgeLCIA, get_available_methods
from edges.filesystem_constants import DATA_DIR
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

    def post_load(lcia: EdgeLCIA):
        if lcia.logger is None:
            lcia.logger = logger
        if lcia.lca:
            if lcia.lca.logger is None:
                lcia.lca.logger = logger
            if lcia.lca.solver is None:
                lcia.lca.solver = factorized(lcia.lca.technosphere_matrix.tocsc())
        if lcia._geo and lcia._geo.logger is None:
            lcia._geo.logger = logger

    def get_state(lcia: EdgeLCIA):
        lcia.logger = None
        if lcia.lca:
            lcia.lca.logger = None
        if lcia._geo:
            lcia._geo.logger = None
        if lcia.lca.solver:
            lcia.lca.solver = None
        try:
            return lcia
        finally:
            # Executed even after the return
            post_load(lcia)

    def set_state(lcia, state):
        lcia.__dict__.update(state)
        post_load(lcia)

    EdgeLCIA.__getstate__ = get_state
    EdgeLCIA.__setstate__ = set_state


class EdgeCache(_CacheDict):
    """Custom cache for EdgeLCIA. We use one separate cache per method, because each pickled file is big (100Mb)"""

    def __init__(self, db_name, method: tuple):
        key = "_".join(item for item in method)
        super().__init__(f"edge_lcia_{key}", db_name)


def setup_edge_lcia(method_key, act: ActivityExtended):
    """Done once then cached"""

    info(f"Edge LCIA not found for {method_key}. Building it once.")

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
            cache.data[MAIN_KEY] = setup_edge_lcia(method, acts[0])

        lcia: EdgeLCIA = cache.data[MAIN_KEY]
        res = dict()
        for act in acts:
            lcia.redo_lcia(demand={act: 1.0})
            res[act] = lcia.score

        return res


# Call once to setup custom serialization
setup_edges_serialization()
