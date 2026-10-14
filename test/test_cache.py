import time

import pytest
from dill import dump

from lca_algebraic import newFloatParam, resetParams, resetDb, newActivity, compute_impacts
from lca_algebraic.cache import EXPR_CACHE, ExprCache, SyncDict, clear_caches
from test.conftest import BG_DB, assert_impacts


def test_cached_param_should_bind_to_live_param(data):
    """Check that parameters are bound to unique param values"""

    KEY = "key"

    with ExprCache(BG_DB) as cache:
        p1 = newFloatParam("p1", default=1.0, min=0.0, max=3.0)
        cache.data[KEY] = 2 * p1

    resetParams()

    # Clear in memory cache, but not on disk (pickle
    clear_caches(disk=False)

    # Override p1. Default value is 2 now.
    p1_bis = newFloatParam("p1", default=2.0, min=0.0, max=3.0)

    with ExprCache(BG_DB) as cache:
        expr = cache.data[KEY]

    p1_from_disk = list(expr.free_symbols)[0]

    assert id(p1_from_disk) == id(p1_bis)
    assert p1_from_disk.default == p1_bis.default


def test_should_invalidate_cache(data):
    # Create a second backgroung db
    BG2 = "bg2"
    resetDb(BG2, True)

    bg_act = newActivity(BG2, "bg_act1", unit="kg", exchanges={data.bio1: 1})

    res = compute_impacts(bg_act, [data.ibio1])
    assert_impacts(res, 1.0)

    # Update same act in background
    newActivity(BG2, "bg_act1", unit="kg", exchanges={data.bio1: 2})

    # Cache should be invcalidated and vamue should change
    res = compute_impacts(bg_act, [data.ibio1])
    assert_impacts(res, 2.0)


@pytest.mark.skipif(not hasattr(time, "tzset"), reason="needs time.tzset")
def test_should_invalidate_cache_across_timezones(data, monkeypatch):
    """bw2data stores 'modified' in naive local time: a change made in another timezone must still clear the cache"""

    def set_tz(tz):
        monkeypatch.setenv("TZ", tz)
        time.tzset()

    BG2 = "bg2"
    resetDb(BG2, True)
    bg_act = newActivity(BG2, "bg_act1", unit="kg", exchanges={data.bio1: 1})

    try:
        set_tz("UTC")
        assert_impacts(compute_impacts(bg_act, [data.ibio1]), 1.0)

        # Changed in New York (UTC-4): 'modified' reads 4 hours old in UTC
        set_tz("America/New_York")
        newActivity(BG2, "bg_act1", unit="kg", exchanges={data.bio1: 2})

        set_tz("UTC")
        assert_impacts(compute_impacts(bg_act, [data.ibio1]), 2.0)
    finally:
        monkeypatch.undo()
        time.tzset()


def test_should_rebuild_old_cache_file(data):
    """A cache file from a previous version (a bare dict) is dropped, not loaded"""

    cache = SyncDict(EXPR_CACHE, BG_DB)
    with open(cache._filename(), "wb") as f:
        dump({"key": 42}, f)

    clear_caches(disk=False)
    with ExprCache(BG_DB) as cache:
        assert "key" not in cache.data
