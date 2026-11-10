from lca_algebraic import compute_impacts, newActivity, resetDb, temp_settings
from lca_algebraic.bw_wrapper import databases
from lca_algebraic.cache import get_last_update
from lca_algebraic.database import _setMeta
from lca_algebraic.settings import PROXY_DB_FLAG
from test.conftest import USER_DB, assert_impacts


def test_get_last_update_tracks_proxy_db(data):
    """A proxy db modification must be visible through get_last_update.

    Proxy dbs are only materialized when background factorization is enabled :
    in that mode the cache of the proxy db must be invalidated when it changes,
    otherwise stale proxy impacts are silently reused (see _CacheDict).
    """
    resetDb("fg-proxy")
    _setMeta("fg-proxy", PROXY_DB_FLAG, True)

    with temp_settings(factorize_static_bg=True):
        assert get_last_update("fg-proxy") > 0


def test_proxy_cache_invalidated_on_update(data):
    """With factorization, updating a static foreground amount must update the
    impact, even though the value is served from the (rebuilt) proxy db cache."""

    with temp_settings(factorize_static_bg=True):
        act = newActivity(USER_DB, "act", "kg", {data.bg_act1: 2})

        assert_impacts(compute_impacts(act, [data.ibio1]), 2.0)

        # Factorization materialized a proxy db for the background amount
        assert "fg-proxy" in databases

        act.updateExchanges({data.bg_act1: 3})

        assert_impacts(compute_impacts(act, [data.ibio1]), 3.0)


def test_no_proxy_db_without_factorization(data):
    """Caches are only proxy-dependent when factorization is enabled."""

    with temp_settings(factorize_static_bg=False):
        act = newActivity(USER_DB, "act", "kg", {data.bg_act1: 2})

        assert_impacts(compute_impacts(act, [data.ibio1]), 2.0)
        assert "fg-proxy" not in databases
