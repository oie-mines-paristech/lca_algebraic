from lca_algebraic import newFloatParam, resetParams
from lca_algebraic.cache import ExprCache, clear_caches


def test_cached_param_should_bind_to_live_param():
    """Check that parameters are bound to unique param values"""

    KEY = "key"

    with ExprCache() as cache:
        p1 = newFloatParam("p1", default=1.0, min=0.0, max=3.0)
        cache.data[KEY] = 2 * p1

    resetParams()

    # Clear in memory cache, but not on disk (pickle
    clear_caches(disk=False)

    # Override p1. Default value is 2 now.
    p1_bis = newFloatParam("p1", default=2.0, min=0.0, max=3.0)

    with ExprCache() as cache:
        expr = cache.data[KEY]

    p1_from_disk = list(expr.free_symbols)[0]

    assert id(p1_from_disk) == id(p1_bis)
    assert p1_from_disk.default == p1_bis.default
