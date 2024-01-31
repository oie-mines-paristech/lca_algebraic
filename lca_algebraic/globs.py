# Global cache, put in separate file to prevent clenup upon autoreload

# Cache of (act, method) => values
_BG_IMPACTS_CACHE = dict()

def _clearLCACache() :
    _BG_IMPACTS_CACHE.clear()