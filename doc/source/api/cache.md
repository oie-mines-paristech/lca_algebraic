# Cache

To keep recomputations fast, *lca_algebraic* caches the developed expressions and the
computed LCIA per database. Caches live both in memory and on disk (inside the Brightway
project), and are invalidated automatically when the databases they depend on change.

Caching is controlled by three [settings](settings.md) :

* `cache_enabled` : enable or disable caching entirely.
* `auto_flush` : write caches to disk automatically at the end of a computation (default `True`).
* `auto_flush_on_exit` : flush the caches when the interpreter stops, on uncaught
  exception, or on `SIGINT` / `SIGTERM` (default `True`). This is a safety net when
  `auto_flush` is off, so that a run — or a crash, or a `Ctrl-C` — does not lose the
  work done so far.

In most cases you do not need to touch the cache. The functions below are useful when
you suspect a stale cache, or when you drive the flushing yourself.

## Disable caching

```{eval-rst}
.. autofunction:: lca_algebraic.disable_cache
```

## Flush caches

```{eval-rst}
.. autofunction:: lca_algebraic.flush_caches
```

## Clear caches

```{eval-rst}
.. autofunction:: lca_algebraic.clear_caches
```

## Exit hooks

```{eval-rst}
.. autofunction:: lca_algebraic.install_exit_hooks
```
