# Settings

*lca_algebraic* exposes a set of global options as class attributes of the `Settings`
class. Set them at the beginning of your script / notebook :

```python
from lca_algebraic import Settings

Settings.units_enabled = True
```

Each setting can also be overridden temporarily, for the duration of a `with` block,
using `temp_settings()`.

## Available settings

| Setting | Default | Description |
|---|---|---|
| `units_enabled` | `False` | Enable automatic checks and conversion of physical units (see [Units](units.md)). |
| `cache_enabled` | `True` | Enable the on-disk / in-memory caches of LCIA and developed expressions. Disable it if you suspect a stale cache (see [Cache](cache.md)). |
| `auto_flush` | `True` | If `True`, caches are written automatically when a computation finishes. If `False`, you must flush them yourself with `flush_caches()`. |
| `auto_flush_on_exit` | `True` | Install hooks that flush the caches on interpreter exit, uncaught exception and `SIGINT` / `SIGTERM`. Safety net when `auto_flush` is off. |
| `strict_mode` | `False` | Forbid any creation or update of activities in a **background** database. |
| `factorize_static_bg` | `False` | Factorize static background exchanges during impact computation. Improves performance on large models, at the cost of a more complex internal representation. |
| `param_overriding_enabled` | `True` | Allow redefining an existing parameter. If `False`, redefining a parameter raises an exception. |
| `scenario_separator` | `"#"` | Character used in database names to separate the base name from a scenario name. |
| `cse_mode` | `CSE.DEFAULT` | Common-subexpression-elimination mode used when lambdifying expressions (see below). |

## Common subexpression elimination mode

`Settings.cse_mode` selects how sympy eliminates recurring sub-expressions when an
expression is compiled to a numpy function. Use the `CSE` enum :

* `CSE.NONE` : disables CSE.
* `CSE.DEFAULT` : use sympy's default CSE pass (default value).
* `CSE.DAG_REFS` : reuse the existing references of the expression DAG instead of
  recomputing common subexpressions.

```python
from lca_algebraic import Settings
from lca_algebraic.settings import CSE

Settings.cse_mode = CSE.NONE
```

## Temporary override

```{eval-rst}
.. autofunction:: lca_algebraic.temp_settings
```
