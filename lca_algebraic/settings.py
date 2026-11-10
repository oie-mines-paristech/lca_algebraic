# This module holds global settings
from contextlib import contextmanager

try:
    from enum import Enum, StrEnum
except ImportError:
    # Simulate StrEnum for python<3.11
    from enum import Enum

    class StrEnum(str, Enum):
        pass


# Custom method that returns "1.0" used for inventories
UNITY_METHOD = ("UNITY", "UNITY", "UNITY")


class CSE(StrEnum):
    NONE = "none"
    DEFAULT = "default"  # Default sympy one
    DAG_REFS = "dag"  # Use existing refs in DAG pexresions


class Settings:
    cache_enabled: bool = True
    units_enabled: bool = False
    param_overriding_enabled: bool = True

    # How to recognize an lca_algebraic symbolic amount in an exchange.
    # The "formula" field is shared : lca_algebraic stores its expressions there
    # on FOREGROUND exchanges, while ecoinvent stores chemical formulas in it on
    # background exchanges.
    #   True  (default) : a "formula" is ours only when tagged with
    #                     "has_lca_algebraic_formula" (introduced in v1.4.1.devX).
    #   False           : legacy detection -- on a FOREGROUND exchange, any present
    #                     "formula" is an lca_algebraic expression; background
    #                     formulas stay ignored. Required to read databases written
    #                     before the tag existed (untagged formulas otherwise fall
    #                     back to "amount", which is 0 for foreground).
    use_lca_formula_tag: bool = True

    # If true, upon computation, static background activities would be
    factorize_static_bg: bool = False

    # When activted, forbids any creation or update in background database
    strict_mode: bool = False

    # Set temporarly to true for internal tools not to trigger warnings
    internals = False

    # Use "Common Sub expression patterns" in sympy lambdify
    cse_mode = CSE.DEFAULT

    # If true, caches are flushed automatically but more frequently
    # If set to flase, you shoudl ensure flushing caches yourself
    auto_flush = True

    # If true (default), install hooks that flush the caches on interpreter exit,
    # on uncaught exception, and on SIGINT/SIGTERM. Safety net so a run leaves no
    # lost work when auto_flush is off (caches are otherwise only written on an
    # explicit flush_caches()).
    auto_flush_on_exit = True

    # Separator considered in database name for scenarios
    # You may chose another one with several safer catacters like "__" or "-scenario-"
    scenario_separator = "#"


# Flag used on Db to set it as proxy
PROXY_DB_FLAG = "isProxy"


@contextmanager
def temp_settings(**kwargs):
    """Temporary switch the value of a setting. Argument can be one of the Settings.
    For instance `with settings(units_enabled=False) : ...`"""

    # Ensure keys are corect
    for key in kwargs.keys():
        if not hasattr(Settings, key):
            raise Exception(f"Invalid setting name {key}")

    old_values = {key: getattr(Settings, key) for key in kwargs}

    try:
        for key, val in kwargs.items():
            setattr(Settings, key, val)
        yield Settings
    finally:
        # Reset old values
        for key, val in old_values.items():
            setattr(Settings, key, val)
