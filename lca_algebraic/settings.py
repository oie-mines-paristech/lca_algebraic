# This module holds global settings
from contextlib import contextmanager


class Settings:
    cache_enabled: bool = True
    units_enabled: bool = False
    param_overriding_enabled: bool = True

    # If true, upon computation, static background activities would be
    factorize_static_bg: bool = False

    # When activted, forbids any creation or update in background database
    strict_mode: bool = False


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
