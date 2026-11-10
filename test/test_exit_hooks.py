from lca_algebraic import temp_settings
from lca_algebraic.cache import flush_caches_on_exit, install_exit_hooks, _exit_hooks_installed


def test_flush_caches_on_exit_respects_setting(monkeypatch):
    calls = []
    monkeypatch.setattr("lca_algebraic.cache.flush_caches", lambda: calls.append(1))

    with temp_settings(auto_flush_on_exit=True):
        flush_caches_on_exit()
    assert calls == [1]

    # Disabled : no flush
    with temp_settings(auto_flush_on_exit=False):
        flush_caches_on_exit()
    assert calls == [1]


def test_flush_caches_on_exit_swallows_errors(monkeypatch):
    def boom():
        raise RuntimeError("boom")

    monkeypatch.setattr("lca_algebraic.cache.flush_caches", boom)
    with temp_settings(auto_flush_on_exit=True):
        # Must not raise
        flush_caches_on_exit()


def test_install_exit_hooks_is_idempotent():
    install_exit_hooks()
    install_exit_hooks()
    assert _exit_hooks_installed
