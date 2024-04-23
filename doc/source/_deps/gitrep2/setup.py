"""
Set up the extension in sphinx

The ``setup()`` function will be called automatically by Sphinx
"""
from pathlib import Path

from .constants import DEFAULT_LABEL_FORMAT
from .git import Repo
from .remote import registry
from .role import gitref


def lookup_remote(app, config):
    """
    Once configuration is complete, build the Remote
    """
    try:
        remote_url = config.gitref_remote_url
        if not remote_url:
            raise AttributeError
    except AttributeError:
        raise ValueError("Could not determine gitref_remote_url, must set explicitly")

    try:
        branch = config.gitref_branch
        if not branch:
            raise AttributeError
    except AttributeError:
        raise ValueError("Could not determine gitref_branch, must set explicitly")

    config.gitref_remote = registry.get_by_url(remote_url, branch)


def setup(app):
    """
    Prepare config values and register role
    """
    # Pick up config defaults from current repo
    doc_root = Path(app.confdir)
    git_root = doc_root.parent / ".git"
    repo = Repo(git_root)

    # Add config variables
    app.add_config_value("gitref_remote_url", repo.get_remote_url(), "html")
    app.add_config_value("gitref_branch", repo.get_local_branch(), "html")
    app.add_config_value("gitref_label_format", DEFAULT_LABEL_FORMAT, "html")

    # Listen for config initialised
    app.connect("config-inited", lookup_remote)

    # Register role
    app.add_role("gitref", gitref)
