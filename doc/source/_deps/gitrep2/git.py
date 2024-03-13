"""
Working with git

This is a lightweight abstraction layer to provide access to necessary data from the
current git repository. At some point we may want the reliability improvement and
additional features available by moving to something like GitPython; that should be
possible with a drop-in replacement, but for now this should suffice.
"""
import re
from configparser import NoSectionError, RawConfigParser
from io import StringIO

from .constants import DEFAULT_REMOTE


class Repo:
    """
    Represent a local git repository
    """

    #: Path to ``.git`` repo
    path = None

    #: Name of the remote in ``.git/config``
    remote_name = None

    BRANCH_PATTERN = re.compile(r"^ref: refs/heads/(?P<branch>.+)$")

    def __init__(self, path, remote_name=DEFAULT_REMOTE):
        """
        Initialise a repo wrapper

        Args:
            path (Path): The path to the .git repo
        """
        if not path.is_dir():
            # If the pass doesn't exist then leave path as None so we can skip
            # unnecessary checks
            return

        self.path = path
        self.remote_name = remote_name

    def get_remote_url(self):
        """
        Find the URL for the current remote
        """
        if self.path is None:
            return None

        config_path = self.path / "config"
        if not config_path.is_file():
            return None

        # Read the config file and strip each line (configparser doesn't like tabs)
        config_src = "\n".join(
            [line.strip() for line in config_path.read_text().splitlines()]
        )
        config_io = StringIO(config_src)

        # Read the config - intentionally don't catch exceptions, we need those reported
        config = RawConfigParser(allow_no_value=True)
        config.readfp(config_io)

        try:
            url = config.get(f'remote "{self.remote_name}"', "url")
        except NoSectionError:
            return None

        return url

    def get_local_branch(self):
        """
        The current branch name
        """
        if self.path is None:
            return None

        # Get the current HEAD
        git_head_path = self.path / "HEAD"
        if not git_head_path.is_file():
            return None
        git_head = git_head_path.read_text().strip()

        # See if we're on a named branch
        matches = self.BRANCH_PATTERN.match(git_head)
        if matches:
            return matches.group("branch")

        # We're on a detached head, likely readthedocs; try to look it up in packed-refs
        packed_refs_path = self.path / "packed-refs"
        if not packed_refs_path.is_file():
            return None

        packed_refs_raw = packed_refs_path.read_text()
        remote_ref_prefix = f"refs/remotes/{self.remote_name}/"
        for line in packed_refs_raw.splitlines():
            # Skip irrelevant lines
            if line.startswith("#") or " " not in line:
                continue

            # Typical line:
            #   1234567890abcdef refs/remotes/origin/develop
            commit_hash, ref = line.split(" ", 1)
            if commit_hash == git_head and ref.startswith(remote_ref_prefix):
                ref = ref[len(remote_ref_prefix) :]
                return ref

        return None
