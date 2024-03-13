"""
Remote definitions
"""
import re


if hasattr(re, "Pattern"):
    Pattern = re.Pattern
else:
    # Python 3.6
    Pattern = re._pattern_type


class Registry:
    """
    Remote registry
    """

    def __init__(self):
        self.remotes = []

    def register(self, cls):
        """
        Register the class
        """
        self.remotes.append(cls)
        return cls

    def get_by_url(self, url, branch):
        """
        Return an instantiated Remote for the given git Repo instance
        """
        for remote in self.remotes:
            for pattern in remote.get_remote_patterns():
                matches = pattern.match(url)
                if matches:
                    return remote(matches.group("repo"), branch)
        raise ValueError(f"Unable to find a match for {url}")


registry = Registry()


class Remote:
    """
    Base class for remotes
    """

    #: RegExp: Compiled regex pattern to match the remote url and extract the repo name.
    #: Must contain the named group ``repo``.
    remote_match = None

    #: str: Pattern of URL to code on remote site. Will be formatted with the values
    #: ``repo``, ``branch``, ``filename`` and ``line``
    url_pattern = None

    #: Pattern for line-based URLs. If a line number is specified, this will be rendered
    #: and passed to ``url_pattern`` as ``line``, otherwise ``line`` will be empty.
    url_pattern_line = None

    #: str: The identifying repo name on the remote site. Usually in the format
    #: ``username/repo``.
    repo = None

    #: str: The branch to link to.
    branch = None

    @classmethod
    def __init_subclass__(cls, **kwargs):
        """
        Automatically register Remote subclasses with the registry
        """
        super().__init_subclass__(**kwargs)
        registry.register(cls)

    @classmethod
    def get_remote_patterns(self):
        """
        Return remote_match as a list, in case it's defined as a single regex
        """
        if isinstance(self.remote_match, Pattern):
            return [self.remote_match]
        return self.remote_match

    def __init__(self, repo, branch):
        """
        Arguments:
            repo (str): The identifying repo name on the remote site. Usually in the
                format ``username/repo``.
            branch (str): The branch to link to.
        """
        self.repo = repo
        self.branch = branch

    def get_url(self, filename, line=None):
        """
        Return URL to file.

        Renders ``url_pattern`` and ``url_pattern_line``.
        """
        line_rendered = self.render_line(line)

        url = self.url_pattern.format(
            repo=self.repo, branch=self.branch, filename=filename, line=line_rendered
        )
        return url

    def render_line(self, line):
        if line is None:
            return ""
        return self.url_pattern_line.format(line=line)


class GitHub(Remote):
    """
    Repositories on https://github.com/
    """

    remote_match = [
        re.compile(r"^git@github.com:(?P<repo>.+?)(\.git)?$"),
        re.compile(r"^https://github.com/(?P<repo>.+?)(\.git)?$"),
    ]
    url_pattern = "https://github.com/{repo}/blob/{branch}/{filename}{line}"
    url_pattern_line = "#L{line}"


class Bitbucket(Remote):
    """
    Repositories on https://bitbucket.org/
    """

    remote_match = [
        re.compile(r"^git@bitbucket.org:(?P<repo>.+?)(\.git)?$"),
        re.compile(r"^https://bitbucket.org/(?P<repo>.+?)(\.git)?$"),
    ]
    url_pattern = "https://bitbucket.org/{repo}/src/{branch}/{filename}{line}"
    url_pattern_line = "#lines-{line}"


class GitLab(Remote):
    """
    Repositories on https://gitlab.com/
    """

    remote_match = [
        re.compile(r"^git@gitlab.com:(?P<repo>.+?)(\.git)?$"),
        re.compile(r"^https://gitlab.com/(?P<repo>.+?)(\.git)?$"),
    ]
    url_pattern = "https://gitlab.com/{repo}/blob/{branch}/{filename}{line}"
    url_pattern_line = "#L{line}"

class PrivateGitlab(Remote):

    remote_match = [
        re.compile(r"^https://git.sophia.mines-paristech.fr/.*/(?P<repo>.+?)(\.git)?$"),
    ]
    url_pattern = "https://git.sophia.mines-paristech.fr/oie/{repo}/-/blob/{branch}/{filename}"
    url_pattern_line = "#L{line}"
