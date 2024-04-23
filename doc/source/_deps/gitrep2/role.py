"""
Sphinx role
"""
from pathlib import Path

from docutils import nodes, utils
from sphinx.util.nodes import split_explicit_title

from .exceptions import ParseError
from .parser import python_to_lineno


def gitref(name, rawtext, text, lineno, inliner, options={}, content=[]):
    """
    Reference a file in git

    Uses the ``gitref_repo`` config value

    Usage:
        :gitref:`target`
        :gitref:`label <target>`

    where ``target`` is a filename with an optional reference to a definition in the
    code:

    * ``path/to/filename.py``
    * ``path/to/filename.py#L6``
    * ``path/to/filename.py::coderef``
    * ``path/to/filename.py::code.ref``

    Renders HTML:
        <a href="url">label</a>
    """
    # Collect config vars
    app = inliner.document.settings.env.app
    try:
        remote = app.config.gitref_remote
        if not remote:
            raise AttributeError
    except AttributeError:
        raise ValueError("Config does not specify gitref_remote")

    # Assume the project root is one dir up from the docs dir
    doc_root = inliner.document.settings.env.srcdir
    project_root = Path(doc_root).parent.parent

    # Process text value
    has_t, title, target = split_explicit_title(text)
    title = utils.unescape(title)
    target = utils.unescape(target)

    # Break target
    if "::" in target:
        filename, coderef = target.split("::", 1)
    else:
        filename = target
        coderef = None

    # Set title if not set
    if title == target:
        if coderef is not None:
            title = app.config.gitref_label_format.format(
                filename=filename,
                coderef=coderef,
            )

    # Ensure the file exists - can be a file or a dir
    filepath = project_root / filename
    ref_start = None
    if not filepath.exists():
        inliner.reporter.error(f"Referenced file does not exist: {filename}")

    # Convert a code ref into a line number
    elif coderef is not None:
        try:
            ref_start = python_to_lineno(filepath, coderef)
        except ParseError as error:
            inliner.reporter.error(
                f'Error resolving code reference "{target}": {error}',
                line=lineno,
            )
            target = filename

    ref = remote.get_url(filename=filename, line=ref_start)

    node = nodes.reference(rawtext, title, refuri=ref, **options)
    return [node], []
