"""
Parse python source
"""
import ast

from .exceptions import NodeNotFound, ParseError


# Prep assignment node types:
#   ast.Assign: assignment, eg ``x = 1``
#   ast.AugAssign: augmented assignment, eg ``x += 1``
#   ast.AnnAssign: annotated assignment, eg ``x: int = 1``
AssignTypes = (ast.Assign, ast.AugAssign, ast.AnnAssign)


def find_name_in_nodes(name, nodes, path=""):
    """
    Given a list of AST nodes, find the named object

    Arguments:
        name: object name, or Python path to object name

    Returns:
        node: AST node that defines the name

    Raises:
        NodeNotFound: no AST node was found matching that name
    """
    # Look for a path
    if "." in name:
        name, rest = name.split(".", 1)
    else:
        rest = ""

    # Track where we came from for reporting purposes
    if not path:
        path = name
    else:
        path = f"{path}.{name}"

    for node in nodes:
        # Check for assignments
        if type(node) in AssignTypes:
            for target in node.targets:
                if not hasattr(target, "id") or target.id != name:
                    continue

                if rest:
                    # Found the object but don't know how to peek inside assignments
                    raise NodeNotFound(f'Found "{path}" but cannot go further')

                return node
            # Didn't find it in that assignment

        # Named code blocks
        elif type(node) in (ast.FunctionDef, ast.ClassDef):
            if node.name != name:
                continue

            # Found the object
            if not rest:
                return node

            # More to go, peek inside
            return find_name_in_nodes(rest, node.body, path)

    # Exhausted all nodes
    raise NodeNotFound(f'Couldn\'t find "{path}"')


def python_to_lineno(path, ref):
    """
    Look up a code reference in the specific python file and return the line number it
    is defined on

    Args:
        path (Path): Path to python file to examine
        ref (str): Variable, function or class name to find

    Returns:
        start (int): Line number the reference starts on
        end (int): Line number the reference ends on

    Raises:
        ParseError: If unable to resovle the reference for some reason
    """
    if path.suffix != ".py":
        raise ParseError("Source file is not Python")

    # Examine the file
    raw = path.read_text()

    return python_string_to_lineno(raw, ref)


def python_string_to_lineno(raw, ref):
    """
    Look up a code reference in the provided python source and return the line number it
    is defined on

    Args:
        filename (str): Path to python file to examine
        ref (str): Variable, function or class name to find

    Returns:
        start (int): Line number the reference starts on

    Raises:
        ParseError: If unable to resolve the reference for some reason
    """
    try:
        module = ast.parse(raw)
    except Exception as e:
        raise ParseError(str(e))

    try:
        node = find_name_in_nodes(ref, module.body)
    except NodeNotFound as e:
        raise ParseError(str(e))

    return node.lineno
