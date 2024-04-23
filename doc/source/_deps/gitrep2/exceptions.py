"""
Exceptions
"""


class ParseError(Exception):
    """
    The parser was unable to complete its work
    """


class NodeNotFound(Exception):
    """
    Node not found when walking AST
    """
