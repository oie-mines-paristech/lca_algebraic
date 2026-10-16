import inspect
from collections import Counter
from functools import wraps
from typing import Any, Callable, Dict

from sympy import Add, Basic
from sympy import Dict as SympyDict
from sympy import Expr, Mul, Symbol, Tuple, symbols


def _get_args(expr, only_values=False):
    if isinstance(expr, Tuple):
        return list(expr)
    if isinstance(expr, SympyDict):
        if only_values:
            return expr.values()
        else:
            return expr.items()

    return expr.args


def _walk_unique(root_expr: Basic, callback_f: Callable):
    """Smart walker on non-tree (DAG) sympy expression that doesn't pass twice at the same node"""
    counts = Counter()

    def f_rec(expr: Basic):
        counts[expr] += 1

        if counts[expr] > 1:
            return

        callback_f(expr)

        for arg in _get_args(expr, only_values=True):
            f_rec(arg)

    f_rec(root_expr)
    return counts


def node_stats(root_expr: Basic):
    """Compute stats for each unique expression.
    Returns dicts of (ref_counts, heights, node_number)"""
    counts = Counter()
    heights = Counter()
    numbers = Counter()

    def f_rec(expr: Basic):
        counts[expr] += 1

        if expr not in heights:
            total = 0
            max_height = 0

            # Recurse on args
            for arg in _get_args(expr):
                height, num = f_rec(arg)
                total += num
                max_height = max(height, max_height)

            heights[expr] = max_height + 1
            numbers[expr] = total + 1

        return heights[expr], numbers[expr]

    f_rec(root_expr)

    return counts, heights, numbers


def _normalize_callback(cb: Callable, n: int) -> Callable:
    """Adapts a callback to an given number of argument : drops the extra ones"""
    sig = inspect.signature(cb)
    params = sig.parameters.values()

    # Nombre d'arguments positionnels acceptés
    positional = [
        p
        for p in params
        if p.kind
        in (
            inspect.Parameter.POSITIONAL_ONLY,
            inspect.Parameter.POSITIONAL_OR_KEYWORD,
        )
    ]

    # Si *args est présent, pas besoin de wrapper
    if any(p.kind == inspect.Parameter.VAR_POSITIONAL for p in params):
        return cb

    accepted = len(positional)

    if accepted > n:
        raise TypeError(f"Callback {accepted} arguments, but only {n} are expected")

    @wraps(cb)
    def wrapped(*args):
        return cb(*args[:accepted])

    return wrapped


def _replace_unique(root_expr: Basic, replace_f: Callable[[Any], Basic]):
    """Smart walker that walk nodes only once and possible replace them (with callback)"""

    results = dict()

    replace_f = _normalize_callback(replace_f, 2)

    def f_rec(expr: Basic):
        key = id(expr)
        if key in results:
            return results[key]

        initial_expr = expr

        # Recurse
        args = _get_args(expr)
        if len(args) > 0:
            args = [f_rec(arg) for arg in args]
            expr = expr.func(*args)

        res = replace_f(expr, initial_expr)
        results[key] = res

        return res

    return f_rec(root_expr)


def find_symbols(root_expr: Basic):
    """Find free symbols recursively. Efficient in DAG expressions"""

    symbols = set()

    def _find_symbols(node):
        nonlocal symbols
        if isinstance(node, Symbol):
            symbols.add(node)

    _walk_unique(root_expr, _find_symbols)
    return symbols


def replace_and_cleanup(root_expr: Basic, replacements: Dict[Symbol, Any]) -> Basic:
    """Replace symbols efficiently in DAG expression and possibly reduce expressions (0x..) or sum of zeros"""

    def replace_node(node: Basic) -> Basic:
        if isinstance(node, Symbol) and node in replacements:
            return replacements[node]

        elif isinstance(node, Mul):
            for arg in node.args:
                if arg == 0:
                    return 0
        elif isinstance(node, Add):
            for arg in node.args:
                if arg != 0:
                    return node
            return 0

        return node

    return _replace_unique(root_expr, replace_node)


def count_unique_nodes(expr: Basic):
    _count = 0

    def count_f(node):
        nonlocal _count
        _count += 1

    _walk_unique(expr, count_f)

    return _count


def find_cses(root_expr: Basic, min_nodes=3):
    """CSE (common sub expression) algorithm that just find the existing CSE in the DAG expression."""

    # First pass, coun references
    ref_counts, heights, numbers = node_stats(root_expr)

    # Expression with more than one ref
    candidates = set(expr for expr, count in ref_counts.items() if count > 1 and numbers[expr] > min_nodes)

    cses: Dict[str, Basic] = dict()

    def replacer(node, initial_node):
        if initial_node in candidates:
            name = f"cse{len(cses)}"
            cses[name] = node
            return symbols(name)
        return node

    res_expr = _replace_unique(root_expr, replacer)

    return (list((symbol, expr) for symbol, expr in cses.items()), res_expr)
