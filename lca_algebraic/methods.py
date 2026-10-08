import builtins
from typing import Dict, Tuple

from pint import Unit

from lca_algebraic.bw_wrapper import Method, methods


def _impact_labels():
    """Dictionnary of custom impact names
    Dict of "method tuple" => string
    """
    # Prevent reset upon auto reload in jupyter notebook
    if "_impact_labels" not in builtins.__dict__:
        builtins._impact_labels = dict()

    return builtins._impact_labels


def set_custom_impact_labels(impact_labels: Dict):
    """Global function to override name of impact method in graphs"""
    _impact_labels().update(impact_labels)


def findMethods(search=None, mainCat=None):
    """
    Find impact method. Search in all methods against a list of match strings.
    Each parameter can be either an exact match, or case-insensitive search, if suffixed by '*'

    Parameters
    ----------
    search :
        String to search (optional)
    mainCat :
        If specified, limits the research to methods of this family, e.g. "EF v3.1".
        It is matched against the first or second element of the method tuple, as
        bw2io.import_ecoinvent_release puts the ecoinvent version first.


    Returns
    -------
    A list of tuples, identifying the methods.


    """
    res = []
    search = (search or "").lower()
    for method in methods:
        text = str(method).lower()
        match = search in text
        if mainCat:
            match = match and mainCat in method[:2]
        if match:
            res.append(method)
    return res


def method_unit(method: Tuple, fu_unit: Unit = None):
    """Get the unit of an impact method"""

    res = Method(method).metadata["unit"]
    if fu_unit is not None:
        res += f" / {fu_unit}"

    return res


def method_name(method):
    """Return name of method, taking into account custom label set via set_custom_impact_labels(...)"""
    if method in _impact_labels():
        return _impact_labels()[method]
    # Last two elements : (category, indicator) for both 3 and 4 element tuples
    return " - ".join(str(part) for part in method[-2:])
