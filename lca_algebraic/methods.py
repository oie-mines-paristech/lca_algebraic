import builtins
from typing import Dict, Tuple

import brightway2 as bw
from pint import Unit


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
        String to search
    mainCat :
        If specified, limits the research for method[0] == mainCat.


    Returns
    -------
    A list of tuples, identifying the methods.


    """
    res = []
    search = search.lower()
    for method in bw.methods:
        text = str(method).lower()
        match = search in text
        if mainCat:
            match = match and (mainCat == method[0])
        if match:
            res.append(method)
    return res


def method_unit(method: Tuple, fu_unit: Unit = None):
    """Get the unit of an impact method"""

    res = bw.Method(method).metadata["unit"]

    if fu_unit is not None:
        res += f" / {fu_unit}"

    return res


def method_name(method):
    """Return name of method, taking into account custom label set via set_custom_impact_labels(...)"""
    if method in _impact_labels():
        return _impact_labels()[method]
    return method[1] + " - " + method[2]
