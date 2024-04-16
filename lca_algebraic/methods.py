import brightway2 as bw


def findMethods(search=None, mainCat=None):
    """
    Find impact method. Search in all methods against a list of match strings.
    Each parameter can be either an exact match match, or case insenstive search, if suffixed by '*'

    Parameters
    ----------
    search : String to search
    mainCat : if specified, limits the research for method[0] == mainCat.
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


def method_unit(method):
    """Get the unit of an impact method"""
    return bw.Method(method).metadata["unit"]
