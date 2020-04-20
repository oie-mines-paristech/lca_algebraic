import builtins
import math
import numpy as np
from typing import Dict, List, Union, Tuple

import brightway2 as bw
from tabulate import tabulate
from IPython.core.display import HTML
from bw2data.parameters import ActivityParameter, ProjectParameter, DatabaseParameter, Group
from scipy.stats import triang, truncnorm
from sympy import Symbol

from .base_utils import _eprint, as_np_array

DEFAULT_PARAM_GROUP = "acv"


def _param_registry():
    # Prevent reset upon auto reload in jupyter notebook
    if not 'param_registry' in builtins.__dict__:
        builtins.param_registry = dict()

    return builtins.param_registry


class ParamType:
    '''Type of parameters'''
    ENUM = "enum"
    BOOL = "bool"
    FLOAT = "float"


class DistributionType:
    '''Type of distribution'''
    LINEAR = "linear"
    NORMAL = "normal"
    TRIANGLE = "triangle"
    FIXED = "fixed"


class ParamDef(Symbol):
    '''Generic definition of a parameter, with name, bound, type, distribution
    This definition will serve both to generate brightway2 parameters and to evaluate.

    This class inherits sympy Symbol, making it possible to use in standard arithmetic python
    while keeping it as a symbolic expression (delayed evaluation).
    '''

    def __new__(cls, name, *karg, **kargs):
        return Symbol.__new__(cls, name)

    def __init__(self, name, type: str, default, min=None, max=None, unit="", description="", label=None, label_fr=None,
                 group=None, distrib=DistributionType.LINEAR, std=None):
        self.name = name
        self.type = type
        self.default = default
        self.description = description
        self.min = min
        self.max = max
        self.unit = unit
        self.label = label
        self.label_fr = label_fr
        self.group = group
        self.distrib = distrib

        if type == ParamType.FLOAT and self.min is None:
            self.distrib = DistributionType.FIXED

        if distrib == DistributionType.NORMAL and std is None:
            raise Exception("Standard deviation is mandatory for normal distribution")
        self.std = std

    def label(self):
        if self.label is not None:
            return self.label
        else:
            return self.name.replace("_", " ")

    def range(self, n):
        '''Used for parametric analysis'''
        step = (self.max - self.min) / (n - 1)
        return list(i * step + self.min for i in range(0, n))

    def rand(self, alpha):
        """Transforms a random number between 0 and 1 to valid value according to the distribution of probability of the parameter"""
        if self.distrib == DistributionType.LINEAR:
            return self.min + alpha * (self.max - self.min)

        else :

            if not hasattr(self, "_distrib"):
                if self.distrib == DistributionType.TRIANGLE:
                    scale = self.max - self.min
                    c = (self.default - self.min) / scale
                    self._distrib = triang(c, loc=self.min, scale=scale)

                elif self.distrib == DistributionType.NORMAL:

                    self._distrib = truncnorm(
                        (self.min - self.default) / self.std,
                        (self.max - self.min) / self.std,
                        loc=self.default,
                        scale=self.std)
                else:
                    raise Exception("Unkown distribution type " + self.distrib)


            return self._distrib.ppf(alpha)



    # Expand parameter (useful for enum param)
    def expandParams(self, value=None) -> Dict[str, float]:
        if value == None:
            value = self.default
        return {self.name: value}

    # Useful for enum param, having several names
    def names(self):
        return [self.name]

    def __repr__(self):
        return self.name


class BooleanDef(ParamDef):
    """Parameter with discrete value 0 or 1"""

    def __init__(self, name, **argv):
        super(BooleanDef, self).__init__(name, ParamType.BOOL, min=0, max=1, **argv)

    def range(self, n):
        return [0, 1]

    def rand(self, alpha):
        return np.around(alpha)


class EnumParam(ParamDef):
    """Enum param is a facility representing a choice / switch as many boolean parameters.
    It is not itself a Sympy symbol. use #symbol("value") to access it.
    Statistics weight can be attached to values by providing a dict.
    """

    def __init__(self, name, values: Union[List[str], Dict[str, float]], **argv):
        super(EnumParam, self).__init__(name, ParamType.ENUM, min=None, max=None, **argv)
        if type(values) == list :
            self.values = values
            self.weights = {key:1 for key in values}
        else :
            self.weights = values
            self.values = list(values)
        self.sum = sum(self.weights.values())

    def expandParams(self, currValue=None):
        values = self.values + [None]
        res = dict()
        for enum_val in values:
            var_name = "%s_%s" % (self.name, enum_val if enum_val is not None else "default")
            res[var_name] = 1.0 if enum_val == currValue else 0.0
        return res

    def symbol(self, enumValue):
        """Access parameter for each enum value : <paramName>_<paramValue>"""
        if enumValue is None:
            return Symbol(self.name + '_default')
        if not enumValue in self.values:
            raise Exception("enumValue should be one of %s. Was %s" % (str(self.values), enumValue))
        return Symbol(self.name + '_' + enumValue)

    def names(self):
        return ["%s_%s" % (self.name, value) for value in (self.values + ["default"])]

    def rand(self, alpha):
        alpha = as_np_array(alpha)
        alpha = alpha * self.sum

        # Build bins
        if not hasattr(self, "_bins"):
            self._bins = [0]
            for i in range(len(self.values)) :
                enumvalue = self.values[i]
                self._bins.append(self._bins[i] + self.weights[enumvalue])

        inds = np.digitize(alpha, self._bins, right=True)
        values = np.asarray(self.values)

        return values[inds - 1]

    def range(self, n):
        return self.values


def newParamDef(name, type, **kwargs):
    """Creates a param and register it into a global registry and as a brightway parameter"""
    if type == ParamType.ENUM:
        param = EnumParam(name, **kwargs)
    elif type == ParamType.BOOL:
        param = BooleanDef(name, **kwargs)
    else:
        param = ParamDef(name, type=type, **kwargs)

    # Put it in local registry (in memory)
    if name in _param_registry():
        _eprint("Param %s was already defined : overriding" % name)
    _param_registry()[name] = param

    # Save in brightway2 project
    bwParams = [dict(name=key, amount=value) for key, value in param.expandParams().items()]
    bw.parameters.new_project_parameters(bwParams)

    return param


def newFloatParam(name, default, **kwargs):
    return newParamDef(name, ParamType.FLOAT, default=default, **kwargs)


def newBoolParam(name, default, **kwargs):
    return newParamDef(name, ParamType.BOOL, default=default, **kwargs)


def newEnumParam(name, default, **kwargs):
    return newParamDef(name, ParamType.ENUM, default=default, **kwargs)


def _variable_params(param_names=None):
    if param_names is None :
        param_names =  _param_registry().keys()
    params = {key : _param_registry()[key] for key in param_names}
    return {key: param for key, param in params.items() if param.distrib != DistributionType.FIXED}


def _fixed_params(param_names=None):
    if param_names is None :
        param_names =  _param_registry().keys()
    params = {key : _param_registry()[key] for key in param_names}
    return {key: param for key, param in params.items() if param.distrib == DistributionType.FIXED}


def _listOfDictToDictOflist(LD):
    return {k: [dic[k] for dic in LD] for k in LD[0]}


def _completeParamValues(params):
    """Check parameters and expand enum params.

    Returns
    -------
        Dict of param_name => float value
    """

    # undef_params = param_registry.keys() - params.keys()
    # if undef_params :
    #    raise Exception("Some model parameters are not set : %s" % undef_params)

    res = dict()
    for key, val in params.items():
        if key in _param_registry():
            param = _param_registry()[key]
        else:
            raise Exception("Parameter not found : %s. Valid parameters : %s" % (key, list(_param_registry().keys())))

        if isinstance(val, list):
            newvals = [param.expandParams(val) for val in val]
            res.update(_listOfDictToDictOflist(newvals))
        else:
            res.update(param.expandParams(val))
    return res


def resetParams(db_name):
    """Reset project and activity parameters"""
    _param_registry().clear()
    ProjectParameter.delete().execute()
    ActivityParameter.delete().execute()
    DatabaseParameter.delete().execute()
    Group.delete().execute()


def list_parameters():
    """ Print a pretty list of all defined parameters """
    params = [[param.group, param.label_fr or param.label or param.name, param.default, param.min, param.max, param.unit] for param in
              _param_registry().values()]
    groups = list({p[0] for p in params})
    sorted_params = sorted(params, key=lambda p: groups.index(p[0]))
    return HTML((tabulate(sorted_params, tablefmt="html", headers=["Phase", "param", "default", "min", "max", "unit"])))
