import builtins
import math
import numpy as np
from typing import Dict, List, Union, Tuple

import brightway2 as bw
from tabulate import tabulate
from IPython.core.display import HTML
from bw2data.parameters import ActivityParameter, ProjectParameter, DatabaseParameter, Group
from scipy.stats import triang, truncnorm, norm, beta
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
    '''
        Type of statistic distribution of a parameter.
        Some type of distribution requires extra parameters, in italic, to be provided in the constructor of **ParamDef**()

        * **LINEAR** : uniform distribution between *min* and *max*
        * **NORMAL** : Normal distribution, centered on *default* value (mean), with deviation of *std* and truncated between *min* and *max*
        * **TRIANGLE** : Triangle distribution between *min* and *max* (set to zero probability), with highest probability at *default* value
        * **BETA** : Beta distribution with extra params *a* and *b*, using *default* value as 'loc' (0 of beta distribution) and *std* as 'scale' (1 of beta distribution).
                See [scipy doc](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.beta.html#scipy.stats.beta)
    '''
    LINEAR = "linear"
    NORMAL = "normal" # requires 'std' param. 'default' is used as the mean
    BETA="beta" # requires a, b 'default' is used as the mean. 'std' is used as 'scale' factor
    TRIANGLE = "triangle"
    FIXED = "fixed"


class FixedParamMode:
    """ Enum describing what value to set for fixed params """
    DEFAULT = "default"
    MEDIAN = "median"
    MEAN = "mean"

class ParamDef(Symbol):
    '''Generic definition of a parameter, with name, bound, type, distribution
    This definition will serve both to generate brightway2 parameters and to evaluate.

    This class inherits sympy Symbol, making it possible to use in standard arithmetic python
    while keeping it as a symbolic expression (delayed evaluation).
    '''

    def __new__(cls, name, *karg, **kargs):
        return Symbol.__new__(cls, name)

    def __init__(self, name, type: str, default, min=None, max=None, unit="", description="", label=None, label_fr=None,
                 group=None, distrib=None, **kwargs):

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

        # Cleanup distribution in case of overriding already existing param (reused because of inheritance of Symbol)
        if hasattr(self, "_distrib") :
            del self._distrib

        if not distrib :
            if type == ParamType.FLOAT and self.min is None:
                self.distrib = DistributionType.FIXED
            else:
                self.distrib = DistributionType.LINEAR

        elif distrib == DistributionType.NORMAL :
            if not 'std' in kwargs:
                raise Exception("Standard deviation is mandatory for normal distribution")
            self.std = kwargs['std']

        elif distrib == DistributionType.BETA :
            if not 'a' in kwargs or not 'b' in kwargs or not 'std' in kwargs :
                raise Exception("Beta distribution requires params 'a' 'b' and 'std' (used as scale)")
            self.a = kwargs['a']
            self.b = kwargs['b']
            self.std = kwargs['std']

    def stat_value(self, mode : FixedParamMode):
        """Method used to compute fixed statistic value to use for fixed variables"""
        if mode == FixedParamMode.DEFAULT :
            return self.default
        else :
            rnd = np.random.rand(1000)
            x = self.rand(rnd)

            if mode == FixedParamMode.MEAN :
                return np.mean(x)
            elif mode == FixedParamMode.MEDIAN :
                return np.median(x)
            else :
                raise Exception("Unkown mode " + mode)


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

        if self.distrib == DistributionType.FIXED :
            return self.default
        
        elif self.distrib == DistributionType.LINEAR:
            return self.min + alpha * (self.max - self.min)

        else :
            if not hasattr(self, "_distrib"):

                if self.distrib == DistributionType.TRIANGLE:
                    scale = self.max - self.min
                    c = (self.default - self.min) / scale
                    self._distrib = triang(c, loc=self.min, scale=scale)

                elif self.distrib == DistributionType.NORMAL:

                    if self.min :
                        # Truncated normal
                        self._distrib = truncnorm(
                            (self.min - self.default) / self.std,
                            (self.max - self.min) / self.std,
                            loc=self.default,
                            scale=self.std)
                    else :
                        # Normal
                        self._distrib = norm(
                            loc=self.default,
                            scale=self.std)

                elif self.distrib == DistributionType.BETA:
                    self._distrib = beta(
                        self.a,
                        self.b,
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

        # A dict of weights was passed
        if isinstance(currValue, dict) :
            res = { "%s_%s" % (self.name, key) : val / self.sum for key, val in currValue.items()}
            res["%s_default" % self.name] = 0
            return res

        # Normal case
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

    def stat_value(self, mode : FixedParamMode):
        if mode == FixedParamMode.DEFAULT :
            return self.default
        else :
            # For statistical analysis we setup enum as its weights of values,
            # This distrib is then expanded as float parameters, for better fit of the distribution
            return self.weights


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


def _completeParamValues(params, required_params : List[str]=None):
    """Check parameters and expand enum params.

    Returns
    -------
        Dict of param_name => float value
    """

    # Add default values for required params
    if required_params :
        for param_name in required_params :
            param = _param_registry()[param_name]
            if not param_name in params :
                params[param_name] = param.default
                _eprint("Required param '%s' was missing, replacing by default value : %s" % (param_name, str(param.default)))

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
    params = [[param.group, param.label_fr or param.label or param.name, param.default, param.min, param.max, param.std, param.distrib, param.unit] for param in
              _param_registry().values()]
    groups = list({p[0] for p in params})
    sorted_params = sorted(params, key=lambda p: groups.index(p[0]))
    return HTML((tabulate(sorted_params, tablefmt="html", headers=["Phase", "param", "default", "min", "max", "std", "distrib", "unit"])))
