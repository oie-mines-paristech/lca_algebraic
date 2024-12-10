import builtins
from collections import defaultdict
from enum import Enum
from typing import Any, Dict, List, Union

import brightway2 as bw
import numpy as np
import pandas as pd
from bw2data.backends.peewee import ExchangeDataset
from bw2data.parameters import (
    ActivityParameter,
    DatabaseParameter,
    Group,
    ProjectParameter,
)
from IPython.core.display import HTML
from pint import Quantity
from scipy.stats import beta, lognorm, norm, triang, truncnorm
from sympy import Basic, Expr, Symbol, lambdify, parse_expr
from tabulate import tabulate

from lca_algebraic.base_utils import ExceptionContext, ValueOrExpression
from lca_algebraic.log import logger

from .base_utils import _snake2camel, _user_functions, as_np_array
from .database import DbContext
from .log import warn
from .settings import Settings
from .units import unit_registry as u

DEFAULT_PARAM_GROUP = "acv"
UNCERTAINTY_TYPE = "uncertainty type"


class ParamType:
    """Type of parameters"""

    ENUM = "enum"
    """ Enum Parameter """

    BOOL = "bool"
    """ Boolean parameter """

    FLOAT = "float"
    """Float parameter """


class DistributionType:
    """
    Type of statistic distribution of a float parameter.
    Some type of distribution requires extra parameters, in italic, to be provided in the constructor of **ParamDef**()
    """

    LINEAR = "linear"
    """ Uniform distribution between *min* and *max*"""

    NORMAL = "normal"
    """ Normal distribution, centered on *default* value (mean), with deviation of *std* and truncated between *min* and *max*"""

    LOGNORMAL = "lognormal"
    """ Lognormal distribution, centered on *default* value (mean), with deviation of *std*, not truncated """

    BETA = "beta"  # requires a, b 'default' is used as the mean. 'std' is used as 'scale' factor
    """ Beta distribution with extra params *a* and *b*,
    using *default* value as 'loc' (0 of beta distribution) and *std* as 'scale' (1 of beta distribution)
    See [scipy doc](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.beta.html#scipy.stats.beta) """

    TRIANGLE = "triangle"
    """ Triangle distribution between *min* and *max* (set to zero probability), with highest probability at *default* value """

    FIXED = "fixed"
    """ Fixed value, not considered as a variable input for monte carlo simulation. """


class _UncertaintyType:
    """Enum of uncertainty types of Brightway.
    See https://stats-arrays.readthedocs.io/en/latest/
    """

    UNDEFINED = 0
    FIXED = 1
    LOGNORMAL = 2
    NORMAL = 3
    UNIFORM = 4
    TRIANGLE = 5
    DISCRETE = 7
    BETA = 10


# Map to Brightay2 / stats-arrays Distribution Type
_DistributionTypeMap = {
    DistributionType.LINEAR: _UncertaintyType.UNIFORM,
    DistributionType.BETA: _UncertaintyType.BETA,
    DistributionType.NORMAL: _UncertaintyType.NORMAL,
    DistributionType.LOGNORMAL: _UncertaintyType.LOGNORMAL,
    DistributionType.TRIANGLE: _UncertaintyType.TRIANGLE,
    DistributionType.FIXED: _UncertaintyType.FIXED,  # I made that up
}

_DistributionTypeMapReverse = {val: key for key, val in _DistributionTypeMap.items()}


class FixedParamMode:
    """Enum describing what value to set for fixed params"""

    DEFAULT = "default"
    """ Use the default value as the parameter """

    MEDIAN = "median"
    """ Use the median of the distribution of the parameter """

    MEAN = "mean"
    """ Use the mean of the distribution of the parameter """


class ParamDef(Symbol):
    """
    Generic definition of a parameter, with name, bound, type, distribution
    This definition will serve both to generate brightway2 parameters and to evaluate.
    This class inherits sympy Symbol, making it possible to use in standard arithmetic python
    while keeping it as a symbolic expression (delayed evaluation).

    Don't instantiate it directly. Use the function **newXXXParam() instead.
    """

    def __new__(cls, name, *karg, **kargs):
        # We use dbname as an "assumption" so that two symbols with same name are not equal if from separate DBs
        assumptions = dict()
        assumptions["real"] = True

        if "dbname" in kargs and kargs["dbname"]:
            assumptions[kargs["dbname"]] = True

        return Symbol.__new__(cls, name, **assumptions)

    def __init__(
        self,
        name,
        type: str,
        default,
        min=None,
        max=None,
        unit="",
        description="",
        label=None,
        group=None,
        distrib: DistributionType = None,
        dbname=None,
        formula=None,
        **kwargs,
    ):
        self.name = name
        self.type = type
        self.default = default
        self.description = description
        self.min = min
        self.max = max
        self.unit = unit
        self.label = label
        self.group = group
        self.distrib = distrib
        self.dbname = dbname
        self.formula = formula

        # if (self.dbname == None) :
        #    error("Warning : param '%s' linked to root project instead of a specific DB" % self.name)

        # Cleanup distribution in case of overriding already existing param (reused because of inheritance of Symbol)
        if hasattr(self, "_distrib"):
            del self._distrib

        if not distrib and type == ParamType.FLOAT:
            if self.min is None:
                raise Exception(f"No 'min/max' provided for {self.name}, distrib should explicitely set to FIXED")
            else:
                self.distrib = DistributionType.LINEAR

        elif distrib in [DistributionType.NORMAL, DistributionType.LOGNORMAL]:
            if "std" not in kwargs:
                raise Exception("Standard deviation is mandatory for normal / lognormal distribution")
            self.std = kwargs["std"]

            if distrib == DistributionType.LOGNORMAL and self.min is not None:
                warn(
                    "Warning : LogNormal does not support min/max boundaries for parameter : ",
                    self.name,
                )

        elif distrib == DistributionType.BETA:
            if "a" not in kwargs or "b" not in kwargs or "std" not in kwargs:
                raise Exception("Beta distribution requires params 'a' 'b' and 'std' (used as scale)")
            self.a = kwargs["a"]
            self.b = kwargs["b"]
            self.std = kwargs["std"]

    def stat_value(self, mode: FixedParamMode):
        """Computes a fixed value for this parameter, either median or mean, according to the requested FixedParamMode and its
        statistical distribution."""
        if mode == FixedParamMode.DEFAULT:
            return self.default
        else:
            # Compute statistical value for replacement
            rnd = np.random.rand(1000)
            x = self.rand(rnd)

            if mode == FixedParamMode.MEAN:
                return np.mean(x)
            elif mode == FixedParamMode.MEDIAN:
                return np.median(x)
            else:
                raise Exception("Unkown mode " + mode)

    def get_label(self):
        if self.label is not None:
            return self.label
        else:
            return self.name.replace("_", " ")

    def range(self, n):
        """Return N uniform values of this parameter, within its range of definition"""
        step = (self.max - self.min) / (n - 1)
        return list(i * step + self.min for i in range(0, n))

    def rand(self, alpha):
        """Transforms a random number between 0 and 1 to valid value according
        to the distribution of probability of the parameter

        Parameters
        ----------

        alpha:
            Can be a *float* or numpy array of floats, between 0 and 1.
        """

        if self.distrib == DistributionType.FIXED:
            return self.default

        elif self.distrib == DistributionType.LINEAR:
            if self.min is None or self.max is None:
                raise Exception("Missing min/max for : " + self.name)
            return self.min + alpha * (self.max - self.min)

        else:
            if not hasattr(self, "_distrib"):
                if self.distrib == DistributionType.TRIANGLE:
                    scale = self.max - self.min
                    c = (self.default - self.min) / scale
                    self._distrib = triang(c, loc=self.min, scale=scale)

                elif self.distrib == DistributionType.NORMAL:
                    if self.min:
                        # Truncated normal
                        self._distrib = truncnorm(
                            (self.min - self.default) / self.std,
                            (self.max - self.min) / self.std,
                            loc=self.default,
                            scale=self.std,
                        )
                    else:
                        # Normal
                        self._distrib = norm(loc=self.default, scale=self.std)

                elif self.distrib == DistributionType.LOGNORMAL:
                    self._distrib = lognorm(self.default, self.std)

                elif self.distrib == DistributionType.BETA:
                    self._distrib = beta(self.a, self.b, loc=self.default, scale=self.std)

                else:
                    raise Exception("Unkown distribution type " + self.distrib)

            return self._distrib.ppf(alpha)

    def __hash__(self):
        return hash((getattr(self, "dbname", ""), self.name))

    def __eq__(self, other):
        if isinstance(other, ParamDef):
            return self.name == other.name and getattr(self, "dbname", None) == getattr(other, "dbname", None)
        else:
            return Symbol.__eq__(self, other)

    def expandParams(self, value=None) -> Dict[str, float]:
        """Abstract function to represent a parameter as a dict : useful to be consistent with Enum params"""
        if value is None:
            value = self.default
        return {self.name: value}

    def names(self, use_label=False):
        """Generic method usefull to be consistent with Enum parameters"""
        if use_label:
            return [self.get_label()]
        else:
            return [self.name]

    def with_unit(self):
        """Returns the symbol together with its unit, as a Pint Quantity"""
        return u.Quantity(self, self.unit)

    def __repr__(self):
        return self.name


class BooleanDef(ParamDef):
    """Parameter with discrete value 0 or 1. Use **newBoolParam()** to instantiate it."""

    def __init__(self, name, **argv):
        if "min" not in argv:
            argv = dict(argv, min=None, max=None)
        super(BooleanDef, self).__init__(name, ParamType.BOOL, **argv)

    def range(self, n):
        return [0, 1]

    def rand(self, alpha):
        return np.around(alpha)


class EnumParam(ParamDef):
    """
    Enum param is a facility wrapping a choice / switch as many boolean parameters.
    It is not itself a Sympy symbol. use #symbol("value") to access it.
    Statistics weight can be attached to values by providing a dict.
    """

    def __init__(self, name, values: Union[List[str], Dict[str, float]], **argv):
        if "min" not in argv:
            argv = dict(argv, min=None, max=None)
        super(EnumParam, self).__init__(name, ParamType.ENUM, **argv)
        if isinstance(values, list):
            self.values = values
            self.weights = {key: 1 for key in values}
        else:
            self.weights = values
            self.values = list(values)
        self.sum = sum(self.weights.values())

    def expandParams(self, currValue=None):
        """
        Return a dictionarry of single enum values as sympy symbols, with only a single one set to 1.
        currValue can be either a single enum value (string) of dict of enum value => weight.
        """

        # A dict of weights was passed
        if isinstance(currValue, dict):
            res = {"%s_%s" % (self.name, key): val / self.sum for key, val in currValue.items()}
            res["%s_default" % self.name] = 0
            return res

        # Normal case
        values = self.values + [None]

        # Bad value ?
        if currValue not in values:
            raise Exception("Invalid value %s for param %s. Should be in %s" % (currValue, self.name, str(self.values)))

        res = dict()
        for enum_val in values:
            var_name = "%s_%s" % (
                self.name,
                enum_val if enum_val is not None else "default",
            )
            res[var_name] = 1.0 if enum_val == currValue else 0.0
        return res

    def symbol(self, choice):
        """Returns the invididual Sympy symbol for a given choice : <paramName>_<choice>"""
        if choice is None:
            return Symbol(self.name + "_default")
        if choice not in self.values:
            raise Exception("enumValue should be one of %s. Was %s" % (str(self.values), choice))
        return Symbol(self.name + "_" + choice)

    def names(self, use_label=False):
        if use_label:
            base_name = self.get_label()
        else:
            base_name = self.name
        return ["%s_%s" % (base_name, value) for value in (self.values + ["default"])]

    def rand(self, alpha):
        alpha = as_np_array(alpha)
        alpha = alpha * self.sum

        # Build bins
        if not hasattr(self, "_bins"):
            self._bins = [0]
            for i in range(len(self.values)):
                enumvalue = self.values[i]
                self._bins.append(self._bins[i] + self.weights[enumvalue])

        inds = np.digitize(alpha, self._bins, right=True)
        values = np.asarray(self.values)

        return values[inds - 1]

    def range(self, n):
        return self.values

    def stat_value(self, mode: FixedParamMode):
        if mode == FixedParamMode.DEFAULT:
            return self.default
        else:
            # For statistical analysis we setup enum as its weights of values,
            # This distrib is then expanded as float parameters, for better fit of the distribution
            return self.weights


def newParamDef(name, type, dbname=None, save=True, **kwargs):
    """
    Creates a parameter and register it into a global registry and as a brightway parameter.

    Parameters
    ----------
    type :
        Type of the parameter (From ParamType)
    save :
        Boolean, persist this into Brightway2 project (True by default)
    dbname :
        Optional name of db. If None, the parameter is a project parameter
    other arguments :
        Refer to the documentation of BooleanDef ParamDef and EnumParam
    """
    if type == ParamType.ENUM:
        param = EnumParam(name, dbname=dbname, **kwargs)
    elif type == ParamType.BOOL:
        param = BooleanDef(name, dbname=dbname, **kwargs)
    else:
        param = ParamDef(name, dbname=dbname, type=type, **kwargs)

    _param_registry()[name] = param

    # Save in brightway2 project
    if save:
        _persistParam(param)

    return param


_BOOLEAN_UNCERTAINTY_ATTRIBUTES = {
    UNCERTAINTY_TYPE: _UncertaintyType.DISCRETE,
    "minimum": 0,
    "maximum": 2,  # upper bound + 1
}


def persistParams():
    """
    Persist parameters into Brightway project, as per : https://stats-arrays.readthedocs.io/en/latest/.

    This is important only in case you want t o:
    - Use parameters outside lca_algebraic (Activity Browser for instance)
    - Use *loadParams()* to init your parameters instead of deleting them and creating them programmaically each time.

    This function is automatically run when creating new parameters. However, it should be called manually after **updating**
    manually some properties of a parameter.
    """
    for param in _param_registry().all():
        _persistParam(param)


def _persistParam(param):
    """Persist parameter into Brightway project"""
    out = []

    # Common attributes for all types of params
    bwParam = dict(
        name=param.name,
        group=param.group,
        label=param.label,
        unit=param.unit,
        formula=str(param.formula) if param.formula else None,
        description=param.description,
    )

    if param.dbname:
        bwParam["database"] = param.dbname

    if param.type == ParamType.ENUM:
        # Enum are not real params but a set of parameters
        for value in param.values:
            enumValueParam = dict(bwParam)
            enumValueParam["name"] = param.name + "_" + value
            enumValueParam.update(_BOOLEAN_UNCERTAINTY_ATTRIBUTES)
            # Use 'scale' as weight for this enum value
            enumValueParam["scale"] = param.weights[value]
            enumValueParam["amount"] = 1 if param.default == value else 0
            out.append(enumValueParam)
    else:
        bwParam["amount"] = param.default

        if param.type == ParamType.BOOL:
            # "Discrete uniform"
            bwParam.update(_BOOLEAN_UNCERTAINTY_ATTRIBUTES)

        elif param.type == ParamType.FLOAT:
            # Save uncertainty
            bwParam[UNCERTAINTY_TYPE] = _DistributionTypeMap[param.distrib]
            bwParam["minimum"] = param.min
            bwParam["maximum"] = param.max
            bwParam["loc"] = param.default

            if param.distrib in [DistributionType.NORMAL, DistributionType.LOGNORMAL]:
                bwParam["scale"] = param.std

            elif param.distrib == DistributionType.BETA:
                bwParam["scale"] = param.std
                bwParam["loc"] = param.a
                bwParam["shape"] = param.b

        else:
            warn("Param type not supported", param.type)

        out.append(bwParam)

    if param.dbname:
        bw.parameters.new_database_parameters(out, param.dbname)
    else:
        bw.parameters.new_project_parameters(out)


def _loadArgs(data):
    """Load persisted data attributes into ParamDef attributes"""
    return {
        "group": data.get("group"),
        "dbname": data.get("database"),
        "default": data.get("amount"),
        "label": data.get("label"),
        "description": data.get("description"),
        "min": data.get("minimum"),
        "max": data.get("maximum"),
        "unit": data.get("unit"),
        "formula": data.get("formula", None),
    }


def loadParams(global_variable=True, dbname=None):
    """
    Load parameters from Brightway database, as per : https://stats-arrays.readthedocs.io/en/latest/

    The recommended way is to delete parameters at the start of your session and to define them programmically each time.

    However, it can be useful if your parameters come from an import or where defined in **Activity Browser**.

    Parameters
    ----------
    global_variable:
        If true, loaded parameters are made available as global variable.
    dbname:
        If provided, load only database parameters for the given db

    Returns
    -------
    A dictionary of parameters
    """

    enumParams = defaultdict(lambda: dict())

    def register(param: ParamDef):
        _param_registry()[param.name] = param

        # Make it available as global var
        if global_variable:
            if param.name in builtins.__dict__:
                warn("Variable '%s' was already defined : overidding it with param." % param.name)

            # If units are activated store param with unit in global variables
            if Settings.units_enabled and param.type == ParamType.FLOAT:
                builtins.__dict__[param.name] = param.with_unit()
            else:
                builtins.__dict__[param.name] = param

    select = DatabaseParameter.select()
    if dbname:
        select = select.where(DatabaseParameter.database == dbname)

    params = list(select)
    if not dbname:
        params += list(ProjectParameter.select())

    for bwParam in params:
        data = bwParam.data
        data["amount"] = bwParam.amount
        name = bwParam.name

        type = data.get(UNCERTAINTY_TYPE, None)

        # print("Data for ", name, data)

        # Common extra args
        args = _loadArgs(data)

        if type == _UncertaintyType.DISCRETE:
            # Boolean or enum

            if data.get("scale") is not None:
                # Enum param : group them by common prefix
                splits = name.split("_")
                enum_value = splits.pop()
                enum_name = "_".join(splits)
                enumParams[enum_name][enum_value] = data
                continue

            elif data["maximum"] == 2:
                del args["max"], args["min"]
                param = newBoolParam(name, save=False, **args)
            else:
                warn(
                    "Non boolean discrete values (max != 2) are not supported for param :",
                    name,
                )
                continue
        else:
            # Float parameter
            if type is None or type == _UncertaintyType.UNDEFINED:
                warn("'Uncertainty type' of param %s not provided. Assuming UNIFORM")
                type = _UncertaintyType.UNIFORM

            # Uncertainty type to distribution type
            args["distrib"] = _DistributionTypeMapReverse[type]

            if type == _UncertaintyType.TRIANGLE:
                args["default"] = data["loc"]

            elif type in [_UncertaintyType.NORMAL, _UncertaintyType.LOGNORMAL]:
                args["default"] = data["loc"]
                args["std"] = data["scale"]

            elif type == _UncertaintyType.BETA:
                args["default"] = data["loc"]
                args["std"] = data["scale"]
                args["a"] = data["loc"]
                args["b"] = data["shape"]

            param = newParamDef(name=name, type=ParamType.FLOAT, save=False, **args)

        # Save it in shared dictionnary
        register(param)

    # Loop on EnumParams
    for param_name, param_values in enumParams.items():
        first_enum_param = list(param_values.values())[0]
        args = _loadArgs(first_enum_param)
        del args["default"]

        # Dictionary of enum values with scale as weight
        args["values"] = {key: data["scale"] for key, data in param_values.items()}

        # Default enum value is the one with amount=1
        defaults = list(key for key, data in param_values.items() if data.get("amount") == 1)
        if len(defaults) == 1:
            default = defaults[0]
        else:
            default = None
            warn("No default enum value found for ", param_name, defaults)

        param = newEnumParam(name=param_name, default=default, save=False, **args)

        # Save it in shared dict
        register(param)

    # Parse formulas
    for param in _param_registry().all():
        with DbContext(param.dbname):
            if isinstance(param.formula, str):
                param.formula = _parse_formula(param.formula)

    return _param_registry()


def newFloatParam(
    name,
    default,
    min: float = None,
    max: float = None,
    unit: str = None,
    description: str = None,
    label: str = None,
    group: str = None,
    distrib: DistributionType = DistributionType.LINEAR,
    formula=None,
    save=True,
    **kwargs,
) -> Union[ParamDef, Quantity]:
    """
    Creates a float (decimal) parameter.

    Parameters
    ----------
    name:
        Name of the parameter
    default:
        Default value
    min:
        Minimum value
    max:
        Maximum value
    unit:
        Unit of the parameter
    description:
        Long description (optional)
    label:
        Extended name (optional)
    group:
        Name of the group (optional). Used to organize parameters.
    distrib:
        Type of the distribution (optional) Linear (uniform) by default
    formula:
        Sympy expression. Optional. If provided the default value of this parameter (if not provided at runtime) will be computed
        from other parameter values.
    kwargs:
        Extra parameters required for advanced distribution types.


    Examples
    --------

    The following code defines a float parameter *p1* of unit *kg*, with a triangle distribution.

    >>> p1 = newFloatParam("p1", min=1.0, max=3.0, default=2.0, distrib=DistributionType.TRIANGLE, unit="kg")

    Returns
    -------
    The newly created parameter

    """

    if Settings.units_enabled and unit is None:
        raise Exception("Unit mode activated : unit is mandatory for parameters")

    param = newParamDef(
        name=name,
        type=ParamType.FLOAT,
        default=default,
        min=min,
        max=max,
        unit=unit,
        description=description,
        label=label,
        group=group,
        distrib=distrib,
        formula=formula,
        save=save,
        **kwargs,
    )

    # If units are enables, wrap float params with their unit
    if Settings.units_enabled:
        return param.with_unit()
    else:
        return param


def newBoolParam(name, default, description: str = None, label: str = None, group: str = None, formula=None, save=True, **kwargs):
    """
    Creates a boolean parameter.

    Parameters
    ----------
    name:
        Name of the parameter
    default:
        Default value
    description:
        Long description (optional)
    label:
        Extended name (optional)
    group:
        Name of the group (optional). Used to organize parameters.
    formula:
        Sympy expression. Optional. If provided the default value of this parameter (if not provided at runtime) will be computed
        from other parameter values.


    Examples
    --------

    >>> p1 = newBoolParam("p1", default=0, group="param group")

    Returns
    -------
    The newly created parameter

    """
    return newParamDef(
        name,
        ParamType.BOOL,
        default=default,
        description=description,
        label=label,
        group=group,
        formula=formula,
        save=save,
        **kwargs,
    )


def newEnumParam(
    name: str,
    default: str,
    values: Union[List[str], Dict[str, float]],
    description: str = None,
    label: str = None,
    group: str = None,
    save=True,
    **kwargs,
):
    """
    Creates an enum parameter : a set of mutually exclusive boolean choices.
    Enum parameters themselves are *not* Sympy symbols. Each of the choice is represented internally
    as a boolean sympy symbol, that can be accessed via **param.symbol("choice_name")**

    Parameters
    ----------

    name:
        Name of the parameter
    default:
        Default
    values:
        Possible choices. The values can be provided as a list of strings, in which case every choice is equiprobable.
        They can also be provided as a python dictionnary of "choice" => value. The proability to be picked is then the
        pro-rata of the value.
    description:
        Long description (optional)
    label:
        Extended name (optional)
    group:
        Name of the group (optional). Used to organize parameters.

    Examples
    --------

    *p1* is an enum param with equiprobable choices "choice_a", "choice_b" and "choice_c"

    >>> p1 = newEnumParam("p1", default="choice_a", values=["choice_a", "choice_b", "choice_c"])

    *p2* is an anum param with "choice_a" and "choice_b" of probability 25% and "choice_c" of probability 50%.

    >>> p2 = newEnumParam("p2", default="choice_a", values={"choice_a":1, "choice_b":1, "choice_c":2})

    """
    return newParamDef(
        name,
        ParamType.ENUM,
        default=default,
        values=values,
        description=description,
        label=label,
        group=group,
        save=save,
        **kwargs,
    )


def _variable_params(param_names=None):
    if param_names is None:
        param_names = _param_registry().keys()
    params = {key: _param_registry()[key] for key in param_names}
    return {key: param for key, param in params.items() if param.distrib != DistributionType.FIXED}


def _fixed_params(param_names=None):
    if param_names is None:
        param_names = _param_registry().keys()
    params = {key: _param_registry()[key] for key in param_names}
    return {key: param for key, param in params.items() if param.distrib == DistributionType.FIXED}


def _listOfDictToDictOflist(LD):
    return {k: [dic[k] for dic in LD] for k in LD[0]}


class DuplicateParamsAndNoContextException(Exception):
    pass


class ParamRegistry:

    """In memory registry of parameters, acting like a dict and maintaining parameters with possibly same names on several DBs"""

    def __init__(self):
        # We store a dict of dict
        # Param Name -> { dbname ->  param}
        self.params: Dict[Dict[ParamDef]] = defaultdict(dict)

    def __len__(self):
        return len(self.params)

    def __getitem__(self, key):
        try:
            params_per_db = self.params[key]

            if len(params_per_db) == 0:
                # Param not found
                raise KeyError("Parameter %s not found" % key)

            if len(params_per_db) == 1:
                return list(params_per_db.values())[0]

            if DbContext.current_db() is None:
                dbs = [key or "<project>" for key in params_per_db.keys()]
                raise DuplicateParamsAndNoContextException(
                    """
                    Found several params with name '%s', linked to databases (%s) . Yet no context is provided.
                    Please embed you code in a DbContext :
                        with DbContext(currentdb) :
                            <code>
                    """
                    % (key, ", ".join(dbs))
                )
            if DbContext.current_db() in params_per_db:
                return params_per_db[DbContext.current_db()]
            else:
                return params_per_db[None]

        except KeyError:
            raise Exception(f"Parameter {key} not found :. Valid parameters : {self.keys()}")

    def as_dict(self):
        return dict(self.items())

    def __setitem__(self, key, param: ParamDef):
        if param.dbname in self.params[key]:
            warn("[ParamRegistry] Param %s was already defined in '%s' : overriding." % (param.name, param.dbname or "<project>"))

        self.params[key][param.dbname] = param

    def __contains__(self, key):
        return key in self.params

    def values(self):
        return [self.__getitem__(key) for key in (self.params)]

    def keys(self):
        return self.params.keys()

    def items(self):
        return [(key, self.__getitem__(key)) for key in self.params.keys()]

    def clear(self, db_name=None):
        if db_name is None:
            self.params.clear()
        else:
            for param_name, db_params in self.params.items():
                if db_name in db_params:
                    del db_params[db_name]

    def all(self):
        """Return list of all parameters, including params with same names and different DB"""
        return list(param for params in self.params.values() for param in params.values())


# Possible param values : either floator string (enum value)
ParamValue = Union[float, str]

# Single value or list of values
ParamValues = Union[List[ParamValue], ParamValue]


def _param_registry() -> ParamRegistry:
    # Prevent reset upon auto reload in jupyter notebook
    if "param_registry" not in builtins.__dict__:
        builtins.param_registry = ParamRegistry()

    return builtins.param_registry


def all_params() -> Dict[str, ParamDef]:
    """Return the dict of all parameters defined in memory"""
    return {param.name: param for param in _param_registry().all()}


def _toSymbolDict(params: Dict[str, Any]):
    """Replace names with actual params as key when possible"""
    all_params = _param_registry().as_dict()
    return {all_params[name] if name in all_params else Symbol(name): val for name, val in params.items()}


def _compute_param_length(params):
    # Check length of parameter values
    param_length = 1
    for key, val in params.items():
        if isinstance(val, (list, np.ndarray)):
            if param_length == 1:
                param_length = len(val)
            elif param_length != len(val):
                raise Exception("Parameters should be a single value or a list of same number of values")
    return param_length


def _expand_params(param_values: Dict[str, ParamValues]):
    res = dict()

    # Expand enum values
    for key, val in list(param_values.items()):
        param = _param_registry()[key]

        if isinstance(val, (list, np.ndarray)):
            newvals = [param.expandParams(val) for val in val]
            res.update(_listOfDictToDictOflist(newvals))
        else:
            res.update(param.expandParams(val))

    # Expand single values to lists of values of same size
    param_length = _compute_param_length(res)
    if param_length > 1:
        for key, val in list(res.items()):
            if not isinstance(val, (list, np.ndarray)):
                val = list([val] * param_length)
            res[key] = np.array(val, float)

    return res


def _complete_params(params: Dict[str, ParamValues], required_params):
    params = params.copy()

    # Add default values for required params
    for param_name in required_params:
        param = _param_registry()[param_name]

        if param_name not in params:
            if param.formula:
                params[param_name] = compute_expr_value(param.formula, params)
                logger.info(f"Param {param_name} was not set. Computing its value from formula :  {params[param_name]}")
            else:
                params[param_name] = param.default
                logger.info("Required param '%s' was missing, replacing by default value : %s" % (param_name, str(param.default)))

    return params


def _complete_and_expand_params(params: Dict[str, ParamValues], required_params: List[str] = None, asSymbols=True):
    """
    Check parameters and expand enum params.
    Also transform single values to list of param values of same size and compute formulas of missing params.

    Returns
    -------
        Dict of param_name => float value or np.arary of float values
    """
    if required_params:
        params = _complete_params(params, required_params)

    # Expand enum values and list of values
    params = _expand_params(params)

    # Replace param name with full symbols
    if asSymbols:
        params = _toSymbolDict(params)

    return params


def resetParams(db_name=None):
    """
    Clear parameters in live memory (registry) and on disk.
    Clear either all params (project and all db params) or db params from a single database (if db_name provided).

    This is a good practice in your code to start fresh, cleaning your foreground database and parameters and redefine all
    programmatically at the start. This ensures the state of the projet / database is always in sync your code and your session
    / in memory.

    """
    _param_registry().clear(db_name)

    if db_name is None:
        ProjectParameter.delete().execute()
        ActivityParameter.delete().execute()
        DatabaseParameter.delete().execute()
    else:
        ActivityParameter.delete().where(ActivityParameter.database == db_name).execute()
        DatabaseParameter.delete().where(DatabaseParameter.database == db_name).execute()
        Group.delete().execute()


class NameType(Enum):
    NAME = "name"
    LABEL = "label"
    CAMEL_NAME = "CAMEL_NAME"


def _param_name(param, name_type: NameType):
    if name_type == NameType.NAME:
        return param.name
    elif name_type == NameType.LABEL:
        return param.get_label()
    else:
        return _snake2camel(param.name)


def list_parameters(name_type=NameType.NAME, as_dataframe=False):
    """Prints a pretty list of all defined parameters

    Parameters
    ----------
    as_dataframe:
        If true, a pandas *Dataframe* is returned. Otherwise, an HTML table is generated.
    """
    params = [
        dict(
            group=param.group or "",
            name=_param_name(param, name_type),
            label=param.get_label(),
            default=param.default,
            min=param.min,
            max=param.max,
            std=getattr(param, "std", None),
            distrib=param.distrib,
            unit=param.unit,
            db=param.dbname or "[project]",
        )
        for param in _param_registry().all()
    ]

    groups = list({p["group"] for p in params})
    groups = sorted(groups)

    # Sort by Group / name
    def keyf(param):
        return (groups.index(param["group"]), param["name"])

    sorted_params = sorted(params, key=keyf)

    if as_dataframe:
        return pd.DataFrame(sorted_params)
    else:
        return HTML(tabulate(sorted_params, tablefmt="html", headers="keys"))


def compute_expr_value(expr: Expr, param_values: Dict):
    """Compute value of an expression for given set of parameter values"""
    free_symbols = [str(symbol) for symbol in expr.free_symbols]
    lambd = lambdify(free_symbols, expr)

    required_params = _expanded_names_to_names(free_symbols)

    values = _complete_and_expand_params(param_values, required_params=required_params, asSymbols=False)

    # Filter only required params
    values = {name: val for name, val in values.items() if name in free_symbols}

    return lambd(**values)
    # return expr.evalf(subs=_completeParamValues(param_values, required_params=required_params))


def freezeParams(db_name, **params: Dict[str, float]):
    """
    Freezes amounts in all exchanges for a given set of parameter values.
    The formulas are computed and the 'amount' attributes are set with the result.

    This enables parametric datasets to be used by standard, non-parametric tools of Brightway2 (like Activities browser).

    Parameters
    ----------
    db_name :
        Name of the database for freeze (your foreground db usually)

    params:
        All other parameters of this function are threated as the values of *lca_algebraic* parameters to be set.
        The default values will be used for the *lca_algebraic* parameters not mentioned here.

    Examples
    --------

    >>> freezeParams("USER_DB", p1=0.1, p2=3.0)

    """

    db = bw.Database(db_name)

    with DbContext(db):
        for act in db:
            for exc in act.exchanges():
                amount = _getAmountOrFormula(exc)

                # Amount is a formula ?
                if isinstance(amount, Expr):
                    val = compute_expr_value(amount, params)

                    with ExceptionContext(val):
                        val = float(val)

                    print("Freezing %s // %s : %s => %0.2f" % (act, exc["name"], amount, val))

                    # Update in DB
                    exc["amount"] = val
                    exc.save()


def _listParams(db_name) -> List[ParamDef]:
    """
    Return a set of all parameters used in activities
    """

    db = bw.Database(db_name)
    res = set()

    with DbContext(db):
        for act in db:
            for exc in act.exchanges():
                amount = _getAmountOrFormula(exc)

                # Amount is a formula ?
                if isinstance(amount, Basic):
                    expanded_names = list(str(symbol) for symbol in amount.free_symbols)
                    param_names = _expanded_names_to_names(expanded_names)
                    params = list(_param_registry()[param_name] for param_name in param_names)
                    res.update(params)
    return res


def _expand_param_names(param_names: List[str]) -> List[str]:
    """Expand parameters names (with enum params)"""
    return [name for key in param_names for name in _param_registry()[key].names()]


def _expanded_names_to_names(param_names):
    """Find params corresponding to expanded names, including enums."""
    param_names = set(param_names)

    # param name => param
    res = dict()

    # Search for param with same name of prefix paramName_enumValue
    for expended_name in param_names:
        for param_name in _param_registry().keys():
            if expended_name.startswith(param_name):
                param = _param_registry()[param_name]
                for name in param.names():
                    if name == expended_name:
                        res[expended_name] = param

    missing = param_names - set(res.keys())
    if len(missing) > 0:
        raise Exception("Unkown params : %s" % missing)

    return {param.name for param in res.values()}


def _parse_formula(formula):
    local_dict = {x[0].name: x[0] for x in _user_functions.values()}
    return parse_expr(formula, local_dict=local_dict | _param_registry().as_dict())


def _getAmountOrFormula(ex: ExchangeDataset) -> Union[Basic, float]:
    """Return either a fixed float value or an expression for the amount of this exchange"""
    if "formula" in ex:
        try:
            # We don't want support for units there
            return _parse_formula(ex["formula"])
        except Exception as e:
            warn(f"Error '{e}' while parsing formula {ex['formula']} : backing to amount")

    return ex["amount"]


def switchValue(param: EnumParam, **values: Dict[str, ValueOrExpression]):
    """

    Helper method defining an expression that returns a different value / formula for each possible choice of an anum param.

    Parameters
    ----------
    param: EnumParam
        The enum param

    values: Dict[str, ValueOrExpression]
        Each param should correspond to a valid choice of the num parameter.


    Examples
    --------

    Given the enum parameter *p1* :

    >>> p1 = newEnumParam("p1", values=["choice1", "choice2", "choice3"])

    The following code defines an expression worth 0.1 for *choice1*, 0.2 for *choice2*  and *4 x p2* for *choice3*

    >>> amount = switchValue(p1, choice1=0.1, choice2=0.2, choice3=4*p2)

    """

    res = 0
    for key, val in values.items():
        res += param.symbol(key) * val
    return res
