import builtins
import enum
from collections import defaultdict
from enum import Enum
from typing import Dict, List, Union, Tuple

import brightway2 as bw
import ipywidgets as widgets
import numpy as np
from IPython.core.display import HTML
from bw2data.backends import LCIBackend
from bw2data.backends.peewee import Activity
from bw2data.database import Database
from bw2data.parameters import ActivityParameter, ProjectParameter, DatabaseParameter, Group
from bw2data.proxies import ActivityProxyBase
from lca_algebraic.base_utils import ExceptionContext
from scipy.stats import triang, truncnorm, norm, beta, lognorm
from sympy import Symbol, Basic
from tabulate import tabulate

from .base_utils import _snake2camel
from .base_utils import error, as_np_array, _getAmountOrFormula, LANG

import lca_algebraic.base_utils
from .base_utils import as_np_array, _getAmountOrFormula, _actName

DEFAULT_PARAM_GROUP = "acv"
UNCERTAINTY_TYPE = "uncertainty type"




class DbContext :
    """
        Context class specifying the current foreground DB in use. in internal
        Used internally to distinguish database parameters with same names

        usage :
        with DbContext("db") :
            <some code>

    """
    stack = []

    @staticmethod
    def current_db() :
        if len(DbContext.stack) == 0:
            return None
        return DbContext.stack[-1]


    def __init__(self, db : Union[str, ActivityProxyBase, LCIBackend]) :
        if isinstance(db, ActivityProxyBase) :
            self.db = db.key[0]
        elif isinstance(db, str) :
            self.db = db
        else :
            self.db = db.name

    def __enter__(self):
        DbContext.stack.append(self.db)

    def __exit__(self, exc_type, exc_value, exc_traceback):
        DbContext.stack.pop()





class ParamType:
    """Type of parameters"""

    ENUM = "enum"
    """ Enum Parameter """

    BOOL = "bool"
    """ Boolean parameter """

    FLOAT = "float"
    """Float parameter """

class DistributionType():
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

    BETA = "beta" # requires a, b 'default' is used as the mean. 'std' is used as 'scale' factor
    """ Beta distribution with extra params *a* and *b*, 
    using *default* value as 'loc' (0 of beta distribution) and *std* as 'scale' (1 of beta distribution)
    See [scipy doc](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.beta.html#scipy.stats.beta) """

    TRIANGLE = "triangle"
    """ Triangle distribution between *min* and *max* (set to zero probability), with highest probability at *default* value """

    FIXED = "fixed"
    """ Fixed value, not considered as a variable input for monte carlo simulation. """


class _UncertaintyType :
    '''Enum of uncertainty types of Brightway.
    See https://stats-arrays.readthedocs.io/en/latest/
    '''
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
    DistributionType.LINEAR : _UncertaintyType.UNIFORM,
    DistributionType.BETA : _UncertaintyType.BETA,
    DistributionType.NORMAL : _UncertaintyType.NORMAL,
    DistributionType.LOGNORMAL: _UncertaintyType.LOGNORMAL,
    DistributionType.TRIANGLE : _UncertaintyType.TRIANGLE,
    DistributionType.FIXED : _UncertaintyType.FIXED, # I made that up
}

_DistributionTypeMapReverse = {val:key for key, val in _DistributionTypeMap.items()}


class FixedParamMode:
    """ Enum describing what value to set for fixed params """
    DEFAULT = "default"
    """ Use the default value as the parameter """

    MEDIAN = "median"
    """ Use the median of the distribution of the parameter """

    MEAN = "mean"
    """ Use the mean of the distribution of the parameter """

class ParamDef(Symbol):
    '''Generic definition of a parameter, with name, bound, type, distribution
    This definition will serve both to generate brightway2 parameters and to evaluate.

    This class inherits sympy Symbol, making it possible to use in standard arithmetic python
    while keeping it as a symbolic expression (delayed evaluation).
    '''

    def __new__(cls, name, *karg, **kargs):
        # We use dbname as an "assumption" so that two symbols with same name are not equal if from separate DBs
        assumptions = dict()
        if 'dbname' in kargs and kargs['dbname'] :
            assumptions[kargs['dbname']] = True

        return Symbol.__new__(cls, name, **assumptions)

    def __init__(self, name, type: str, default, min=None, max=None, unit="", description="", label=None, label_fr=None,
                 group=None, distrib:DistributionType=None, dbname=None, **kwargs):

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
        self.dbname = dbname

        if (self.dbname == None) :
            error("Warning : param '%s' linked to root project instead of a specific DB" % self.name)

        # Cleanup distribution in case of overriding already existing param (reused because of inheritance of Symbol)
        if hasattr(self, "_distrib") :
            del self._distrib

        if not distrib and type == ParamType.FLOAT  :
            if self.min is None:
                error("No 'min/max' provided, param %s marked as FIXED" % self.name)
                self.distrib = DistributionType.FIXED
            else:
                self.distrib = DistributionType.LINEAR

        elif distrib in [DistributionType.NORMAL, DistributionType.LOGNORMAL] :
            if not 'std' in kwargs:
                raise Exception("Standard deviation is mandatory for normal / lognormal distribution")
            self.std = kwargs['std']

            if distrib == DistributionType.LOGNORMAL and self.min is not None :
                error("Warning : LogNormal does not support min/max boundaries for parameter : ", self.name)

        elif distrib == DistributionType.BETA :
            if not 'a' in kwargs or not 'b' in kwargs or not 'std' in kwargs :
                raise Exception("Beta distribution requires params 'a' 'b' and 'std' (used as scale)")
            self.a = kwargs['a']
            self.b = kwargs['b']
            self.std = kwargs['std']

    def stat_value(self, mode : FixedParamMode):
        """ Compute a fixed value for this parameter, according to the requested FixedParamMode """
        if mode == FixedParamMode.DEFAULT :
            return self.default
        else :
            # Compute statistical value for replacement
            rnd = np.random.rand(1000)
            x = self.rand(rnd)

            if mode == FixedParamMode.MEAN :
                return np.mean(x)
            elif mode == FixedParamMode.MEDIAN :
                return np.median(x)
            else :
                raise Exception("Unkown mode " + mode)


    def get_label(self):
        if LANG == "fr " and self.label_fr is not None :
            return self.label_fr
        elif self.label is not None:
            return self.label
        else:
            return self.name.replace("_", " ")

    def range(self, n):
        """Return N uniform values of this parameter, within its range of definition"""
        step = (self.max - self.min) / (n - 1)
        return list(i * step + self.min for i in range(0, n))

    def rand(self, alpha):
        """Transforms a random number between 0 and 1 to valid value according to the distribution of probability of the parameter"""

        if self.distrib == DistributionType.FIXED :
            return self.default
        
        elif self.distrib == DistributionType.LINEAR:
            if self.min is None or self.max is None :
                raise Exception("Missing min/max for : " + self.name)
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

                elif self.distrib == DistributionType.LOGNORMAL:

                    self._distrib = lognorm(self.default, self.std)

                elif self.distrib == DistributionType.BETA:
                    self._distrib = beta(
                        self.a,
                        self.b,
                        loc=self.default,
                        scale=self.std)

                else:
                    raise Exception("Unkown distribution type " + self.distrib)


            return self._distrib.ppf(alpha)

    def __hash__(self):
        return hash((self.dbname, self.name))

    def __eq__(self, other):
        if isinstance(other, ParamDef) :
            return self.name == other.name and self.dbname == other.dbname
        else :
            return Symbol.__eq__(self, other)

    # Expand parameter (useful for enum param)
    def expandParams(self, value=None) -> Dict[str, float]:
        if value == None:
            value = self.default
        return {self.name: value}

    # Useful for enum param, having several names
    def names(self, use_label=False):
        if use_label :
            return [self.get_label()]
        else:
            return [self.name]

    def __repr__(self):
        return self.name


class BooleanDef(ParamDef):
    """Parameter with discrete value 0 or 1"""

    def __init__(self, name, **argv):
        if not "min" in argv:
            argv = dict(argv, min=None, max=None)
        super(BooleanDef, self).__init__(name, ParamType.BOOL, **argv)


    def range(self, n):
        return [0, 1]

    def rand(self, alpha):
        return np.around(alpha)


class EnumParam(ParamDef):
    """Enum param is a facility wrapping a choice / switch as many boolean parameters.
    It is not itself a Sympy symbol. use #symbol("value") to access it.
    Statistics weight can be attached to values by providing a dict.
    """

    def __init__(self, name, values: Union[List[str], Dict[str, float]], **argv):

        if not "min" in argv :
            argv = dict(argv, min=None, max=None)
        super(EnumParam, self).__init__(name, ParamType.ENUM, **argv)
        if type(values) == list :
            self.values = values
            self.weights = {key:1 for key in values}
        else :
            self.weights = values
            self.values = list(values)
        self.sum = sum(self.weights.values())

    def expandParams(self, currValue=None):
        """
        Return a dictionarry of single enum values as sympy symbols, with only a single one set to 1.
        currValue can be either a single enum value (string) of dict of enum value => weight.
        """

        # A dict of weights was passed
        if isinstance(currValue, dict) :
            res = { "%s_%s" % (self.name, key) : val / self.sum for key, val in currValue.items()}
            res["%s_default" % self.name] = 0
            return res



        # Normal case
        values = self.values + [None]

        # Bad value ?
        if currValue not in values :
            raise Exception("Invalid value %s for param %s. Should be in %s" %
                            (currValue, self.name, str(self.values)))

        res = dict()
        for enum_val in values:
            var_name = "%s_%s" % (self.name, enum_val if enum_val is not None else "default")
            res[var_name] = 1.0 if enum_val == currValue else 0.0
        return res

    def symbol(self, enumValue):
        """Return the invididual symbol for a given enum value : <paramName>_<paramValue>"""
        if enumValue is None:
            return Symbol(self.name + '_default')
        if not enumValue in self.values:
            raise Exception("enumValue should be one of %s. Was %s" % (str(self.values), enumValue))
        return Symbol(self.name + '_' + enumValue)

    def names(self, use_label=False):
        if use_label :
            base_name = self.get_label()
        else :
            base_name = self.name
        return ["%s_%s" % (base_name, value) for value in (self.values + ["default"])]

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


def newParamDef(name, type, dbname=None, save=True, **kwargs):
    """
        Creates a parameter and register it into a global registry and as a brightway parameter.

        Parameters
        ----------

        type : Type of the parameter (From ParamType)
        save : Boolean, persist this into Brightway2 project (True by default)
        dbname : Optional name of db. If None, the parameter is a project parameter
        other arguments : Refer to the documentation of BooleanDef ParamDef and EnumParam

    """
    if type == ParamType.ENUM:
        param = EnumParam(name, dbname=dbname, **kwargs)
    elif type == ParamType.BOOL:
        param = BooleanDef(name, dbname=dbname, **kwargs)
    else:
        param = ParamDef(name, dbname=dbname, type=type, **kwargs)

    _param_registry()[name] = param

    # Save in brightway2 project
    if save :
        _persistParam(param)

    return param

_BOOLEAN_UNCERTAINTY_ATTRIBUTES = {
    UNCERTAINTY_TYPE: _UncertaintyType.DISCRETE,
    "minimum" : 0,
    "maximum" : 2 # upper bound + 1
}

def persistParams() :
    """ Persist parameters into Brightway project, as per :
     https://stats-arrays.readthedocs.io/en/latest/
    """

    for param in _param_registry().all() :
        _persistParam(param)

def _persistParam(param):
    """ Persist parameter into Brightway project """
    out = []

    # Common attributes for all types of params
    bwParam = dict(
        name=param.name,
        group=param.group,
        label=param.label,
        unit=param.unit,
        description=param.description)

    if (param.dbname) :
        bwParam["database"] = param.dbname

    if param.type == ParamType.ENUM :
        # Enum are not real params but a set of parameters
        for value in param.values :
            enumValueParam = dict(bwParam)
            enumValueParam["name"] = param.name + '_' + value
            enumValueParam.update(_BOOLEAN_UNCERTAINTY_ATTRIBUTES)
            # Use 'scale' as weight for this enum value
            enumValueParam['scale'] = param.weights[value]
            enumValueParam['amount'] = 1 if param.default == value else 0
            out.append(enumValueParam)
    else :

        bwParam["amount"] = param.default

        if param.type == ParamType.BOOL :
            # "Discrete uniform"
            bwParam.update(_BOOLEAN_UNCERTAINTY_ATTRIBUTES)

        elif param.type == ParamType.FLOAT :

            # Save uncertainty
            bwParam[UNCERTAINTY_TYPE] = _DistributionTypeMap[param.distrib]
            bwParam["minimum"] = param.min
            bwParam["maximum"] = param.max
            bwParam["loc"] = param.default

            if param.distrib in [DistributionType.NORMAL, DistributionType.LOGNORMAL] :
                bwParam["scale"] = param.std

            elif param.distrib == DistributionType.BETA:

                bwParam["scale"] = param.std
                bwParam["loc"] = param.a
                bwParam["shape"] = param.b

        else :
            error("Param type not supported", param.type)

        out.append(bwParam)

    if param.dbname :
        bw.parameters.new_database_parameters(out, param.dbname)
    else:
        bw.parameters.new_project_parameters(out)

def _loadArgs(data) :
    """Load persisted data attributes into ParamDef attributes"""
    return {
        "group": data.get("group"),
        "dbname" : data.get("database"),
        "default": data.get("amount"),
        "label": data.get("label"),
        "description": data.get("description"),
        "min": data.get("minimum"),
        "max": data.get("maximum"),
        "unit" : data.get("unit")
    }

def loadParams(global_variable=True, dbname=None):
    """
    Load parameters from Brightway database, as per : https://stats-arrays.readthedocs.io/en/latest/

    Parameters
    ----------
    global_variable If true, loaded parameters are made available as global variable.
    dbname : None. By default load all project and database parameters. If provided, only load DB params
    """

    enumParams=defaultdict(lambda : dict())

    def register(param) :
        _param_registry()[param.name] = param

        # Make it available as global var
        if global_variable:
            if param.name in builtins.__dict__ :
                error("Variable '%s' was already defined : overidding it with param." % param.name)
            builtins.__dict__[param.name] = param

    select = DatabaseParameter.select()
    if dbname :
        select = select.where(DatabaseParameter.database == dbname)

    params = list(select)
    if not dbname :
        params += list(ProjectParameter.select())

    for bwParam in params:
        data = bwParam.data
        data["amount"] = bwParam.amount
        name = bwParam.name

        type = data.get(UNCERTAINTY_TYPE, None)

        # print("Data for ", name, data)

        # Common extra args
        args = _loadArgs(data)

        if type == _UncertaintyType.DISCRETE :
            # Boolean or enum

            if data.get('scale') is not None :
                # Enum param : group them by common prefix
                splits = name.split("_")
                enum_value = splits.pop()
                enum_name = "_".join(splits)
                enumParams[enum_name][enum_value] = data
                continue

            elif data["maximum"] == 2 :
                del args["max"], args["min"]
                param = newBoolParam(name, save=False, **args)
            else:
                error("Non boolean discrete values (max != 2) are not supported for param :", name)
                continue
        else :
            # Float parameter
            if type is None:
                error("'Uncertainty type' of param %s not provided. Assuming UNIFORM")
                type = _UncertaintyType.UNIFORM

            # Uncertainty type to distribution type
            args["distrib"] = _DistributionTypeMapReverse[type]

            if type == _UncertaintyType.TRIANGLE :
                args["default"] = data["loc"]

            elif type in [_UncertaintyType.NORMAL, _UncertaintyType.LOGNORMAL]:
                args["default"] = data["loc"]
                args["std"] = data["scale"]

            elif type == _UncertaintyType.BETA:
                args["default"] = data["loc"]
                args["std"] = data["scale"]
                args["a"] = data["loc"]
                args["b"] = data["shape"]

            param = newFloatParam(name, save=False, **args)

        # Save it in shared dictionnary
        register(param)


    # Loop on EnumParams
    for param_name, param_values in enumParams.items() :
        first_enum_param = list(param_values.values())[0]
        args = _loadArgs(first_enum_param)
        del args["default"]

        # Dictionary of enum values with scale as weight
        args["values"] = {key : data["scale"] for key, data in param_values.items()}

        # Default enum value is the one with amount=1
        defaults = list(key for key, data in param_values.items() if data.get("amount") == 1)
        if len(defaults) == 1 :
            default = defaults[0]
        else :
            default= None
            error("No default enum value found for ", param_name, defaults)

        param = newEnumParam(param_name, default, save=False, **args)

        # Save it in shared dict
        register(param)


def newFloatParam(name, default, dbname=None, **kwargs):
    """ Create a FLOAT parameter. See the documentation of arguments for #newParamDef()."""
    return newParamDef(name, ParamType.FLOAT, dbname=dbname, default=default, **kwargs)

def newBoolParam(name, default, dbname=None, **kwargs):
    """ Create a BOOL parameter. See the documentation of arguments for #newParamDef()."""
    return newParamDef(name, ParamType.BOOL, dbname=dbname, default=default, **kwargs)

def newEnumParam(name, default, dbname=None, **kwargs):
    """ Create a ENUM parameter. See the documentation of arguments for #newParamDef()."""
    return newParamDef(name, ParamType.ENUM, dbname=dbname, default=default, **kwargs)

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

class DuplicateParamsAndNoContextException(Exception) :
    pass

class ParamRegistry :
    """ In memory registry of parameters, acting like a dict and maintaining parameters with possibly same names on several DBs"""
    def __init__(self) :
        # We store a dict of dict
        # Param Name -> { dbname ->  param}
        self.params : Dict[Dict[ParamDef]] = defaultdict(dict)

    def __len__(self):
        return len(self.params)

    def __getitem__(self, key):
        params_per_db = self.params[key]
        if len(params_per_db) == 1 :
            return list(params_per_db.values())[0]

        if DbContext.current_db() == None :
            dbs = [key or '<project>' for key in params_per_db.keys()]
            raise DuplicateParamsAndNoContextException(

                """
                Found several params with name '%s', linked to databases (%s) . Yet no context is provided. 
                Please embed you code in a DbContext :
                    with DbContext(currentdb) :
                        <code>
                """ % (key, ", ".join(dbs)))
        if DbContext.current_db() in params_per_db :
            return params_per_db[DbContext.current_db()]
        else :
            return params_per_db[None]

    def __setitem__(self, key, param : ParamDef):
        if param.dbname in self.params[key]:
            error("[ParamRegistry] Param %s was already defined in '%s' : overriding." %
                  (param.name, param.dbname or '<project>'))

        self.params[key][param.dbname] = param

    def __contains__(self, key):
        return key in self.params

    def values(self) :
        return [self.__getitem__(key) for key in (self.params)]

    def keys(self) :
        return self.params.keys()

    def items(self) :
        return [(key, self.__getitem__(key)) for key in self.params.keys()]

    def clear(self, db_name=None):
        if db_name is None :
            self.params.clear()
        else :
            for param_name, db_params in self.params.items() :
                if db_name in db_params :
                    del db_params[db_name]

    def all(self) :
        """Return list of all parameters, including params with same names and different DB"""
        return list(param for params in self.params.values() for param in params.values())

# Possible param values : either floator string (enum value)
ParamValue = Union[float, str]

# Single value or list of values
ParamValues = Union[List[ParamValue],  ParamValue]


def _param_registry() -> Dict[str, ParamDef] :
    # Prevent reset upon auto reload in jupyter notebook
    if not 'param_registry' in builtins.__dict__:
        builtins.param_registry = ParamRegistry()

    return builtins.param_registry

def _completeParamValues(params: Dict[str, ParamValues], required_params : List[str]=None, setDefaults=False) :
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
                error("Required param '%s' was missing, replacing by default value : %s" % (param_name, str(param.default)))

    # Set default variables for missing values
    if setDefaults :
        for name, param in _param_registry().items() :
            if not name in params :
                params[name] = param.default

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


def resetParams(db_name=None):
    """Clear parameters in live memory (registry) and on disk.
    Clear either all params (project and all db params) or db params from a single database (if db_name provided)"""
    _param_registry().clear(db_name)

    if db_name is None :
        ProjectParameter.delete().execute()
        ActivityParameter.delete().execute()
        DatabaseParameter.delete().execute()
    else :
        ActivityParameter.delete().where(ActivityParameter.database == db_name).execute()
        DatabaseParameter.delete().where(DatabaseParameter.database == db_name).execute()
        Group.delete().execute()

class NameType(Enum) :
    NAME="name"
    LABEL="label"
    CAMEL_NAME = "CAMEL_NAME"

def _param_name(param, name_type:NameType) :
    if name_type == NameType.NAME :
        return param.name
    elif name_type == NameType.LABEL :
        return param.get_label()
    else :
        return _snake2camel(param.name)

def list_parameters(name_type=NameType.NAME):

    """ Print a pretty list of all defined parameters """
    params = [dict(
        group=param.group or "",
        name=_param_name(param, name_type),
        label=param.get_label(),
        default=param.default,
        min=param.min,
        max=param.max,
        std=getattr(param, "std", None),
        distrib=param.distrib,
        unit=param.unit,
        db=param.dbname or "[project]") for  param in _param_registry().all()]

    groups = list({p["group"] for p in params})
    groups = sorted(groups)

    # Sort by Group / name
    def keyf(param) :
        return (groups.index(param["group"]), param["name"])

    sorted_params = sorted(params, key=keyf)

    return HTML(tabulate(sorted_params, tablefmt="html", headers="keys"))


def freezeParams(db_name, **params) :
    """
    Freeze parameters values in all exchanges amounts of a DB.
    The formulas are computed and the 'amount' attributes are set with the result.
    This enables parametric datasets to be used by standard, non parametric tools of Brightway2.
    """

    db = bw.Database(db_name)

    with DbContext(db) :
        for act in db :
            for exc in act.exchanges():

                amount = _getAmountOrFormula(exc)

                # Amount is a formula ?
                if isinstance(amount, Basic):

                    replace = [(name, value) for name, value in _completeParamValues(params, setDefaults=True).items()]
                    val = amount.subs(replace).evalf()

                    with ExceptionContext(val) :
                        val = float(val)

                    print("Freezing %s // %s : %s => %d" % (act, exc['name'], amount, val))

                    # Update in DB
                    exc["amount"] = val
                    exc.save()


def _listParams(db_name) -> List[ParamDef]:
    """
    Return a set of all parameters used in activities
    """

    db = bw.Database(db_name)
    res = set()

    with DbContext(db) :
        for act in db :
            for exc in act.exchanges():

                amount = _getAmountOrFormula(exc)

                # Amount is a formula ?
                if isinstance(amount, Basic):
                    expanded_names = list(str(symbol) for symbol in amount.free_symbols)
                    param_names = _expanded_names_to_names(expanded_names)
                    params = list(_param_registry()[param_name] for param_name in param_names)
                    res.update(params)
    return res


def _expand_param_names(param_names):
    '''Expand parameters names (with enum params) '''
    return [name for key in param_names for name in _param_registry()[key].names()]


def _expanded_names_to_names(param_names):
    """Find params corresponding to expanded names, including enums."""
    param_names = set(param_names)

    # param name => param
    res = dict()

    # Search for param with same name of prefix paramName_enumValue
    for expended_name in param_names :
        for param_name in _param_registry().keys():
            if expended_name.startswith(param_name) :
                param = _param_registry()[param_name]
                for name in param.names():
                    if name == expended_name:
                        res[expended_name] = param

    missing = param_names - set(res.keys())
    if len(missing) > 0:
        raise Exception("Unkown params : %s" % missing)

    return {param.name for param in res.values()}
