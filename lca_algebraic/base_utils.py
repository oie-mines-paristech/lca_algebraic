from contextlib import AbstractContextManager
from sys import stderr

import brightway2 as bw
from bw2data.backends.peewee import Activity
from six import raise_from
import ipywidgets as widgets
from IPython.core.display import display
import numpy as np

DEBUG=False
LANG="fr"
UNIT_OVERRIDE = dict()


def set_debug(value=True) :
    """ Activate debug logs """
    global DEBUG
    DEBUG=value

def set_lang(lang):
    """Set language """
    global LANG
    LANG=lang

def debug(*args, **kwargs) :
    if DEBUG :
        print(*args, **kwargs)

def error(*args, **kwargs):
    """Print message on stderr """
    print(*args, **kwargs, file=stderr)


def _isOutputExch(exc) :
    return exc.get("type") == "production"


def _isnumber(value):
    return isinstance(value, int) or isinstance(value, float)

dbs = dict()
def _getDb(dbname) -> bw.Database:
    """Pool of Database instances"""
    if not dbname in dbs:
        dbs[dbname] = bw.Database(dbname)
    return dbs[dbname]


def interpolate(x, x1, x2, y1, y2):
    """Build an expression for linear interpolation between two points.
    If x is not within [x1, x2] the corresponding bound Y values are returned"""
    x = Min(Max(x, x1), x2)
    return y1 + (y2 - y1) * (x - x1) / (x2 - x1)


def Max(a, b) :
    """Max define as algrebraic forumal with 'abs' for proper computation on vectors """
    return (a + b + abs(a - b)) / 2


def Min(a, b) :
    """Max define as algrebraic forumal with 'abs' for proper computation on vectors """
    return (a + b - abs(b - a)) / 2


def _actDesc(act: Activity):
    """Generate pretty name for activity + basic information """
    name = _actName(act)
    amount = act.getOutputAmount()

    return "%s (%f %s)" % (name, amount, act['unit'])


def _method_unit(method) :
    if method in UNIT_OVERRIDE :
        return UNIT_OVERRIDE[method]
    return bw.Method(method).metadata['unit']


def _actName(act: Activity):
    
    """Generate pretty name for activity, appending location if not 'GLO' """
    res = act['name']
    if 'location' in act and act['location'] != 'GLO':
        res += "[%s]" % act["location"]
    return res


def displayWithExportButton(df):
    '''Display dataframe with option to export'''

    button = widgets.Button(description="Export data")
    button.style.button_color = "lightgray"
    def click(e) :
        df.to_csv("out.csv")
        button.description = "exported as 'out.csv'"
    dfout = widgets.Output()
    with dfout :
        display(df)

    button.on_click(click)

    display(widgets.VBox([button, dfout]))


def as_np_array(a) :
    if type(a) == list :
        return np.asarray(a)
    else :
        return a

def r_squared(y, y_hat):
    y_bar = y.mean()
    ss_tot = ((y - y_bar) ** 2).sum()
    ss_res = ((y - y_hat) ** 2).sum()
    return 1 - (ss_res / ss_tot)


class  ExceptionContext(AbstractContextManager) :
    def __init__(self, context):
        self.context = context

    def __exit__(self, exc_type, exc_val, exc_tb):
        if exc_val != None :
            raise_from(Exception("Context : %s" % str(self.context)), exc_val)
        return True

def _snake2camel(val):
    return ''.join(word.title() for word in val.split('_'))


class SymDict :
    """ This class acts like a dict with arithmetic operations. It is useful to process 'axes' LCA computations """

    def __init__(self, values) :
        self.dict = values

    def _apply_op(self, other, fop, null_val):

        # None is the key for non flagged values
        if not isinstance(other, SymDict) :
            dic = dict()
            dic[None] = other
            other = SymDict(dic)

        all_keys = set(other.dict.keys()) | set(self.dict.keys())
        return SymDict({
            key: fop(
                self.dict.get(key, null_val),
                other.dict.get(key, null_val)) for key in all_keys
        })

    def _apply_self(self, fop):
        return SymDict({key: fop(val) for key, val in self.dict.items()})

    def __add__(self, other):
        return self._apply_op(other, lambda a, b: a + b, 0)

    def __radd__(self, other):
        return self._apply_op(other, lambda a, b: b + a, 0)

    def __mul__(self, other):
        return self._apply_self(lambda a : a*other)

    def __rmul__(self, other):
        return self._apply_self(lambda a : other*a)

    def __truediv__(self, other):
        return self._apply_self(lambda a: a / other)

    def __rtruediv__(self, other):
        return NotImplemented

    def __repr__(self):
        return "{" + \
            "; ".join("%s: %s" % (key, str(val)) for key, val in self.dict.items()) + \
        "}"

    def _defer(self, funcname, args, kwargs):

        return SymDict({
            key: val if not hasattr(val, funcname) else getattr(val, funcname)(*args, **kwargs)
            for key, val in self.dict.items()
        })

    def xreplace(self, *args, **kwargs):
        return self._defer("xreplace", args, kwargs)

    @property
    def free_symbols(self):
        res = set()
        for key, val in self.dict.items() :
            if hasattr(val, "free_symbols") :
                res |= val.free_symbols
        return list(res)