from contextlib import AbstractContextManager
from sys import stderr
from typing import Union

import brightway2 as bw
from bw2data.backends.peewee import Activity, ExchangeDataset
from six import raise_from
from sympy import Basic
from sympy.parsing.sympy_parser import parse_expr
import ipywidgets as widgets
from IPython.core.display import display
import numpy as np
import os

DEBUG=False
LANG="fr"
UNIT_OVERRIDE = dict()

_user_functions = dict()

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
    return exc.get('input') == exc.get('output') or exc.get("type") == "production"


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
    amount = 1
    for ex in act.exchanges() :
        if _isOutputExch(ex):
            amount = ex['amount']

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


def _getAmountOrFormula(ex: ExchangeDataset) -> Union[Basic, float]:
    """ Return either a fixed float value or an expression for the amount of this exchange"""
    if 'formula' in ex:
        try:
            local_dict = {x[0].name: x[0] for x in _user_functions.values()}
            return parse_expr(ex['formula'], local_dict=local_dict)
        except:
            error("Error while parsing formula '%s' : backing to amount" % ex['formula'])

    return ex['amount']


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


def eprint(msg) :
    """Print directly to os, bypassing jupyter"""
    os.write(2, (msg+"\n").encode())



def with_output(out):
    """Decorator for wrapping any output. Useful for callback function to print errors """

    def decorator(func):
        def wrapper(*args, **kwargs):
            with out:
                return func(*args, **kwargs)
        return wrapper
    return decorator