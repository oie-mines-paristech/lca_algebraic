from contextlib import AbstractContextManager
from functools import cache
from inspect import isfunction
from typing import Dict, Iterable, Tuple, Union

import sympy
import ipywidgets as widgets
import numpy as np
import pandas as pd
from IPython.display import display
from six import raise_from
from sympy import Basic
from pint import Quantity

from lca_algebraic.settings import Settings
from lca_algebraic.bw_wrapper import Activity, Database, is_output_exchange

_user_functions = dict()


def _isOutputExch(exc):
    return is_output_exchange(exc)


def _isnumber(value):
    return isinstance(value, int) or isinstance(value, float)


dbs = dict()


def invalidate_db(dbname=None):
    """Drop the pooled Database instances so a fresh one is built on next access.

    bw25 keeps an internal (session/ORM) view inside a Database object that is not
    invalidated when a database is deleted or recreated. Keep the pool in sync so that
    activities are not resolved to stale (deleted) node ids.
    """
    global dbs
    if dbname is None:
        dbs = dict()
    else:
        dbs.pop(dbname, None)


def _getDb(dbname):
    """Pool of Database instances"""
    if dbname not in dbs:
        dbs[dbname] = Database(dbname)
    return dbs[dbname]


def _MinMax(op, arg0, *args):
    from lca_algebraic.units import unit_registry as u
    if not Settings.units_enabled:
        return op(arg0, *args)
    same_unit = all(x.units == arg0.units for x in args)
    if not same_unit:
        raise Exception("MinMax argments must have the same units")
    return u.Quantity(op(arg0.magnitude, *[x.magnitude for x in args]), arg0.units)


def Max(arg0, *args):
    return _MinMax(sympy.Max, arg0, *args)


def Min(arg0, *args):
    return _MinMax(sympy.Min, arg0, *args)


def _actDesc(act: Activity):
    """Generate pretty name for activity + basic information"""
    name = _actName(act)
    amount = act.getOutputAmount()

    return "%s (%f %s)" % (name, amount, act["unit"])


def _actName(act: Activity):
    """Generate pretty name for activity, appending location if not 'GLO'"""
    res = act["name"]
    if "location" in act and act["location"] != "GLO":
        res += "[%s]" % act["location"]
    return res


def displayWithExportButton(df):
    """Display dataframe with option to export"""

    button = widgets.Button(description="Export data")
    button.style.button_color = "lightgray"

    def click(e):
        df.to_csv("out.csv")
        button.description = "exported as 'out.csv'"

    dfout = widgets.Output()
    with dfout:
        display(df)

    button.on_click(click)

    display(widgets.VBox([button, dfout]))


def as_np_array(a):
    if isinstance(a, list):
        return np.asarray(a)
    else:
        return a


def r_squared(y, y_hat):
    y_bar = y.mean()
    ss_tot = ((y - y_bar) ** 2).sum()
    ss_res = ((y - y_hat) ** 2).sum()
    return 1 - (ss_res / ss_tot)


class ExceptionContext(AbstractContextManager):
    def __init__(self, context):
        self.context = context

    def __exit__(self, exc_type, exc_val, exc_tb):
        if exc_val is not None:
            raise_from(Exception("Context : %s" % str(self.context)), exc_val)
        return True


def _snake2camel(val):
    return "".join(word.title() for word in val.split("_"))


class TabbedDataframe:
    """This class holds a dictionnary of dataframes and can display and saved them awith 'tabs'/'sheets'"""

    def __init__(self, metadata=dict(), **dataframes):
        self.dataframes = dataframes
        self.metadata = metadata

    def __str__(self):
        res = ""
        for name, df in self.dataframes.items():
            res += f"\n{name} : \n"
            res += df.__str__() + "\n"
        return res

    def _repr_html_(self):
        display(_mk_tabs(self.dataframes))

    def to_excel(self, filename):
        assert filename.endswith(".xlsx")

        with pd.ExcelWriter(filename, engine="xlsxwriter") as writer:
            for itab, (name, df) in enumerate(self.dataframes.items()):
                if itab == 0:
                    df.to_excel(writer, sheet_name=name, startrow=len(self.metadata) + 1)

                    # Write metadata in header
                    worksheet = writer.sheets[name]
                    for imeta, (key, val) in enumerate(self.metadata.items()):
                        worksheet.write_string(imeta, 0, str(key))
                        worksheet.write_string(imeta, 1, str(val))

                else:
                    df.to_excel(writer, sheet_name=name)


def _mk_tabs(titlesAndContent: Dict):
    """Generate iPywidget tabs"""
    tabs = []
    titles = []
    for title, content in titlesAndContent.items():
        titles.append(title)

        tab = widgets.Output()
        with tab:
            if isfunction(content):
                content()
            else:
                display(content)
        tabs.append(tab)

    res = widgets.Tab(children=tabs)
    for i, title in enumerate(titles):
        res.set_title(i, title)
    return res


def _display_tabs(titlesAndContent: Dict):
    display(_mk_tabs(titlesAndContent))


def one(it: Iterable):
    """Expect a list with single value a returns it"""
    it = list(it)
    if len(it) != 1:
        raise Exception(f"Expected a single value but got {len(it)}")
    return it[0]


@cache
def getActByCode(db_name, code):
    """Get activity by code"""
    return _getDb(db_name).get(code)


# Types
ValueOrExpression = Union[float, Basic, Quantity]
MethodKey = Tuple[str, ...]
