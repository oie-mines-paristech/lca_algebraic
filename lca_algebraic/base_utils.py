from contextlib import AbstractContextManager
from inspect import isfunction
from typing import Dict, Iterable, Tuple, Union

import brightway2 as bw
import ipywidgets as widgets
import numpy as np
import pandas as pd
from bw2data.backends.peewee import Activity
from IPython.display import display
from six import raise_from
from sympy import Basic
from sympy.physics.units import Quantity

_user_functions = dict()


def _isOutputExch(exc):
    return exc.get("type") == "production"


def _isnumber(value):
    return isinstance(value, int) or isinstance(value, float)


dbs = dict()


def _getDb(dbname) -> bw.Database:
    """Pool of Database instances"""
    if dbname not in dbs:
        dbs[dbname] = bw.Database(dbname)
    return dbs[dbname]


def Max(a, b):
    """Max define as algrebraic forumal with 'abs' for proper computation on vectors"""
    return (a + b + abs(a - b)) / 2


def Min(a, b):
    """Max define as algrebraic forumal with 'abs' for proper computation on vectors"""
    return (a + b - abs(b - a)) / 2


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


def getActByCode(db_name, code):
    """Get activity by code"""
    return _getDb(db_name).get(code)


# Types
ValueOrExpression = Union[float, Basic, Quantity]
MethodKey = Tuple[str]
