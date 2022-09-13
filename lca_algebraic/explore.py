from collections import defaultdict
from enum import Enum
from functools import partial
from types import SimpleNamespace

import brightway2 as bw
import ipysheet
import ipywidgets as widgets
import xlsxwriter
from IPython.display import display


from ipyevents import Event
from ipywidgets import VBox, Button, Output, Combobox, Layout, ToggleButton, HBox, HTML
from sympy import Basic
from werkzeug.utils import secure_filename

from lca_algebraic import ActivityExtended, multiLCAAlgebric, getActByCode, method_name
from lca_algebraic.base_utils import _actName, eprint
from lca_algebraic.lca import _createTechProxyForBio
from lca_algebraic.params import _completeParamValues
from lca_algebraic.stats import _round_expr
from ipytree import Tree, Node

def _uniq_name(name, names) :
    i = 1
    res = name
    while res in names:
        res = "%s#%d" % (name, i)
        i += 1
    return res


def _ellipsis(value, max=40) :
    if len(value) < max :
        return value
    else :
        short = value[0:max]
        return "<span title='%s'>%s…</span>" % (value, short)


def _sci_repr(float_val) :
    if isinstance(float_val, float):
        return "%02.3g" % float_val
    else :
        return str(float_val)


def _none_float(float_val) :
    return 0.0 if float_val is None else float_val


_SMALL_DIFF_COLOR = "#FFFACD" # LemonChiffon
_HIGH_DIFF_COLOR = "#F08080" # Light coral
_RIGHT_BORDER = {'border-right-width':'2px', 'border-right-color':'black'}
_DB_COLORS = ["lightgreen", "lightyellow", "lightsteelblue", "lightpink", "NavajoWhite", "khaki", "LightCoral"]


def _db_shortname(db_name) :
    """Return 'shortname' metadata if Db has one. Otherwize, return DB name"""
    db = bw.Database(db_name)
    if "shortname" in db.metadata :
        return db.metadata["shortname"]
    return db_name


def _db_idx(db) :
    """Return index if DB, in list of all DBs"""
    dbs = list(bw.databases.keys())
    return dbs.index(db)


class ActNode :
    def __init__(self, data):
        self.data = data
        self.parent = None
        self.children = []

class keydefaultdict(defaultdict):
    def __missing__(self, key):
        if self.default_factory is None:
            raise KeyError(key)
        else:
            ret = self[key] = self.default_factory(key)
            return ret


def act_tree(db_name) :
    """ Build trees of activities for a given DB """

    nodes_dict = keydefaultdict(lambda act: ActNode(act))

    for act in bw.Database(db_name) :
        if "isProxy" in act :
            continue

        parent_node = nodes_dict[act]
        for exch in act.listExchanges() :
            input = exch.input

            if input.key[0] != db_name or "isProxy" in act :
                continue

            if input is nodes_dict :
                continue

            child_node = nodes_dict[input]
            child_node.parent = parent_node

    for node in nodes_dict.values() :
        if node.parent is not None :
            node.parent.children.append(node)

    # Return list of all root nodes
    return [node for node in nodes_dict.values() if node.parent is None]

def print_tree(db_name) :

    def print_rec(node, indent=0) :
        print(indent * " " + _actName(node.data))
        for child in node.children:
            print_rec(child, indent+2)

    roots = act_tree(db_name)
    for root in roots:
        print_rec(root)

def show_tree(db_name, callback=None) :
    res = Tree()
    roots = act_tree(db_name)

    def transform(act_node) :
        children = list(transform(child) for child in act_node.children)
        node = Node(
            _actName(act_node.data),
            children,
            opened=act_node.parent is not None)

        node.data = act_node.data.key[1]

        if callback:
            def handle_click(event) :
                id = event['owner'].data
                act = getActByCode(db_name, id)
                callback(act)

            node.observe(handle_click, 'selected')

        return node

    nodes = list(transform(root) for root in roots)

    for node in sorted(nodes, key=lambda node : node.name):
        res.add_node(node)

    return res



def compare_impacts(
        activities : ActivityExtended,
        method,
        default_params=True,
        return_data=False,
        **params):
    """
    Advanced version of #printAct()

    Displays all exchanges of one or several activities and their impacts.
    If parameter values are provided, formulas will be evaluated accordingly.
    If two activities are provided, they will be shown side by side and compared.

    :param return_data: If true return attributes (including impacts). If false, formats and returns a live sheet
    :param activities: One activity, or couple of activities
    :param method: Impact method
    :param default_params: If true, replace params with default values in formulas
    :param params: Extra parameter values to replace in formulas
    """

    if not isinstance(activities, (list, tuple)):
        activities = [activities]

    # List of dict[exchange=>attributes]
    # As many items as source activities
    datas = []

    for main_act in activities:

        inputs_by_ex_name = dict()

        # Dict of exchange name => attributes
        data = dict()
        datas.append(data)

        for ex in main_act.listExchanges():

            amount = ex.amount
            formula = None

            # Evaluate formulas
            if isinstance(ex.amount, Basic) :

                # Keep formula
                formula = _round_expr(amount, 3)

                # Replace params in to compute a float amount
                new_params = [(name, value) for name, value in _completeParamValues(params, setDefaults=default_params).items()]
                amount = ex.amount.subs(new_params)

            ex_name = _uniq_name(ex.name, data)
            inputs_by_ex_name[ex_name] = _createTechProxyForBio(ex.input.key, main_act.key[0])

            # Save data for this exchange
            data[ex_name] = SimpleNamespace(
                input=ex.input["name"],
                formula=formula,
                amount=amount,
                db= ex.input.key[0],
                code= ex.input.key[1],
                loc="GLO" if not "location" in ex.input else ex.input["location"],
                unit=ex.unit)

        # Provide impact calculation if impact provided
        all_acts = list(set(inputs_by_ex_name.values()))

        impacts = multiLCAAlgebric(all_acts, [method], raw=True, **params)

        # Add impacts to data
        for key, attrs in data.items() :
            amount = attrs.amount
            act = inputs_by_ex_name[key]
            impact = impacts[(act, method)]

            attrs.impact = amount * impact

    if return_data :
        return datas

    sheet = output_sheet(activities, datas)
    return sheet


class SheetFormat(str, Enum) :
    EXCEL = "xls"
    NOTEBOOK = "notebook"


_LIGHT_GREY = "#808080"
_LIGHT_GREEN = "#AFE1AF"

def output_sheet(
        activities,
        datas,
        diff_impact=0.1,
        act_callback=None,
        out_format: SheetFormat = SheetFormat.NOTEBOOK,
        xls_filename= "out.xlsx"):

    """
    :param activities: List of activities
    :param datas: Impact data as output by compare_impacts
    :param diff_impact: Min minimum impact difference to highlight, as part of total impact
    :param act_callback: Callback when clicking on activity link
    :param out_format: Excel or live / notebook
    """

    if not isinstance(activities, (list, tuple)):
        activities = [activities]

    nact = len(activities)

    # Compute sum of impacts
    sum_impacts = [0.0, 0.0]
    for idata, data in enumerate(datas) :
        for attrs in data.values() :
            sum_impacts[idata] += attrs.impact

    # List of all exchange names
    all_ex_names = list(set(ex for data in datas for ex in data.keys()))
    headers = ["exchange", "unit"]
    # Exchange, Unit, Inputs, Amounts, Impacts
    widths = [170, 50] \
             + [210] * nact \
             + [70] * nact \
             + [70] * nact

    def _numerate(str_val):
        if nact == 1:
            return [str_val]
        else:
            return ["%s#%s" % (str_val, _db_shortname(act.key[0])) for act in activities]

    headers += _numerate("activity") + _numerate("amount") + _numerate("impact")
    sheet = None

    def draw_cell(irow, column, value, background_color=None, style=None, bold=False):
        nonlocal sheet

        if out_format == SheetFormat.NOTEBOOK:
            ipysheet.cell(
                row=irow,
                column=column,
                style=style.copy() if style else None,
                background_color=background_color,
                value=value,
                read_only=True)
        else:

            format = workbook.add_format()
            if background_color is not None:
                format.set_bg_color(background_color)
            if bold:
                format.set_bold()
            sheet.write(irow + 1, column, value, format)

    if out_format == SheetFormat.NOTEBOOK :

        # To sheet
        sheet = ipysheet.sheet(
            rows=len(all_ex_names)+1,
            columns=len(headers),
            column_headers=headers,
            row_headers=False,
            stretch_headers="none",
            column_width=widths)
    else :
        workbook = xlsxwriter.Workbook(xls_filename)
        sheet = workbook.add_worksheet()

        # Write header
        [draw_cell(-1, i, header, bold=True, background_color=_LIGHT_GREY) for i, header in enumerate(headers)]


    def draw_name(irow, name, color="#FFFFFF"):

        if out_format == SheetFormat.NOTEBOOK :
            value = HTML("<b class='trunc' title='%s'>%s</b>" % (name, name))
        else:
            value = name

        draw_cell(
            irow=irow, column=0,
            background_color=color,
            bold=True,
            value=value)

    def draw_unit(irow, unit, color="#FFFFFF"):
        draw_cell(irow=irow, column=1,
                      background_color=color,
                      value=unit, style=_RIGHT_BORDER)

    def draw_col(irow, col_name, offset, cell_value, color="#FFFFFF"):

        icol = 2 + nact * (0 if col_name == "input" else 1 if col_name == "amount" else 2)

        style = _RIGHT_BORDER if (offset == nact - 1) else None

        if isinstance(cell_value, float) and out_format == SheetFormat.NOTEBOOK :
            cell_value == _sci_repr(cell_value)

        draw_cell(
            irow=irow,
            column=icol + offset,
            style=style.copy() if style else None,
            background_color=color,
            value=cell_value)

    # First line : draw Total / output

    draw_name(0, "Total output", _LIGHT_GREEN)
    draw_unit(0, activities[0]["unit"], _LIGHT_GREEN)
    for iact, (act, sum_impact) in enumerate(zip(activities, sum_impacts)):
        draw_col(0, "impact", iact, sum_impact, _LIGHT_GREEN)
        draw_col(0, "amount", iact, act.getOutputAmount(), _LIGHT_GREEN)
        draw_col(0, "input", iact, _actName(act), _LIGHT_GREEN)

    # Loop on exchanges
    for i_ex, ex_name in enumerate(all_ex_names):

        row = i_ex + 1

        draw_name(row, ex_name)

        # Loop on type of columns
        for column in ["input", "amount", "impact"]:

            # On or two cells, one per activity
            values = []
            cells = []

            # Loop on activity (1 or 2)
            for i_act, (act, data) in enumerate(zip(activities, datas)):

                # Empty attr for this exchange ?
                attrs = data[ex_name] if ex_name in data else SimpleNamespace(
                    unit="-",
                    input="_",
                    amount="_",
                    formula=None,
                    impact=None,
                    db="-",
                    loc="-")

                # Override unit (common column)
                if attrs.unit != "-":
                    draw_unit(row, attrs.unit)

                # Switch on column type
                if column == "input":  # Input activity

                    if attrs.input == "_":
                        input = "-"
                        values.append(input)
                    else:

                        if out_format == SheetFormat.NOTEBOOK :

                            input_name = "<a href='#'>%s</a>" % attrs.input
                            loc_html = "<span class='loc'> %s </span>" % attrs.loc
                            db_html = "&nbsp;<span class='db%d'> %s </span>" % (_db_idx(attrs.db), _db_shortname(attrs.db))
                            html = input_name + "<br/>" + loc_html + db_html

                            values.append(html)

                            input = HTML(html)

                        else : # to Excel : text only
                            sub_act = getActByCode(attrs.db, attrs.code)
                            input = _actName(sub_act) + " #" + _db_shortname(attrs.db)
                            values.append(input)
                            cells.append(input)


                        # Link the callback to click on it
                        if act_callback is not None:
                            sub_act = getActByCode(attrs.db, attrs.code)

                            def cb(sub_act, event):
                                act_callback(sub_act)

                            event = Event(source=input, watched_events=['click'])
                            event.on_dom_event(partial(cb, sub_act))

                    cells.append(input)

                elif column == "amount":  # Amount and formula

                    values.append(attrs.amount)

                    amount_repr = attrs.amount
                    if attrs.formula and out_format == SheetFormat.NOTEBOOK :
                        amount_repr = HTML("<abbr title='%s'><b>%s</b></abbr>" % (str(attrs.formula), _sci_repr(attrs.amount)))

                    cells.append(amount_repr)

                else: # Impact

                    values.append(attrs.impact)
                    cells.append(attrs.impact)

            # Colorize background according to Diff
            color = None

            if (nact > 1) and values[0] != values[1]:

                # For impact, highlight differently difference higher than a given share of total impact
                if column == "impact":
                    diff = abs(
                        _none_float(values[1]) -
                        _none_float(values[0]))
                    rel_diff = diff / sum_impacts[0]

                    color = _HIGH_DIFF_COLOR if rel_diff > diff_impact else _SMALL_DIFF_COLOR

                else:
                    color = _SMALL_DIFF_COLOR

            # Display cells for this column
            for icell, cell in enumerate(cells):
                draw_col(row, column, icell, cell, color)

    if out_format == SheetFormat.NOTEBOOK :
        # Styling
        style = """
        .trunc {
               white-space: nowrap;
               overflow: hidden;
               text-overflow: ellipsis;
        }"""
        # Add colors for each DB number
        style += "".join(""".db%d {
                background-color : %s    
            }""" % (i, color) for i, color in enumerate(_DB_COLORS))
        style = HTML("""<style>%s</style>""" % style)
        display(style)
        return sheet

    else :
        workbook.close()



def _acts_by_name(db_name) :
    # Filter out proxy /dummy activities
    return {_actName(act):act for act in bw.Database(db_name) if not "isProxy" in act}


def explore_impacts(
        main_db,
        base_db,
        method):
    """
    Interactive interface to explore activities, compare them to other activities of other DB and compare impacts

    :param main_db: Name of main DB
    :param base_db: Name of base DB : to search reference activities to compare to
    :param method: Method for impact calculation
    :return: Interactive widget ready to be displayed
    """
    message = HTML()
    table_container = VBox()
    tree_container = VBox()
    history = []

    home_button = Button(icon="home")
    tree_button = Button(icon="folder")
    back_button = Button(icon="arrow-left")
    save_button = Button(icon="save")

    output = Output()

    # Current state
    main_act = None
    base_act = None
    current_data = None

    # List activities by name in other DB
    main_acts_by_name = _acts_by_name(main_db)
    base_acts_by_name = _acts_by_name(base_db)


    main_combo = Combobox(
        description=_db_shortname(main_db),
        layout=widgets.Layout(width='700px'),
        options=list(main_acts_by_name.keys()),
        placeholder="Select activity",
        ensure_option=True)

    base_combo = Combobox(
        description=_db_shortname(base_db),
        layout=widgets.Layout(width='700px'),
        options=list(base_acts_by_name.keys()),
        placeholder="Select activity",
        ensure_options=True)

    # Comments
    main_comment = HTML(layout=Layout(display="none"))
    base_comment = HTML(layout=Layout(display="none"))

    main_button_switch = ToggleButton(
        value=False,
        icon="chevron-down",
        description="Comment")

    base_button_switch = ToggleButton(
        value=False,
        icon="chevron-down",
        description="Comment")

    with output :

        def update_main(act) :
            """Update main activity. Try to update base activity with same name"""

            nonlocal main_act, base_act

            if main_act == act :
                return

            if main_act is not None or base_act is not None :
                history.append((main_act, base_act))

            main_act = act


            name = _actName(act)

            # Try to find reference activity in base
            base_act = None
            if act.key[0] == main_db :

                # The activity already references where it comes from
                if "inherited_from" in main_act :
                    base_act = getActByCode(*main_act["inherited_from"])

                # Otherwize, try to find another activity with same name
                elif name in base_acts_by_name:
                    base_act = base_acts_by_name[name]

            update_view()

        def update_base(act) :

            nonlocal base_act

            if base_act == act:
                return

            base_act = act

            history.append((main_act, base_act))
            update_view()

        def pop_and_update(event) :
            nonlocal main_act, base_act
            if len(history) > 0 :
                main_act, base_act = history.pop()
                update_view()

        def reset(event) :
            nonlocal main_act, base_act
            if len(history) > 0 :
                main_act, base_act = history[0]
                history.clear()
                update_view()

        def update_view():
            nonlocal current_data

            # Update names in combo
            main_combo.value = "" if main_act is None else _actName(main_act)
            base_combo.value = "" if base_act is None else _actName(base_act)

            # Update comments
            main_comment.value = "" if main_act is None else (main_act["comment"] if "comment" in main_act else "")
            base_comment.value = "" if base_act is None else (base_act["comment"] if "comment" in base_act else "")

            message.value = "<h4>Computing ...</h4>"

            try:
                acts = main_act if base_act is None else [main_act, base_act]

                # Prepare data for both notebook and excel output
                current_data = compare_impacts(acts, method, return_data=True)
                table = output_sheet(acts, current_data, act_callback=update_main, out_format=SheetFormat.NOTEBOOK)
                table_container.children = [table]

            finally:
                message.value = ""

        def main_combo_change(change):
            val = change.new
            if not val in main_acts_by_name :
                return
            act = main_acts_by_name[val]
            update_main(act)

        def base_combo_change(change):
            val = change.new
            if not val in base_acts_by_name :
                return
            act = base_acts_by_name[val]
            update_base(act)

        def comment_visibility_switch(container, change) :
            val = change.new
            container.layout.display = "block" if val else "none"

        def save_click(event) :
            acts = [main_act] if base_act is None else [main_act, base_act]

            filename = "_".join(set(_actName(act) for act in acts))
            filename += "-" + "_".join(_db_shortname(act.key[0]) for act in acts)
            filename += "_" + method_name(method)

            filename = secure_filename(filename) + ".xlsx"

            output_sheet(acts, current_data, out_format=SheetFormat.EXCEL, xls_filename=filename)

            message.value = "Saved to %s" % filename

        def tree_click(event):

            # Draw tree
            tree = show_tree(main_db, update_main)
            tree_container.children = [tree]


    # Bind events
    home_button.on_click(reset)
    back_button.on_click(pop_and_update)
    save_button.on_click(save_click)
    tree_button.on_click(tree_click)

    main_combo.observe(main_combo_change, "value")
    base_combo.observe(base_combo_change, "value")

    main_button_switch.observe(partial(comment_visibility_switch, main_comment), "value")
    base_button_switch.observe(partial(comment_visibility_switch, base_comment), "value")

    display(VBox([
        HBox([home_button, tree_button, back_button, save_button]),
        HBox([main_combo, main_button_switch]),
        main_comment,
        HBox([base_combo, base_button_switch]),
        base_comment,
        message,
        tree_container,
        table_container,
        output]))
