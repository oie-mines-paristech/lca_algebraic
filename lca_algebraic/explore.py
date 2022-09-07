from functools import partial
from types import SimpleNamespace

import brightway2 as bw
import ipysheet
import ipywidgets as widgets
from IPython.display import display
from ipyevents import Event
from ipywidgets import VBox, Button, Output, Combobox, Layout, ToggleButton, HBox, HTML
from sympy import Basic

from lca_algebraic import ActivityExtended, multiLCAAlgebric, getActByCode
from lca_algebraic.base_utils import _actName
from lca_algebraic.lca import _createTechProxyForBio
from lca_algebraic.params import _completeParamValues
from lca_algebraic.stats import _round_expr


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


_SMALL_DIFF_COLOR = "LemonChiffon"
_HIGH_DIFF_COLOR = "Coral"
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


def compare_impacts(
        activities : ActivityExtended,
        method,
        default_params=True,
        diff_impact=0.1,
        act_callback=None,
        **params):
    """
    Advanced version of #printAct()

    Displays all exchanges of one or several activities and their impacts.
    If parameter values are provided, formulas will be evaluated accordingly.
    If two activities are provided, they will be shown side by side and compared.

    :param activities: One activity, or couple of activities
    :param method: Impact method
    :param default_params: If true, replace params with default values in formulas
    :param diff_impact: Min minimum impact difference to highlight, as part of total impact
    :param act_callback: Callback when clicking on an activity link
    :param params: Extra parameter values to replace in formulas
    """

    if not isinstance(activities, (list, tuple)):
        activities = [activities]

    # List of dict[exchange=>attributes]
    # As many items as source activities
    datas = list()
    nact = len(activities)

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

        res = multiLCAAlgebric(all_acts, [method], **params)
        impacts = res[res.columns.tolist()[0]].to_list()

        impact_by_act = {act : value for act, value in zip(all_acts, impacts)}

        # Add impacts to data
        for key, attrs in data.items() :
            amount = attrs.amount
            act = inputs_by_ex_name[key]
            impact = impact_by_act[act]
            attrs.impact = amount * impact

    # All exch
    all_ex_names = list(set(ex for data in datas for ex in data.keys()))

    headers = ["exchange", "unit"]

    # Exchange, Unit, Inputs, Amounts, Impacts
    widths = [170, 50] \
             + [210] * nact \
             + [70] * nact \
             + [70] * nact

    def _numerate(str_val) :
        return [str_val] + ["%s #%d" % (str_val, i) for i in range(1, nact)]

    headers += _numerate("input") + _numerate("amount") + _numerate("impact")

    # Compute sum of impacts for first activity
    sum_impact = sum(attrs.impact for attrs in datas[0].values())

    # To sheet
    sheet  = ipysheet.sheet(
        rows=len(all_ex_names),
        columns=len(headers),
        column_headers=headers,
        row_headers=False,
        stretch_headers="none",
        column_width=widths)

    # Loop on exchanges
    for i_ex, ex_name in enumerate(all_ex_names) :

        ipysheet.cell(row=i_ex, column=0, value=HTML("<b class='trunc' title='%s'>%s</b>" % (ex_name, ex_name)))

        # Loop on type of columns
        for column in ["input", "amount", "impact"] :

            # On or two cells, one per activity
            values = []
            cells = []

            # Loop on activity (1 or 2)
            for i_act, (act, data) in enumerate(zip(activities, datas)) :

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
                if attrs.unit != "-" :
                    ipysheet.cell(row=i_ex, column=1, value=attrs.unit, style=_RIGHT_BORDER, read_only=True)

                # Switch on column type

                if column == "input" : # Input activity

                    if attrs.input == "_":
                        input = "-"
                        values.append(input)
                    else:
                        input_name = "<a href='#'>%s</a>" % attrs.input
                        loc_html = "<span class='loc'> %s </span>" % attrs.loc
                        db_html = "&nbsp;<span class='db%d'> %s </span>" % (_db_idx(attrs.db), _db_shortname(attrs.db))
                        html = input_name + "<br/>" + loc_html + db_html

                        values.append(html)

                        input = HTML(html)

                        # Link the callback to click on it
                        if act_callback is not None :

                            sub_act = getActByCode(attrs.db, attrs.code)

                            def cb(sub_act, event):
                                act_callback(sub_act)

                            event = Event(source=input, watched_events=['click'])
                            event.on_dom_event(partial(cb, sub_act))

                    cells.append(input)

                elif column == "amount": # Amount and formula

                    values.append(attrs.amount)

                    amount_repr = _sci_repr(attrs.amount)
                    if attrs.formula :

                        amount_repr = HTML("<abbr title='%s'><b>%s</b></abbr>" % (str(attrs.formula), amount_repr))

                    cells.append(amount_repr)

                else:

                    values.append(attrs.impact)
                    cells.append(_sci_repr(attrs.impact))

            # Colorize background according to Diff
            color = None

            if (nact > 1) and values[0] != values[1] :

                # For impact, highlight differently difference higher than a given share of total impact
                if column == "impact":
                    diff = abs(
                        _none_float(values[1]) -
                        _none_float(values[0]))
                    rel_diff = diff / sum_impact

                    color = _HIGH_DIFF_COLOR if rel_diff > diff_impact else _SMALL_DIFF_COLOR

                else:
                    color = _SMALL_DIFF_COLOR


            # Display cells for this column
            for icell, cell in enumerate(cells) :
                icol = 2 + nact * (0 if column == "input" else 1 if column == "amount" else 2)

                style= _RIGHT_BORDER if (icell == nact - 1) else None

                ipysheet.cell(
                    row=i_ex,
                    column=icol + icell,
                    style=style.copy() if style else None,
                    background_color=color,
                    value=cell,
                    read_only=True)

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
    history = []

    home_button = Button(icon="home")
    back_button = Button(icon="arrow-left")
    save_button = Button(icon="save")

    output = Output()

    main_act = None
    base_act = None

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

            # Update names in combo
            main_combo.value = "" if main_act is None else _actName(main_act)
            base_combo.value = "" if base_act is None else _actName(base_act)

            # Update comments
            main_comment.value = "" if main_act is None else (main_act["comment"] if "comment" in main_act else "")
            base_comment.value = "" if base_act is None else (base_act["comment"] if "comment" in base_act else "")

            message.value = "<h4>Computing ...</h4>"

            try:
                acts = main_act if base_act is None else [main_act, base_act]
                table = compare_impacts(acts, method, act_callback=update_main)
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
            table = table_container.children[0]
            df = ipysheet.to_dataframe(table)
            filename = "out.xls"
            df.to_excel(filename)
            print("Saved to %s" % filename)


    # Bind events
    home_button.on_click(reset)
    back_button.on_click(pop_and_update)
    save_button.on_click(save_click)

    main_combo.observe(main_combo_change, "value")
    base_combo.observe(base_combo_change, "value")

    main_button_switch.observe(partial(comment_visibility_switch, main_comment), "value")
    base_button_switch.observe(partial(comment_visibility_switch, base_comment), "value")

    display(VBox([
        HBox([home_button, back_button, save_button]),
        HBox([main_combo, main_button_switch]),
        main_comment,
        HBox([base_combo, base_button_switch]),
        base_comment,
        message,
        table_container,
        output]))
