from collections import defaultdict
from typing import Dict

import pint
from sympy import Piecewise, simplify

from lca_algebraic import ParamDef, newActivity, warn

from .settings import Settings
from .units import parse_db_unit


def _segments_to_piecewise(act, param, segments):
    conds = []
    for start, end, val in segments:
        val = act._transform_unit(val, act["unit"])
        # BUG ? Cannot store quantity in Piecewise.
        if isinstance(val, pint.Quantity):
            val = val.m
        cond = True
        if start is not None:
            cond = cond & (param >= start)
        if end is not None:
            cond = cond & (param < end)
        conds.append((val, simplify(cond)))

    return Piecewise(*conds, (0, True))


def interpolate_activities(
    db_name,
    act_name,
    param: ParamDef,
    act_per_value: Dict,
    add_zero=False,
):
    """
    Creates a virtual activity being a linear interpolation between several activities,
    based on the value of a parameter.

    The result activity is a piecewize linear function of input activities.

    This is useful to produce a continuous parametrized activity based on the scale of the system,
    Given that you have discrete activities corresponding to discrete values of the parameter.

    Parameters
    ----------
    db_name:
        Name of user DB (string)
    act_name:
        Name of the new activity
    param:
        Parameter controlling the interpolation
    act_per_value:
        Dictionnary of value => Activitiy [Dict]
        Notes:
            * Acticity may be None, it's equivalent to an activity without exchanges.
            * If parameter is lower than the lowest bound in act_per_value
                        then this activity is equal to the activity at lower bound
            * If parameter is higher than the highest bound in act_per_value
                        then this activity is equal to the activity at highest bound

    add_zero:
        If True add the "Zero" point to the data, i.e. add act_per_value[0.0] = None
        Useful for linear interpolation of a single activity / point

    Returns
    -------
    The new activity

    Examples
    --------
    >>> interpolated_inverter = interpolate_activities(
    >>>     db_name=USER_DB,
    >>>     act_name="interpolated_inverter",
    >>>     param=power_param, #power parameter, in kW
    >>>     act_per_value={ # Those are activities found in the background database
    >>>         .5: inverter_500W,
    >>>         2: inverter_2KW,
    >>>         50:inverter_50KW,
    >>>         100:inverter_100KW},
    >>>     add_zero=True):
    """

    # Add "Zero" to the list
    act_per_value = act_per_value.copy()

    if add_zero:
        act_per_value[0.0 * next(iter(act_per_value))] = None

    # Find unit
    units = [act["unit"] for act in act_per_value.values() if act is not None]
    same_unit = all(x == units[0] for x in units)

    if not same_unit:
        warn("Warning : units of activities should be the same : %s" % str(units))

    # List of segments : triplet of (start, end, expression)
    segments = defaultdict(list)

    # Transform to sorted list of value => activity
    sorted_points = list(sorted(act_per_value.items(), key=lambda item: item[0]))
    sorted_points = [(None, sorted_points[0][1])] + sorted_points + [(None, sorted_points[-1][1])]
    for (l_val, l_act), (r_val, r_act) in zip(sorted_points[0:-1], sorted_points[1:]):
        # Add segment for current activity

        # Left bound, right bound or same activity on left or right
        if l_act == r_act:
            unit_amount = 1.0
            if Settings.units_enabled:
                unit_amount |= parse_db_unit(units[0])
            segments[l_act].append([l_val, r_val, unit_amount])
            continue

        segments[l_act].append(
            [l_val, r_val, (param - r_val) / (l_val - r_val)]
        )  # Will equal 1 at current point and 0 at next point

        # Add segment for next activity
        segments[r_act].append(
            [l_val, r_val, (param - l_val) / (r_val - l_val)]
        )  # Will equal 0 at current point and 1 at next point

    # Transform segments into piecewize expressions
    exchanges = {act: _segments_to_piecewise(act, param, segs) for act, segs in segments.items() if act is not None}

    if Settings.units_enabled:
        exchanges = {act: amount | parse_db_unit(act["unit"]) for act, amount in exchanges.items()}

    # Create act
    new_act = newActivity(db_name=db_name, name=act_name, unit=units[0], exchanges=exchanges)

    return new_act
