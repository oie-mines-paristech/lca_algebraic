from collections import defaultdict
from typing import Dict

from sympy import Piecewise, simplify

from lca_algebraic import ParamDef, newActivity, warn


def _segments_to_piecewise(param, segments):
    conds = []
    for start, end, val in segments:
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
     Creates a linear virtual activity being a linear interpolation between several activities,
     based on the values of a given parameter.

     This is useful to produce a continuous parametrized activity based on the scale of the system,
     given that you have dicrete activities coresponding to
     discrete values of the parameter.

    :param db_name: Name of user DB (string)
    :param act_name: Name of the new activity
    :param param: Parameter to use [ParamDef]
    :param act_per_value : Dictionnary of value => activitiy [Dict]
    :param add_zero: If True add the "Zero" point to the data. Usefull for linear interoplation of a single activity / point
    :return: the new activity
    """

    # Add "Zero" to the list
    act_per_value = act_per_value.copy()

    if add_zero:
        act_per_value[0.0] = None

    # List of segments : triplet of (start, end, expression)
    segments = defaultdict(list)

    # Transform to sorted list of value => activity
    sorted_points = sorted(act_per_value.items(), key=lambda item: item[0])
    for i, (curr_val, curr_act) in enumerate(sorted_points):
        if i >= len(sorted_points) - 1:
            continue

        # Next val and act
        next_val, next_act = sorted_points[i + 1]

        # Boundaries of segment : none if first / last point
        start = curr_val if i > 0 else None
        end = next_val if i < (len(sorted_points) - 2) else None

        # Add segment for current activity
        segments[curr_act].append(
            [start, end, (param - next_val) / (curr_val - next_val)]
        )  # Will equal 1 at current point and 0 at next point

        # Add segment for next activity
        segments[next_act].append(
            [start, end, (param - curr_val) / (next_val - curr_val)]
        )  # Will equal 0 at current point and 1 at next point

    # Transform segments into piecewize expressions
    exchanges = {act: _segments_to_piecewise(param, segs) for act, segs in segments.items() if act is not None}

    # Find unit
    units = list(act["unit"] for act in exchanges.keys())
    same_unit = all(x == units[0] for x in units)

    if not same_unit:
        warn("Warning : units of activities should be the same : %s" % str(units))

    # Create act
    new_act = newActivity(db_name=db_name, name=act_name, unit=units[0], exchanges=exchanges)

    return new_act
