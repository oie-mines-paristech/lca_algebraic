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
        Dictionnary of value => activitiy [Dict]

    add_zero:
        If True add the "Zero" point to the data.
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
