from lca_algebraic import resetDb, newActivity, newFloatParam, compute_impacts
from test.conftest import DataFixture, USER_DB, assert_impacts


def test_compute_impact_with_scenario(data: DataFixture):
    # Scnerio database are separated with '#'
    BG_DB1 = "bg#scen1"
    BG_DB2 = "bg#scen2"

    # Build two databases
    resetDb(BG_DB1, foreground=False)
    resetDb(BG_DB2, foreground=False)

    # Add twin BG activities
    bg_act = newActivity(BG_DB1, name="bg_act", unit="kg", exchanges={data.bio1: 1.0})

    # Same activitiy with different exchange in scenario 2
    bg_act2 = newActivity(BG_DB2, name="bg_act", unit="kg", exchanges={data.bio1: 2.0})

    # Build parametrized model on BG1
    p1 = newFloatParam("p1", default=1, min=0, max=1)

    fg_act = newActivity(USER_DB, "fg", "kg", {bg_act: p1})

    # Raw impacts without scenario : targeting scen1
    res = compute_impacts(models=fg_act, methods=data.ibio1, p1=1.0)

    assert_impacts(res, 1.0)

    # Impact for scenario 2
    res = compute_impacts(models=fg_act, methods=data.ibio1, scenario="scen2", p1=1.0)

    assert_impacts(res, 2.0)
