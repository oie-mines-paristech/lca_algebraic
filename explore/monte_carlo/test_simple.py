from monte_carlo_utils import TEST_DB, MonteCarlo, init

import lca_algebraic as agb


def setup_test():
    init()
    agb.resetDb(TEST_DB)

    co2 = agb.findBioAct("Carbon dioxide, fossil", categories=("air",))
    co2_bis = agb.findBioAct("Carbon dioxide, fossil", categories=("air", "urban air close to ground"))

    # B with no direct emission but emsision via C
    C = agb.newActivity(db_name=TEST_DB, name="c", exchanges={co2: 1.0}, unit="kg")
    B = agb.newActivity(db_name=TEST_DB, name="b", exchanges={C: 1.0}, unit="kg")

    # D has neglibile impacts
    D = agb.newActivity(db_name=TEST_DB, name="d", exchanges={co2_bis: 0.000001}, unit="kg")

    A = agb.newActivity(
        db_name=TEST_DB,
        name="a",
        unit="kg",
        exchanges={
            B: 1,
            D: 1,
        },
    )

    return A


def test_example():
    A = setup_test()

    climate = ("ecoinvent-3.11", "EF v3.1", "climate change", "global warming potential (GWP100)")

    res = MonteCarlo(demand={A: 1.0}, method=climate, iterations=1, use_jacobi=True, use_guess=True)

    print(res)
