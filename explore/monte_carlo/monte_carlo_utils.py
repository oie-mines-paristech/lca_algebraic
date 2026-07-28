#!/usr/bin/env python
import bw2calc as bc
import bw2data as bd
import numpy as np
from bw2calc import JacobiGMRESLCA
from scipy import sparse
from scipy.sparse.csgraph import breadth_first_tree
from tqdm import tqdm

import lca_algebraic as agb

TEST_DB = "test"


class MonteCarloLCA(JacobiGMRESLCA):
    def __init__(self, *args, threshold=1e-4, **kwargs):
        super().__init__(*args, **kwargs)
        self.major_contributors = None
        self.suppliers = None
        self.threshold = threshold

    def _prepare_matrix(self) -> None:
        # Sparse cleanup is done once per matrix build, then reused.
        if getattr(self, "_matrix_prepared", False) and getattr(self, "_prepared_technosphere_matrix", None) is not None:
            return
        if not sparse.isspmatrix(self.technosphere_matrix):
            raise TypeError("technosphere_matrix must be a SciPy sparse matrix")

        # GMRES works best with canonical sparse structure.
        matrix = self.technosphere_matrix.tocsc(copy=False)

        # if self.major_contributors is not None :
        #    # print("filtering")
        #    # If major contributors already computed (done once) then filter it
        #    print("Size matrix before", matrix.nnz)
        #    matrix = self.major_contributors @ matrix @ self.major_contributors
        #    print("Size matrix after", matrix.nnz)

        matrix.sum_duplicates()
        matrix.eliminate_zeros()
        matrix.sort_indices()

        self._prepared_technosphere_matrix = matrix
        self._matrix_prepared = True

    def compute_suppliers(self):
        """Compute list of suppliers whose net supply > threshold * tortal-impact"""

        self.lci()
        self.load_lcia_data()

        g = self.characterization_matrix.diagonal() * self.biosphere_matrix

        # Initial demand
        impacts = self.supply_array * g

        total_impact = np.sum(impacts)

        # Sort impacts
        sorted_idx = np.argsort(-impacts)
        impacts = impacts[sorted_idx]

        cum_rel_impacts = np.cumsum(impacts) / total_impact

        # Only keep acitvities accounting for (1-threshold)*total_impact
        mask = cum_rel_impacts < 1 - self.threshold

        suppliers = sorted_idx[mask]

        demand_idx = list(self.dicts.product[demand] for demand in self.demand)

        reachable_demand = multi_source_bfs(self.technosphere_matrix, suppliers)
        reachable_supply = multi_source_bfs(self.technosphere_matrix.T, demand_idx)
        intersec = np.intersect1d(reachable_demand, reachable_supply)

        print("Main suppliers", len(suppliers))
        print("Reachable demand", len(reachable_demand))
        print("Reachable supply", len(reachable_supply))
        print("Intersect", len(intersec))

    def compute_major_contributors(self, max_iter=50):
        """
        Recusrively compute contribution of each activity for a the current supply.
        Then return a mask of cumulative impacts reaching more than a threshold (realtive to total impact)


        A : technological matrix (n x n) creuse (format CSC recommandé)
        x : supply (vecteur n)
        g : direct unitaary impacts per activity (vecteur n)
        """

        # Solve supply once
        self.lci()
        self.load_lcia_data()

        n = self.technosphere_matrix.shape[0]
        A = sparse.identity(n, format="csr") - self.technosphere_matrix

        # Net impact of each activity
        g = self.characterization_matrix.diagonal() * self.biosphere_matrix

        # Initial demand
        x = self.demand_array.copy()

        total_impact = np.sum(self.supply_array * g)  # impact total final (g^T x)

        for k in range(max_iter):
            # Add incrementally adds a supply level
            delta = A @ x  # A^T * lambda_k

            # Residual impacts aded at this step is negligible ?
            if np.max(delta * g) < self.threshold * total_impact:
                break

            x += delta

        # Recursive impacts per activity
        impacts = x * g

        mask = (impacts > (self.threshold * total_impact)).astype(np.uint8)

        print(f"Major contributors {np.sum(mask)} / {self.technosphere_matrix.shape[0]}")
        print(f"Nb iterations {k}")

        self.major_contributors = sparse.diags(mask, format="csr")

        return mask


def MonteCarlo(
    demand: dict[bd.Node, float],
    method: tuple,
    iterations: int = 20,
    use_jacobi=False,
    use_guess=True,
) -> dict:
    # Create all possible demands
    common_kwargs = {"demand": demand, "method": method, "use_distributions": True, "seed_override": 42}

    # LcaClass = bc.PartitionedMonteCarloLCA
    if use_jacobi:
        lca = MonteCarloLCA(use_guess=use_guess, rtol=1e-8, threshold=1e-2, **common_kwargs)
    else:
        # Create a single LCA object and use uncertainty distributions
        lca = bc.LCA(**common_kwargs)

    # Compute major contributors once
    if True:
        lca.lci()
    else:
        suppliers = lca.compute_suppliers()
        print("Suppliers", len(suppliers))

    # Create a list of characterization matrices
    lca.switch_method(method)
    c_matrix = lca.characterization_matrix.copy()

    res = np.empty(iterations, dtype=float)

    # Do one monte carlo iteration for all functional units and impact categories
    for i in tqdm(range(iterations)):
        # Resample all matrices
        next(lca)

        # Convert to integer ids instead of `bd.Node` objects
        lca.lci()

        res[i] = (c_matrix * lca.inventory).sum()

    return res


def init():
    bd.projects.set_current("test_monte_carlo")


def perf_test():
    init()

    pv_act = agb.findTechAct("electricity production, photovoltaic, 570kWp open ground installation, multi-Si", loc="CN-FJ")

    print(pv_act)

    climate = ("ecoinvent-3.11", "EF v3.1", "climate change", "global warming potential (GWP100)")

    # Required to fix the bug in trying to clear pypardiso first
    res = MonteCarlo(
        demand={pv_act: 1.0},
        method=climate,
        iterations=100,
        use_jacobi=False)

    res = MonteCarlo(
        demand={pv_act: 1.0},
        method=climate,
        iterations=100,
        use_jacobi=True,
        use_guess=True)

    print(np.mean(res), np.std(res))


def add_super_source(G, seeds):
    n = G.shape[0]

    # Ligne : super-source -> seeds
    row = sparse.csr_array(
        (
            np.ones(len(seeds), dtype=G.dtype),
            (np.zeros(len(seeds), dtype=np.int32), seeds),
        ),
        shape=(1, n),
    )

    # Colonne vide
    col = sparse.csr_array((n, 1), dtype=G.dtype)

    # Coin inférieur droit (0)
    zero = sparse.csr_array((1, 1), dtype=G.dtype)

    G_ext = sparse.block_array(
        [
            [G, col],
            [row, zero],
        ],
        format="csr",
    )

    return G_ext, n


def multi_source_bfs(A, seeds):
    """Find reacheable nodes from multiple nodes"""

    G = (A != 0).astype(int)

    G_ext, super_node = add_super_source(G, seeds)

    tree = breadth_first_tree(G_ext, i_start=super_node, directed=True)

    reachable = np.unique(tree.indices)

    # enlever la super-source
    reachable = reachable[reachable != super_node]

    return reachable


if __name__ == "__main__":
    perf_test()
