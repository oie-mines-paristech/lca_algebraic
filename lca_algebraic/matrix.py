from collections import defaultdict
from multiprocessing import Pool

import numpy as np
import pandas as pd
from scipy.sparse import csr_array
from scipy.sparse.csgraph import connected_components
from sympy import (
    ImmutableMatrix,
    ImmutableSparseMatrix,
    Matrix,
    MatrixBase,
    MutableSparseMatrix,
)


def _to_structural_matrix(mat: MatrixBase):
    """Convert a Sympy matrix to a structural numpy matrix (0 and 1)"""
    n, m = mat.shape
    structural = np.zeros((n, m), dtype=int)
    for i in range(n):
        for j in range(m):
            elem = mat[i, j]
            if not (elem == 0 or elem == 0.0):
                structural[i, j] = 1
    return structural


def _find_sccs(mat: MatrixBase):
    """Find Strongly Connected Components (SCCs) in the matrix graph"""
    structural = _to_structural_matrix(mat)
    sparse_graph = csr_array(structural)
    n_components, labels = connected_components(csgraph=sparse_graph, directed=True, connection="strong", return_labels=True)

    # Group nodes by SCC
    sccs = {}
    for node_idx, label in enumerate(labels):
        if label not in sccs:
            sccs[label] = []
        sccs[label].append(node_idx)

    return sccs, labels


def _extract_submatrix(mat: MatrixBase, indices):
    """Extract a submatrix given row/column indices"""
    n = len(indices)
    submat = MutableSparseMatrix.zeros(n, n)
    for i, row_idx in enumerate(indices):
        for j, col_idx in enumerate(indices):
            submat[i, j] = mat[row_idx, col_idx]
    return submat


def _invert_scc_worker(args):
    """Worker function to invert a single SCC (for multiprocessing)"""
    label, nodes, mat_data = args
    # Reconstruct the submatrix from the data
    scc_submat = Matrix(mat_data)
    scc_size = len(nodes)
    identity = Matrix.eye(scc_size)
    scc_inv = scc_submat.LUsolve(identity)
    return (label, nodes, scc_inv)


def invert(mat: MatrixBase, min_scc_size_for_parallel=5):
    """

    Fast invert of large sparse matrix, handling SCCs (Strongly Connected Components).

    Args:
        mat: Matrix to invert
        min_scc_size_for_parallel: Minimum SCC size to invert in parallel (0 = no parallelism)
    """
    n, _ = mat.shape

    # Find SCCs
    sccs, labels = _find_sccs(mat)

    # Check if entire matrix is a single SCC
    if len(sccs) == 1 and len(list(sccs.values())[0]) == n:
        # Entire matrix is one SCC, use LUsolve directly
        identity = Matrix.eye(n)
        return mat.LUsolve(identity)

    # Identify which nodes are in SCCs of size > 1
    scc_nodes = set()
    scc_inverses = {}  # Map from SCC label to (node_list, inverse_matrix)
    scc_node_to_scc = {}  # Map from node index to its SCC label
    scc_node_to_pos = {}  # Map from node index to its position in SCC

    # Separate SCCs by size for potential parallel processing
    sccs_to_invert = [(label, nodes) for label, nodes in sccs.items() if len(nodes) > 1]

    if sccs_to_invert and min_scc_size_for_parallel > 0:
        # Use parallel processing for SCCs >= min_scc_size_for_parallel
        parallel_sccs = [(label, nodes) for label, nodes in sccs_to_invert if len(nodes) >= min_scc_size_for_parallel]
        sequential_sccs = [(label, nodes) for label, nodes in sccs_to_invert if len(nodes) < min_scc_size_for_parallel]

        # Process parallel SCCs
        if parallel_sccs:
            # Prepare data for workers (convert to list for pickling)
            parallel_args = []
            for label, nodes in parallel_sccs:
                scc_submat = _extract_submatrix(mat, nodes)
                mat_data = [[scc_submat[i, j] for j in range(len(nodes))] for i in range(len(nodes))]
                parallel_args.append((label, nodes, mat_data))

            # Invert in parallel
            with Pool() as pool:
                parallel_results = pool.map(_invert_scc_worker, parallel_args)

            for label, nodes, scc_inv in parallel_results:
                scc_inverses[label] = (nodes, scc_inv)
                for pos, node in enumerate(nodes):
                    scc_nodes.add(node)
                    scc_node_to_scc[node] = label
                    scc_node_to_pos[node] = pos

        # Process sequential SCCs
        for label, nodes in sequential_sccs:
            scc_submat = _extract_submatrix(mat, nodes)
            scc_size = len(nodes)
            identity = Matrix.eye(scc_size)
            scc_inv = scc_submat.LUsolve(identity)
            scc_inverses[label] = (nodes, scc_inv)
            for pos, node in enumerate(nodes):
                scc_nodes.add(node)
                scc_node_to_scc[node] = label
                scc_node_to_pos[node] = pos
    else:
        # Sequential processing
        for label, nodes in sccs_to_invert:
            scc_submat = _extract_submatrix(mat, nodes)
            scc_size = len(nodes)
            identity = Matrix.eye(scc_size)
            scc_inv = scc_submat.LUsolve(identity)
            scc_inverses[label] = (nodes, scc_inv)
            for pos, node in enumerate(nodes):
                scc_nodes.add(node)
                scc_node_to_scc[node] = label
                scc_node_to_pos[node] = pos

    # Initialize result
    res = MutableSparseMatrix.zeros(*mat.shape)
    done_rows = set()

    def _process_row(row_idx):
        """Recursively process a row"""
        if row_idx in done_rows:
            return res[row_idx, :]

        # If this node is in a non-trivial SCC, use the precomputed inverse
        if row_idx in scc_nodes:
            scc_label = scc_node_to_scc[row_idx]
            scc_nodes_list, scc_inv = scc_inverses[scc_label]
            scc_pos = scc_node_to_pos[row_idx]

            # First, copy the row from the SCC inverse for columns within the SCC
            for j, col_node in enumerate(scc_nodes_list):
                res[row_idx, col_node] = scc_inv[scc_pos, j]

            # Now handle dependencies from this SCC to nodes outside the SCC
            # For each node outside the SCC that this SCC depends on
            for col_idx in range(n):
                if col_idx in scc_nodes_list:
                    continue

                # Check if any node in the SCC depends on col_idx
                # If so, we need to propagate this dependency
                for scc_node_idx, scc_node in enumerate(scc_nodes_list):
                    elem = mat[scc_node, col_idx]
                    if not (elem == 0 or elem == 0.0):
                        # Get the row for col_idx recursively
                        col_row = _process_row(col_idx)
                        contribution = scc_inv[scc_pos, scc_node_idx] * mat[scc_node, col_idx] * col_row
                        res[row_idx, :] -= contribution

            done_rows.add(row_idx)
            return res[row_idx, :]

        # Standard recursive processing for non-SCC nodes (or SCC of size 1)
        # Factor is 1 / self_value
        factor = 1 / mat[row_idx, row_idx]

        res[row_idx, row_idx] = factor

        for col_idx in range(n):
            # Self col or null col => exit
            if col_idx == row_idx:
                continue
            elem = mat[row_idx, col_idx]
            if elem == 0 or elem == 0.0:
                continue

            # Recursively get the row to add, multiply by current cell and factor
            row_to_add = mat[row_idx, col_idx] * _process_row(col_idx)
            if (factor != 1) and (factor != 1.0):
                row_to_add = row_to_add * factor

            # Accumulate to output (substract)
            res[row_idx, :] -= row_to_add

        done_rows.add(row_idx)
        return res[row_idx, :]

    for i in range(n):
        _process_row(i)

    return res


class ActMatrix(defaultdict):
    def __init__(self):
        super().__init__(lambda: 0.0)
        self._col_acts = list()
        self._row_acts = list()

    def __getitem__(self, key):
        row_act, col_act = key
        self.add_row(row_act)
        self.add_col(col_act)
        return super().__getitem__(key)

    def __setitem__(self, key, value):
        row_act, col_act = key
        self.add_row(row_act)
        self.add_col(col_act)
        return super().__setitem__(key, value)

    def add_col(self, col_act):
        if col_act not in self._col_acts:
            self._col_acts.append(col_act)
            self._col_acts.sort()

    def add_row(self, row_act):
        if row_act not in self._row_acts:
            self._row_acts.append(row_act)
            self._row_acts.sort()

    def demand_vector(self, act, value=1.0):
        """Generate vector of len (rows) with zeros and 1 only at the index of act"""
        act_idx = self.row_acts().index(act)
        res = [0.0] * len(self._row_acts)
        res[act_idx] = value
        return ImmutableMatrix([res])

    def row_acts(self):
        return self._row_acts

    def cols_acts(self):
        return self._col_acts

    def shape(self):
        return len(self._row_acts), len(self._col_acts)

    def to_sympy(self):
        """Return an immutable sympy matrix"""
        rows = list()
        for row_act in self.row_acts():
            row = list()
            rows.append(row)
            for col_act in self.cols_acts():
                row.append(self.get((row_act, col_act), 0.0))
        return ImmutableSparseMatrix(rows)

    def to_dataframe(self):
        res = dict()
        for row_act in self.row_acts():
            res[str(row_act)] = {str(col_act): self[(row_act, col_act)] for col_act in self.cols_acts()}
        return pd.DataFrame(res)

    def __repr__(self):
        return self.to_dataframe().__repr__()
