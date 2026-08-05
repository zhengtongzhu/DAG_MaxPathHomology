from __future__ import annotations

import numpy as np
import networkx as nx
import sympy as sp

from scipy import sparse as sps
from maxph_matrix import (node_symbol, basis_mat_sparse, 
                          hstack_sparse, identity_sparse,
                          null_identity_sparse, rref_null_sparse,
                          row_del_sparse, sparse_density,
                          sparse_memory_mib)
from dag_preprocess import dag_process

def _initial_cycle_basis_sparse(dim: int) -> sps.csr_matrix:
    """
    Return the sparse matrix [-1; I_dim] used at the first layer.

    The returned matrix has shape (dim + 1, dim) and exactly 2 * dim
    nonzero entries when dim > 0.
    """
    if dim < 0:
        raise ValueError("dim must be nonnegative.")

    if dim == 0:
        return sps.csr_matrix((1, 0), dtype=np.int32)

    columns = np.arange(dim, dtype=np.int32)
    rows = np.concatenate(
        (
            np.zeros(dim, dtype=np.int32),
            np.arange(1, dim + 1, dtype=np.int32),
        )
    )
    repeated_columns = np.concatenate((columns, columns))
    values = np.concatenate(
        (
            -np.ones(dim, dtype=np.int32),
            np.ones(dim, dtype=np.int32),
        )
    )

    return sps.coo_matrix(
        (values, (rows, repeated_columns)),
        shape=(dim + 1, dim),
        dtype=np.int32,
    ).tocsr()

def _initial_symbolic_basis(
    nodes: list,
) -> sp.SparseMatrix:
    """Return [x_1-x_0, ..., x_m-x_0] as a 1-by-(m-1) sparse matrix."""
    if len(nodes) <= 1:
        return sp.SparseMatrix(1, 0, {})

    symbols = {node: node_symbol(node) for node in nodes}
    first_symbol = symbols[nodes[0]]

    entries = {
        (0, column): symbols[node] - first_symbol
        for column, node in enumerate(nodes[1:])
    }

    return sp.SparseMatrix(1, len(nodes) - 1, entries)

def _print_sparse_state(
    label: str,
    matrix,
) -> None:
    """Print shape, nnz, density, and CSR storage for diagnostics."""
    matrix_csr = matrix.tocsr()
    print(
        f"{label}: shape={matrix_csr.shape}, "
        f"nnz={matrix_csr.nnz}, "
        f"density={sparse_density(matrix_csr):.6e}, "
        f"storage={sparse_memory_mib(matrix_csr):.3f} MiB"
    )

def max_path_homology(
    G: nx.DiGraph,
    calculate_basis: bool = False,
    report_sparsity: bool = False,
) -> tuple[int, int, list[list[sp.Expr]] | None]:
    """
    Compute maximal path homology using sparse numerical matrices.

    Parameters
    ----------
    G:
        A finite NetworkX directed acyclic graph.

    calculate_basis:
        Whether to construct an explicit symbolic homology basis.

    report_sparsity:
        Print sparse-matrix diagnostics at every recursive layer.

    Returns
    -------
    lp:
        Maximum directed path length of G.

    betti:
        The lp-dimensional Betti number of path homology.

    basis:
        Explicit basis if requested; otherwise None.

    Notes
    -----
    The numerical recursion uses SciPy CSR matrices throughout. No large
    dense identity, zero, hstack, vstack, or null-space matrix is created.
    Symbolic basis tracking uses SymPy SparseMatrix, but explicit symbolic
    generators can still be expensive when the output basis itself is large.
    """
    if not isinstance(G, nx.DiGraph):
        raise TypeError("G must be a NetworkX DiGraph.")

    (subgraph_dict, N, lp, num_subgraphs, graph_list) = dag_process(G)

    betti = 0
    basis: list[list[sp.Expr]] | None = ([] if calculate_basis else None)

    if lp == 0:
        isolated_nodes = list(G.nodes())
        betti = len(isolated_nodes) - 1

        if not calculate_basis or len(isolated_nodes) <= 1:
            return lp, betti, basis

        symbols = {node: node_symbol(node) for node in isolated_nodes}
        first_node = isolated_nodes[0]
        basis = [[
            symbols[node] - symbols[first_node]
            for node in isolated_nodes[1:]
        ]]
        return lp, betti, basis

    if num_subgraphs == 0:
        return lp, betti, basis

    for subgraph_index in range(num_subgraphs):
        if any(layer_size == 1 for layer_size in N[subgraph_index]):
            continue

        dim = N[subgraph_index][0] - 1
        V = _initial_cycle_basis_sparse(dim)
        partition = {node: 1 for node in subgraph_dict[subgraph_index][0]}

        if calculate_basis:
            basis_iteration = _initial_symbolic_basis(list(partition.keys()))

        if report_sparsity:
            _print_sparse_state(f"subgraph={subgraph_index}, layer=0, V", V)

        for layer_index in range(1, lp + 1):
            if dim == 0:
                break

            column_blocks: list[sps.csr_matrix] = []
            next_partition: dict[str, int] = {}
            has_identity_block = False
            identity_block_index = -1
            V_row_count = V.shape[0]

            # Every full local cycle space uses the same -I_dim block.
            negative_identity = -identity_sparse(dim)

            for node_index, node in enumerate(subgraph_dict[subgraph_index][layer_index]):
                predecessors = list(graph_list[subgraph_index].predecessors(node))
                predecessor_row_count = sum(partition.get(predecessor, 0) for predecessor in predecessors)

                if predecessor_row_count == V_row_count:
                    next_partition[node] = dim
                    column_blocks.append(negative_identity)

                    if not has_identity_block:
                        has_identity_block = True
                        identity_block_index = node_index
                    continue

                reduced_matrix = row_del_sparse(V, partition, predecessors)
                A_x = rref_null_sparse(reduced_matrix)

                # A_x has shape (dim, local_nullity).
                local_nullity = A_x.shape[1]

                if local_nullity == 0:
                    next_partition[node] = 0
                    continue

                next_partition[node] = local_nullity

                if local_nullity == dim:
                    column_blocks.append(negative_identity)

                    if not has_identity_block:
                        has_identity_block = True
                        identity_block_index = node_index
                else:
                    column_blocks.append(-A_x)

            if not column_blocks:
                dim = 0
                break

            A = hstack_sparse(column_blocks)

            if report_sparsity:
                _print_sparse_state(f"subgraph={subgraph_index}, layer={layer_index}, A", A)

            if has_identity_block:
                V, dim = null_identity_sparse(next_partition, identity_block_index, A, dim)
            else:
                V = rref_null_sparse(A)
                dim = V.shape[1]

            if report_sparsity:
                _print_sparse_state(f"subgraph={subgraph_index}, layer={layer_index}, V", V)

            if dim == 0:
                break

            if calculate_basis:
                basis_update = basis_mat_sparse(column_blocks, next_partition, V)
                basis_iteration = (basis_iteration * basis_update)

            partition = next_partition

        betti += dim

        if calculate_basis and dim != 0:
            assert basis is not None
            basis.extend(basis_iteration.tolist())

    return lp, betti, basis

# Example usage
if __name__ == "__main__":
    import random
    import time

    nodes_per_layer = [4, 10, 10, 10]
    density = 0.9
    seed = 1
    calculate_basis = True

    random_generator = random.Random(seed)

    graph = nx.DiGraph()
    layers: list[list[str]] = []

    # Construct the vertices layer by layer.
    for layer_index, layer_size in enumerate(nodes_per_layer):
        layer_nodes = [
            f"K{layer_index}_{node_index}"
            for node_index in range(layer_size)
        ]

        layers.append(layer_nodes)
        graph.add_nodes_from(layer_nodes)

    edge_counts: list[int] = []

    # Add floor(density * n_p * n_(p+1)) randomly selected
    # edges between every pair of consecutive layers.
    for layer_index in range(len(layers) - 1):
        source_layer = layers[layer_index]
        target_layer = layers[layer_index + 1]

        possible_edges = [(source, target) for source in source_layer for target in target_layer]

        random_generator.shuffle(possible_edges)

        number_of_edges = int(density * len(source_layer) * len(target_layer))
        selected_edges = possible_edges[:number_of_edges]

        graph.add_edges_from(selected_edges)
        edge_counts.append(number_of_edges)

    expected_depth = len(nodes_per_layer) - 1
    actual_depth = nx.dag_longest_path_length(graph)

    print("=" * 72)
    print("Recursive-algorithm example 1")
    print("-" * 72)
    print(f"nodes_per_layer = {nodes_per_layer}")
    print(f"density         = {density}")
    print(f"seed            = {seed}")
    print(f"edge_counts     = {edge_counts}")
    print(f"total_nodes     = {graph.number_of_nodes()}")
    print(f"total_edges     = {graph.number_of_edges()}")
    print(f"expected_depth  = {expected_depth}")
    print(f"actual_depth    = {actual_depth}")
    print(f"calculate_basis = {calculate_basis}")
    print("-" * 72)

    start_time = time.perf_counter()
    lp, betti, basis = max_path_homology(graph, calculate_basis=calculate_basis, report_sparsity=True)
    end_time = time.perf_counter()
    elapsed_time = end_time - start_time

    print("-" * 72)
    print("Results:")
    print(f"lp              = {lp}")
    print(f"betti           = {betti}")
    print(f"elapsed_seconds = {elapsed_time:.6f}")
    print("=" * 72)
    print("Example 2")
    print("-" * 72)
    edgelist = [('a0', 'b2'), ('a0', 'b3'), ('a1', 'b2'), ('a1', 'b3'), 
                ('a1', 'c1'), ('b0', 'c0'), ('b0', 'c1'), ('b1', 'c0'), 
                ('b1', 'c1'), ('b2', 'c2'), ('b2', 'c3'), ('b3', 'c2'), 
                ('b3', 'c3'), ('b4', 'c4'), ('b4', 'c5'), ('b5', 'c4'), 
                ('b5', 'c5'), ('b0', 'c2'), ('b0', 'c3'), ('b1', 'c2'), 
                ('b1', 'c3'), ('b4', 'c2'), ('b4', 'c3'), ('b5', 'c2'), 
                ('b5', 'c3'), ('c0', 'd0'), ('c0', 'd1'), ('c1', 'd0'), 
                ('c1', 'd1'), ('c4', 'd2'), ('c4', 'd3'), ('c5', 'd2'), 
                ('c5', 'd3')]
    G = nx.DiGraph(edgelist)
    lp, betti, basis = max_path_homology(G, calculate_basis = True, report_sparsity = True)
    print(f"lp              = {lp}")
    print(f"betti           = {betti}")
    print(f"basis           = {basis}")