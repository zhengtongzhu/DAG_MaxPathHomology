import numpy as np
import networkx as nx
import sympy as sp

from maxph_matrix import rref_null, row_del, node_symbol, basis_mat, null_identity
from dag_preprocess import dag_process

def max_path_homology(G: nx.DiGraph, calculate_basis: bool = False) -> tuple[int, int, list[list[sp.Expr]] | None]:
    """
    Compute the maximal path homology of a finite DAG.

    Args:
        G: 
            A finite NetworkX directed acyclic graph.
        calculate_basis: 
            Whether to compute an explicit homology basis.

    Returns:
        lp:
            The maximum directed path length of G.
        betti:
            The lp-dimensional Betti number.
        basis:
            A basis if requested; otherwise None.
    """

    (subgraph_dict, N, lp, num_subgraphs, graph_list) = dag_process(G)
    betti = 0
    basis = [] if calculate_basis else None

    if lp == 0:
        iso_nodes = list(G.nodes())
        betti = len(iso_nodes) - 1

        if not calculate_basis:
            return lp, betti, basis
        else:
            if len(iso_nodes) <= 1:
                return lp, betti, basis
            else:
                symbols = {node: node_symbol(node) for node in iso_nodes}
                first_node = iso_nodes[0]
                basis = [[symbols[node] - symbols[first_node] for node in iso_nodes[1:]]]
                return lp, betti, basis
    elif num_subgraphs == 0:
        return lp, betti, basis

    assert len(N) == num_subgraphs

    for idx_subgraphs in range(num_subgraphs):
        if any(num == 1 for num in N[idx_subgraphs]):
            continue

        dim = N[idx_subgraphs][0] - 1
        V = np.vstack((-np.ones((1, dim), dtype = np.int8), np.identity(dim, dtype = np.int8)))    
        partition = {node: 1 for node in subgraph_dict[idx_subgraphs][0]}

        if calculate_basis:
            keys = list(partition.keys())
            symbols = {key: node_symbol(key) for key in keys}
            basis_iteration = sp.Matrix([symbols[key]- symbols[keys[0]] for key in keys[1:]]).T

        for layer_idx in range(1, lp + 1):  
            col_list, P, Status = [], {}, False
            V_rows = V.shape[0]
            for idx, node in enumerate(subgraph_dict[idx_subgraphs][layer_idx]):
                predecessors = list(graph_list[idx_subgraphs].predecessors(node))
                sum_pred = sum(partition.get(pred, 0) for pred in predecessors)
                
                if sum_pred == V_rows:
                    P[node] = dim
                    col_list.append(-np.identity(dim, dtype = np.int8))
                    if not Status:
                        Status, index = True, idx
                else:
                    Reduced_matrix = row_del(V, partition, predecessors)
                    A_x = rref_null(Reduced_matrix)

                    if A_x.size == 0:
                        P[node] = 0
                        continue

                    dim_A_x = A_x.shape[1]
                    P[node] = dim_A_x
                    if dim_A_x == dim:
                        col_list.append(-np.identity(dim, dtype = np.int8))
                        if not Status:
                            Status, index = True, idx
                    else:
                        col_list.append(-A_x)

            if not col_list:
                dim = 0
                break
            
            A = np.hstack(col_list)

            if Status:
                V, dim = null_identity(P, index, A, dim)                 
            else:
                V = rref_null(A)
                dim = V.shape[1]

            if dim == 0:
                break
            if calculate_basis:
                basis_iteration = basis_iteration * basis_mat(col_list, P, V)

            partition = P
            
        betti += dim

        if calculate_basis and dim != 0:
            basis.extend(basis_iteration.tolist())

    return lp, betti, basis

if __name__ == "__main__":

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
    lp, betti, basis = max_path_homology(G, calculate_basis = True)
    print(f"lp = {lp},\nbetti = {betti},\nbasis = {basis}\n")