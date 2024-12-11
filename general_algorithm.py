# The following functions are based on an implementation from the following open-source project:
# Carranza, D., Doherty, B., Kapulkin, K., Opie, M., Sarazola, M., & Wong, L. Z. (2022). Python script for computing path homology of digraphs (Version 1.0.0) [Computer software]
# https://github.com/sheaves/path_homology, under the GNU General Public License.

import networkx as nx
import numpy as np
from sympy import Matrix
from maxph_matrix import rref_null

def R_path(G, cutoff):
    """
    Generate regular allowed paths in G, up to length cutoff.
    """
    R = []
    R.append([(v,) for v in G.nodes()])
    R.append([e for e in G.edges()])

    for n in range(2, cutoff + 1):
        R_n = []
        for path in R[n-1]:
            last_node = path[-1]
            for next_node in G.successors(last_node):
                R_n.append(path + (next_node,))
        R.append(R_n)
    return R

def d(path):
    """
    Compute the differential of a path, assuming all paths are allowed.
    Returns a dictionary mapping subpaths to coefficients.
    """
    coeffs = {}
    for i in range(len(path)):
        head = path[:i]
        tail = path[i+1:]
        if (len(head) == 0) or (len(tail) == 0) or head[-1] != tail[0]: 
            subpath = head + tail
            coeffs[subpath] = (-1)**i
    return coeffs

def D_full(R, n):
    """
    Compute the differentials of all paths of length n.
    Returns the differential matrix, and the row and column indices.
    """
    if len(R[n]) == 0:
        return None, [], []

    paths_n = sorted(R[n])
    subpaths = set()

    for path in paths_n:
        subpaths.update(d(path).keys())

    paths_n_minus1 = sorted(subpaths)
    row_idx = {path: idx for idx, path in enumerate(paths_n_minus1)}
    col_idx = {path: idx for idx, path in enumerate(paths_n)}
    D_matrix = [[0 for _ in paths_n] for _ in paths_n_minus1]

    for col_path in paths_n:
        d_col = d(col_path)
        for subpath, coeff in d_col.items():
            if subpath in row_idx:
                i = row_idx[subpath]
                j = col_idx[col_path]
                D_matrix[i][j] = coeff

    return D_matrix, paths_n_minus1, paths_n

def O_n(R):
    """
    Generate Omega chain complex for path homology 
    from list of list of paths, R

    O[n] contains generators for Omega_n,
    expressed as lists of linear combinations of paths in R[n]
    
    Linear combinations are expressed as dictionaries:
        - keys are elements of R[n]
        - values are the coefficients in the linear combination
    """
    Omega = []
    # O_0
    Omega.append([{(v,): 1} for v in R[0]])
    # O_1
    Omega.append([{e: 1} for e in R[1]])
    # O_n
    for n in range(2, len(R)):
        D_matrix, row_paths, col_paths = D_full(R, n)
        if D_matrix:
            D_sympy = Matrix(D_matrix)
            row_indices = {path: idx for idx, path in enumerate(row_paths)}
            rows_to_remove = [row_indices[path] for path in R[n-1] if path in row_indices]
            D_sympy_reduced = D_sympy

            for idx in sorted(rows_to_remove, reverse=True):
                D_sympy_reduced.row_del(idx)

            basis = D_sympy_reduced.nullspace()
            O_n_basis = []

            for vec in basis:
                vec = np.array(vec).astype(int).flatten()
                coeffs = {col_paths[i]: vec[i] for i in range(len(vec)) if vec[i] != 0}
                O_n_basis.append(coeffs)
            Omega.append(O_n_basis)
        else:
            Omega.append([{p: 1} for p in R[n]])
    return Omega

def D(R, Omega, n):
    """
    Generate the n-th differential matrix in the Omega chain complex, from O[n] to O[n-1].
    Returns the differential matrix (SymPy matrix).
    """
    if n == 0:
        return None

    D_matrix, _, col_paths = D_full(R, n)
    if D_matrix is None:
        return None

    D_sympy = Matrix(D_matrix)
    O_n_indices = []

    for generator in Omega[n]:
        vec = [0] * len(col_paths)
        for path, coeff in generator.items():
            idx = col_paths.index(path)
            vec[idx] = coeff
        O_n_indices.append(vec)

    O_n_matrix = Matrix(O_n_indices).T
    if O_n_matrix.cols == 0:
        return None
    D_result = D_sympy * O_n_matrix

    return D_result

def H_path_R(R):
    """
    Compute the reduced path homology of a path set R.    
    Returns:
        H : dimensions of path homology
        C : dimensions of Omega chain complex
        Diffs : differentials
        R : regular allowed paths
        Omega : generators of Omega 
    """
    Omega = O_n(R)
    cutoff = len(R) - 1
    Diffs = [D(R, Omega, i) for i in range(cutoff + 1)]
    C = [len(Omega[i]) for i in range(cutoff + 1)]
    H = []

    # H_0
    if Diffs[1] is not None:
        ker_dim = len(Omega[0]) - 1 #reduced path homology
        D_img = np.array(Diffs[1].tolist(), dtype = int)
        img_dim = D_img.shape[1] - rref_null(D_img).shape[1]
        H.append(ker_dim - img_dim)
    else:
        H.append(max(0, len(Omega[0]) - 1))

    # H_n，n > 0
    for n in range(1, cutoff):
        if Diffs[n] is not None:
            D_ker = np.array(Diffs[n].tolist(), dtype = int)
            ker_dim = rref_null(D_ker).shape[1]
            if Diffs[n + 1] is not None:
                D_img = np.array(Diffs[n + 1].tolist(), dtype = int)
                img_dim = D_img.shape[1] - rref_null(D_img).shape[1]
            else:
                img_dim = 0
            H.append(ker_dim - img_dim)
        else:
            H.append(0)
    return H, C, Diffs, R, Omega

def edgelist_to_graph(edgelist, sep = ":", vertices = None):
    """
    Generates a networkx DiGraph from a list of directed edges.
    Self-loops and parallel edges will be ignored.
    
    Examples of edgelists:
        - ['ab', 'ac', 'bd', 'bc'] (sep = '', i.e. no separator)
        - ['10:11', '10:12', '11:13', '12:13']  (sep = ":")
    
    sep is required if some vertices have length > 1
    
    Vertices will be inferred from the edgelist, but you may specify them:
        - have more control over the order of vertices
        - have "floating" vertices that don't belong to any edge
    
    Examples of vertices:
        - ['a', 'b', 'c', 'd']
        - ['10', '11', '12', '13']
        
    """
    G = nx.DiGraph()
    # Add vertices, if specified
    if vertices is not None:
        G.add_nodes_from(vertices)
    
    # Add edges
    for e in edgelist:
        e0, e1 = e.split(sep)

        # Check for self-loops. Only add non-self-loops
        if e0 != e1:
            G.add_edge(e0, e1)
        else:
            G.add_node(e0)

    return G