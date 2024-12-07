import numpy as np
import sympy as sp

def rref_null(matrix: np.ndarray[np.int8]) -> np.ndarray[np.int8]:
    """
    This function calculates the null space of a given matrix by first converting it to its 
    Reduced Row Echelon Form (RREF). The input matrix consists of elements in {-1, 0, 1},
    and three matrix row operations are used: exchange any two rows, multiply a row by -1, and
    add one row to another to cancel the nonzero entries in the pivot column.

    Args:
        matrix (np.ndarray[np.int8]): The input matrix for which the null space is to be computed.

    Returns:
        null_matrix (np.ndarray[np.int8]): The null space of the input matrix.
    """
    if matrix.size == 0:
        return np.array([], dtype = np.int8)

    rref = matrix.copy()
    rows, cols = rref.shape
    r, pivot_col_indices = 0, {}
    for c in range(cols):
        if r >= rows:
            break        
        pivot = np.nonzero(rref[r:, c])[0]
        if len(pivot) == 0:
            continue
        pivot_1st = pivot[0] + r
        if pivot_1st != r:
            rref[[r, pivot_1st]] = rref[[pivot_1st, r]] 
        if rref[r, c] < 0:
            rref[r] *= -1
        for i in range(r + 1, rows):
            if rref[i, c] != 0:
                rref[i] -= rref[i, c] * rref[r]  
        pivot_col_indices[c] = r
        r += 1

    for c, r in reversed(pivot_col_indices.items()):
        for i in range(r):
            if rref[i, c] != 0:
                rref[i] -= rref[i, c] * rref[r]

    null_basis = []
    for c in range(cols):
        if c not in pivot_col_indices:
            vector = np.zeros(cols, dtype = np.int8)
            vector[c] = 1
        else:
            vector = -rref[pivot_col_indices[c]] 
        null_basis.append(vector)

    null_space = np.array(null_basis, dtype = np.int8)
    null_matrix = np.delete(null_space, list(pivot_col_indices.keys()), axis = 1) 

    return null_matrix

def row_del(matrix: np.ndarray[np.int8], partition_rule: dict[str, int], keys: list[str]) -> np.ndarray[np.int8]:
    """
    Deletes specific rows from the input matrix according to a given partition rule and a list of keys.

    Args:
        matrix (np.ndarray[np.int8]):    The input matrix.
        partition_rule (dict[str, int]): A dictionary defining the partitioning of the rows, where keys 
                                         are strings and values are integers specifying the number of rows to delete.
        keys (list[str]):                A list of keys specifying which rows to delete from the matrix.

    Returns:
        reduced_matrix (np.ndarray[np.int8]): A reduced matrix with the specified rows removed.
    """
    start_index, rows_to_delete = 0, []

    for key, count in partition_rule.items():
        if count >=1:
            end_index = start_index + count
        else:
            continue
        if key in keys:
            rows_to_delete.extend(range(start_index, end_index))
        start_index = end_index

    return np.delete(matrix, rows_to_delete, axis=0) if rows_to_delete else matrix

def basis_mat(A_x_list: np.ndarray[np.int8], partition: dict[str, int], null_matrix: np.ndarray[np.int8]) -> sp.Matrix:
    """
    Generate a symbolic matrix for updating the basis of the null space.

    Args:
        A_x_list (np.ndarray[np.int8]):    The input matrix.
        partition (dict[str, int]):        A dictionary where the keys are the column indices and
                                           the values are the number of columns to replace.
        null_matrix (np.ndarray[np.int8]): The null space of the input matrix.

    Returns:
        sym_matrix (sp.Matrix): The symbolic matrix.
    """
    symbol_matrix = sp.zeros(A_x_list[0].shape[0], null_matrix.shape[1])
    current_row_idx = 0
    A_idx = 0
    for key, count in partition.items():
        if count != 0:
            mat_mul = A_x_list[A_idx] @ null_matrix[current_row_idx: current_row_idx + count, :] 
            symbol_matrix += sp.Symbol(key) * sp.Matrix(mat_mul)
            current_row_idx += count
            A_idx += 1

    return symbol_matrix

def null_identity(partition: dict[str, int], index: int, M: np.ndarray[np.int8], dim: int) -> tuple[np.ndarray[np.int8], int]:
    """
    Construct the null space of a block matrix M when a specific block is an identity matrix.

    Args:
        partition (dict[str, int]): A dictionary that partitions the matrix into blocks, where each key 
                                    corresponds to a block matrix and the associated value represents the column number.
        index (int):                The index of the first identity matrix in M.
        M (np.ndarray[np.int8]):    The matrix for which the null space is being computed.
        dim (int):                  The dimension of the null space of matrix M.

    Returns:
        tuple[np.array[int], int]:
            block (np.ndarray[np.int8]): A block diagonal matrix representing the null space of M.
            total_size - dim (int):           The column number of the block_diag matrix.
    """
    total_size = sum(partition.values())
    block = np.zeros((total_size, total_size - dim), dtype = np.int8)
    current_row, current_col = 0, 0

    for idx, (_, count) in enumerate(partition.items()):
        if idx == index:
            if M.shape != (count, total_size):
                raise ValueError("M's size must match the partition.")
            block[current_row:current_row + count, :current_col] = M[:, :current_col]
            block[current_row:current_row + count, current_col:] = M[:, current_col + dim:]
        else:
            block[current_row:current_row + count, current_col:current_col + count] = np.identity(count, dtype = np.int8)
            current_col += count
        current_row += count
    return block, total_size - dim