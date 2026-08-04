from __future__ import annotations

import numpy as np
import sympy as sp

from collections.abc import Mapping, Sequence
from scipy import sparse as sps

np.set_printoptions(threshold=np.inf, linewidth=np.inf, precision=6, suppress=False)
SparseMatrix = sps.csr_matrix

def node_symbol(node) -> sp.Symbol:
    return sp.Symbol(str(node), commutative=False)

def to_csr_sparse(matrix, *, dtype=np.int32, copy: bool = False) -> SparseMatrix:
    """
    Convert a dense or sparse matrix to canonical CSR format.

    Duplicate entries are summed, explicit zeros are removed, and column
    indices are sorted.
    """
    if sps.issparse(matrix):
        result = matrix.tocsr(copy=copy).astype(dtype, copy=False)
    else:
        result = sps.csr_matrix(np.asarray(matrix, dtype=dtype), copy=copy)

    result.sum_duplicates()
    result.eliminate_zeros()
    result.sort_indices()
    return result

def sparse_density(matrix) -> float:
    """Return nnz / (number of rows * number of columns)."""
    matrix_csr = to_csr_sparse(matrix)
    total_entries = matrix_csr.shape[0] * matrix_csr.shape[1]

    if total_entries == 0:
        return 0.0

    return float(matrix_csr.nnz / total_entries)

def sparse_memory_mib(matrix) -> float:
    """Return the storage used by the CSR data structures in MiB."""
    matrix_csr = to_csr_sparse(matrix)
    total_bytes = (matrix_csr.data.nbytes + matrix_csr.indices.nbytes + matrix_csr.indptr.nbytes)

    return float(total_bytes / (1024**2))

def _checked_int32_array(values: Sequence[int]) -> np.ndarray:
    """Convert Python integers to int32"""
    lower = np.iinfo(np.int32).min
    upper = np.iinfo(np.int32).max

    for value in values:
        if value < lower or value > upper:
            raise OverflowError("An elimination coefficient does not fit in int32. Use an exact-arithmetic backend for this matrix.")

    return np.asarray(values, dtype=np.int32)

def _swap_sparse_rows(rows: list[dict[int, int]], column_rows: list[set[int]], first: int, second: int) -> None:
    """Swap two dictionary rows and update the column-to-row index."""
    if first == second:
        return

    first_columns = set(rows[first])
    second_columns = set(rows[second])

    for column in first_columns - second_columns:
        column_rows[column].discard(first)
        column_rows[column].add(second)

    for column in second_columns - first_columns:
        column_rows[column].discard(second)
        column_rows[column].add(first)

    rows[first], rows[second] = rows[second], rows[first]

def _add_scaled_sparse_row(rows: list[dict[int, int]], column_rows: list[set[int]], target_index: int, source_row: dict[int, int], scale: int) -> None:
    """
    Perform target <- target + scale * source using Python integers.

    The column-to-row index is updated whenever a nonzero entry is created
    or removed.
    """
    if scale == 0:
        return

    target_row = rows[target_index]

    for column, source_value in source_row.items():
        old_value = target_row.get(column, 0)
        new_value = old_value + scale * source_value

        if new_value == 0:
            if old_value != 0:
                target_row.pop(column, None)
                column_rows[column].discard(target_index)
        else:
            target_row[column] = new_value

            if old_value == 0:
                column_rows[column].add(target_index)

def rref_null_sparse(matrix) -> SparseMatrix:
    """
    Compute a basis of the null space using sparse row elimination.

    Parameters
    ----------
    matrix:
        A dense or sparse integer matrix of shape (rows, columns).

    Returns
    -------
    scipy.sparse.csr_matrix
        A matrix of shape (columns, nullity). Its columns form a basis of
        the null space.

    Markdown
    ---------
    At every pivot position, a pivot equal to +1 or -1 must be available.
    If a nonzero column has no unit pivot, an ArithmeticError is raised
    instead of silently returning an incorrect result.

    The elimination itself uses Python integers, so intermediate row
    operations do not overflow. The returned sparse matrix is stored as
    int32 and an OverflowError is raised if a coefficient does not fit.
    """
    matrix_csr = to_csr_sparse(matrix, dtype=np.int32)
    row_count, column_count = matrix_csr.shape

    if column_count == 0:
        return sps.csr_matrix((0, 0), dtype=np.int32)

    if row_count == 0 or matrix_csr.nnz == 0:
        return sps.eye(column_count, dtype=np.int32, format="csr")

    sparse_rows: list[dict[int, int]] = []
    column_rows: list[set[int]] = [set() for _ in range(column_count)]

    for row_index in range(row_count):
        start = matrix_csr.indptr[row_index]
        end = matrix_csr.indptr[row_index + 1]

        row_map: dict[int, int] = {}

        for column, value in zip(matrix_csr.indices[start:end], matrix_csr.data[start:end],):
            integer_value = int(value)

            if integer_value == 0:
                continue

            integer_column = int(column)
            row_map[integer_column] = integer_value
            column_rows[integer_column].add(row_index)

        sparse_rows.append(row_map)

    pivot_column_to_row: dict[int, int] = {}
    current_pivot_row = 0

    # Forward elimination.
    for column in range(column_count):
        if current_pivot_row >= row_count:
            break

        candidate_rows = [row_index for row_index in column_rows[column] if row_index >= current_pivot_row]

        if not candidate_rows:
            continue

        unit_candidates = [row_index for row_index in candidate_rows if sparse_rows[row_index].get(column, 0) in (-1, 1)]

        if not unit_candidates:
            encountered = sorted(
                {sparse_rows[row_index].get(column, 0) for row_index in candidate_rows}
            )
            raise ArithmeticError(
                "A nonzero pivot column has no +1 or -1 pivot. "
                "The original elimination rule is not valid for this case. "
                f"Column: {column}; encountered values: {encountered[:10]}"
            )

        selected_row = min(unit_candidates)

        _swap_sparse_rows(sparse_rows, column_rows, current_pivot_row, selected_row)

        pivot_value = sparse_rows[current_pivot_row][column]

        if pivot_value == -1:
            sparse_rows[current_pivot_row] = {index: -value for index, value in sparse_rows[current_pivot_row].items()}

        pivot_row = sparse_rows[current_pivot_row]

        rows_to_eliminate = sorted(row_index for row_index in column_rows[column] if row_index > current_pivot_row)

        for row_index in rows_to_eliminate:
            factor = sparse_rows[row_index].get(column, 0)

            if factor != 0:
                _add_scaled_sparse_row(sparse_rows, column_rows, row_index, pivot_row, -factor)

        pivot_column_to_row[column] = current_pivot_row
        current_pivot_row += 1

    # Backward elimination.
    for column, pivot_row_index in reversed(tuple(pivot_column_to_row.items())):
        pivot_row = sparse_rows[pivot_row_index]
        rows_to_eliminate = sorted(row_index for row_index in column_rows[column] if row_index < pivot_row_index)

        for row_index in rows_to_eliminate:
            factor = sparse_rows[row_index].get(column, 0)

            if factor != 0:
                _add_scaled_sparse_row(sparse_rows, column_rows, row_index, pivot_row, -factor,)

    pivot_columns = set(pivot_column_to_row)
    free_columns = [column for column in range(column_count) if column not in pivot_columns]

    nullity = len(free_columns)

    if nullity == 0:
        return sps.csr_matrix((column_count, 0), dtype=np.int32)

    free_column_position = {column: basis_index for basis_index, column in enumerate(free_columns)}

    output_rows: list[int] = []
    output_columns: list[int] = []
    output_values: list[int] = []

    # Identity entries associated with free variables.
    for basis_index, free_column in enumerate(free_columns):
        output_rows.append(free_column)
        output_columns.append(basis_index)
        output_values.append(1)

    # Pivot-variable coefficients.
    for pivot_column, pivot_row_index in pivot_column_to_row.items():
        pivot_row = sparse_rows[pivot_row_index]

        for column, value in pivot_row.items():
            basis_index = free_column_position.get(column)

            if basis_index is None:
                continue

            coefficient = -value

            if coefficient != 0:
                output_rows.append(pivot_column)
                output_columns.append(basis_index)
                output_values.append(coefficient)

    null_matrix = sps.coo_matrix(
        (
            _checked_int32_array(output_values),
            (
                np.asarray(output_rows, dtype=np.int32),
                np.asarray(output_columns, dtype=np.int32),
            ),
        ),
        shape=(column_count, nullity),
        dtype=np.int32,
    ).tocsr()

    null_matrix.sum_duplicates()
    null_matrix.eliminate_zeros()
    null_matrix.sort_indices()

    return null_matrix

def row_del_sparse(matrix, partition_rule: Mapping[str, int], keys: Sequence[str]) -> SparseMatrix:
    """
    Delete row blocks specified by keys and return a CSR matrix.
    """
    matrix_csr = to_csr_sparse(matrix)
    keys_to_delete = set(keys)

    keep_rows = np.ones(matrix_csr.shape[0], dtype=bool)

    start_index = 0

    for key, count in partition_rule.items():
        if count < 0:
            raise ValueError("Partition counts must be nonnegative.")

        end_index = start_index + count

        if key in keys_to_delete:
            keep_rows[start_index:end_index] = False

        start_index = end_index

    if start_index != matrix_csr.shape[0]:
        raise ValueError("The partition sizes do not match the number of matrix rows.")

    if np.all(keep_rows):
        return matrix_csr

    result = matrix_csr[keep_rows, :].tocsr()
    result.eliminate_zeros()
    result.sort_indices()
    return result

def basis_mat_sparse(A_x_list: Sequence, partition: Mapping[str, int], null_matrix) -> sp.SparseMatrix:
    """
    Generate the symbolic basis-update matrix using sparse arithmetic.

    The numerical products are computed with SciPy sparse matrices. The
    symbolic result is stored as a SymPy SparseMatrix.

    Parameters
    ----------
    A_x_list:
        Sequence of numerical matrices, one for each nonempty partition block.

    partition:
        Ordered mapping from node labels to block widths.

    null_matrix:
        Sparse null-space basis whose rows follow the partition.

    Returns
    -------
    sympy.SparseMatrix
        The symbolic matrix without explicitly storing zero entries.
    """
    if len(A_x_list) == 0:
        raise ValueError("A_x_list must not be empty.")

    null_csr = to_csr_sparse(null_matrix)

    if sum(partition.values()) != null_csr.shape[0]:
        raise ValueError("The partition sizes do not match the rows of null_matrix.")

    active_partition = [(key, count) for key, count in partition.items() if count != 0]

    if len(active_partition) != len(A_x_list):
        raise ValueError(
            "The number of nonempty partition blocks must match "
            "the number of matrices in A_x_list."
        )

    output_row_count = A_x_list[0].shape[0]
    output_column_count = null_csr.shape[1]

    symbolic_entries: dict[tuple[int, int], sp.Expr] = {}

    current_row_index = 0

    for A_x, (key, count) in zip(A_x_list, active_partition):
        A_csr = to_csr_sparse(A_x)

        if A_csr.shape[0] != output_row_count:
            raise ValueError("All matrices in A_x_list must have the same row count.")

        if A_csr.shape[1] != count:
            raise ValueError(
                f"The A_x block for node {key!r} has "
                f"{A_csr.shape[1]} columns, but the partition specifies "
                f"{count}."
            )

        null_block = null_csr[current_row_index: current_row_index + count, :,]

        product = (A_csr @ null_block).tocoo()
        product.sum_duplicates()

        symbol = node_symbol(key)

        for row, column, value in zip(product.row, product.col, product.data):
            integer_value = int(value)

            if integer_value == 0:
                continue

            location = (int(row), int(column))
            term = symbol * sp.Integer(integer_value)
            updated_expression = (symbolic_entries.get(location, sp.S.Zero) + term)

            if updated_expression == 0:
                symbolic_entries.pop(location, None)
            else:
                symbolic_entries[location] = updated_expression

        current_row_index += count

    return sp.SparseMatrix(output_row_count, output_column_count, symbolic_entries)

def null_identity_sparse(partition: Mapping[str, int], index: int, M, dim: int) -> tuple[SparseMatrix, int]:
    """
    Construct the  nullspace basis.

    Parameters
    ----------
    partition:
        Ordered mapping from block labels to block widths.

    index:
        Position of the selected identity block.

    M:
        Matrix of shape (dim, total_size).

    dim:
        Width of the selected identity block.

    Returns
    -------
    (block, nullity):
        block is a CSR matrix of shape
        (total_size, total_size - dim).
    """
    counts = list(partition.values())

    if any(count < 0 for count in counts):
        raise ValueError("Partition counts must be nonnegative.")

    if not 0 <= index < len(counts):
        raise IndexError("index is outside the partition.")

    if counts[index] != dim:
        raise ValueError("The selected block size must equal dim.")

    offsets = np.zeros(len(counts) + 1, dtype=np.int32,)
    offsets[1:] = np.cumsum(counts)

    total_size = int(offsets[-1])
    nullity = total_size - dim

    pivot_start = int(offsets[index])
    pivot_end = int(offsets[index + 1])

    M_csr = to_csr_sparse(M)

    if M_csr.shape != (dim, total_size):
        raise ValueError(f"M must have shape {(dim, total_size)}, but received {M_csr.shape}.")

    if nullity == 0:
        return (sps.csr_matrix((total_size, 0), dtype=np.int32,), 0)

    free_source_columns = np.concatenate(
        (
            np.arange(0, pivot_start, dtype=np.int32),
            np.arange(pivot_end, total_size, dtype=np.int32),
        )
    )

    # Identity entries in all non-pivot coordinate blocks.
    row_parts: list[np.ndarray] = [free_source_columns]
    column_parts: list[np.ndarray] = [np.arange(nullity, dtype=np.int32)]
    value_parts: list[np.ndarray] = [np.ones(nullity, dtype=np.int32)]

    # Remaining columns of M fill the selected block row.
    M_free = M_csr[:, free_source_columns,].tocoo()

    if M_free.nnz:
        row_parts.append(pivot_start + M_free.row.astype(np.int32, copy=False))
        column_parts.append(M_free.col.astype(np.int32, copy=False))
        value_parts.append(M_free.data.astype(np.int32, copy=False))

    block = sps.coo_matrix(
        (
            np.concatenate(value_parts),
            (
                np.concatenate(row_parts),
                np.concatenate(column_parts),
            ),
        ),
        shape=(total_size, nullity),
        dtype=np.int32,
    ).tocsr()

    block.sum_duplicates()
    block.eliminate_zeros()
    block.sort_indices()

    return block, nullity

def hstack_sparse(blocks: Sequence) -> SparseMatrix:
    """Sparse replacement for np.hstack."""
    if len(blocks) == 0:
        return sps.csr_matrix((0, 0), dtype=np.int32)

    return sps.hstack([to_csr_sparse(block) for block in blocks], format="csr", dtype=np.int32)

def vstack_sparse(blocks: Sequence) -> SparseMatrix:
    """Sparse replacement for np.vstack."""
    if len(blocks) == 0:
        return sps.csr_matrix((0, 0), dtype=np.int32)

    return sps.vstack(
        [to_csr_sparse(block) for block in blocks],
        format="csr",
        dtype=np.int32,
    )

def identity_sparse(size: int) -> SparseMatrix:
    """Sparse replacement for np.identity(size)."""
    return sps.eye(size, dtype=np.int32, format="csr")

def zeros_sparse(row_count: int, column_count: int) -> SparseMatrix:
    """Sparse replacement for np.zeros((row_count, column_count))."""
    if row_count < 0 or column_count < 0:
        raise ValueError("Matrix dimensions must be nonnegative.")

    return sps.csr_matrix((row_count, column_count), dtype=np.int32)