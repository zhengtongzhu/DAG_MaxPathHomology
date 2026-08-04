from __future__ import annotations

import csv
import gc
import random
import time
import traceback
import tracemalloc
import matplotlib
matplotlib.use("Agg", force = True)

import matplotlib.pyplot as plt
import multiprocessing as mp
import networkx as nx
import numpy as np

from tqdm import trange
from typing import Any
from recursive_algorithm import max_path_homology
from general_algorithm import R_path, H_path_R, edgelist_to_graph


MIB = 1024 ** 2

VALID_ALGORITHMS = {"general", "recursive_basis", "recursive_no_basis",}

def build_general_graph(nodes: list[str], edges: list[tuple[str, str]],):
    """
    Build the general implementation's graph while preserving isolated nodes.

    The edge-list converter may omit vertices that are not incident to an
    edge. Explicitly restoring them ensures that all implementations receive
    the same mathematical graph.
    """
    general_edgelist = [f"{source}:{target}" for source, target in edges]
    graph = edgelist_to_graph(general_edgelist)

    if not hasattr(graph, "nodes"):
        raise TypeError("The graph returned by edgelist_to_graph does not expose nodes.")

    nodes_attribute = graph.nodes
    existing_nodes = set(nodes_attribute() if callable(nodes_attribute) else nodes_attribute)

    missing_nodes = [node for node in nodes if node not in existing_nodes]

    if missing_nodes:
        if hasattr(graph, "add_nodes_from"):
            graph.add_nodes_from(missing_nodes)
        elif hasattr(graph, "add_node"):
            for node in missing_nodes:
                graph.add_node(node)
        else:
            raise TypeError(
                "The general graph omits isolated nodes and does not support "
                "add_nodes_from() or add_node()."
            )

    return graph

def run_general_algorithm(G_com: nx.DiGraph, max_L: int,) -> list[int]: 
    """
    Run the general path-homology algorithm.
    """
    R = R_path(G_com,cutoff=max_L + 1,)
    H, _, _, _, _ = H_path_R(R)

    return H

def _tracemalloc_worker(
    algorithm: str,
    nodes: list[str],
    edges: list[tuple[str, str]],
    max_depth: int,
    connection: Any,
) -> None:
    """
    Measure one algorithm call with tracemalloc in a fresh process.

    The child process imports this module before entering this worker. The
    algorithm-specific input is also constructed before tracemalloc starts.
    Consequently, Python startup, package imports, argument deserialization,
    and input construction are excluded from the reported traced allocation
    peak.

    No warm-up call is performed. The measurement therefore represents the
    first formal call on the target graph after all imports and input
    construction have completed.
    """
    result: Any = None
    algorithm_input: Any = None

    try:
        if algorithm not in VALID_ALGORITHMS:
            raise ValueError(f"Unknown algorithm: {algorithm!r}.")

        # Construct the input before tracing begins.
        if algorithm == "general":
            algorithm_input = build_general_graph(nodes, edges,)
        else:
            algorithm_input = nx.DiGraph()
            algorithm_input.add_nodes_from(nodes)
            algorithm_input.add_edges_from(edges)

        del nodes
        del edges
        gc.collect()

        if tracemalloc.is_tracing():
            tracemalloc.stop()

        tracemalloc.start(1)
        tracemalloc.reset_peak()

        start_time = time.perf_counter()

        if algorithm == "general":
            result = run_general_algorithm(algorithm_input, max_depth,)
            betti_number = int(result[max_depth])
        else:
            result = max_path_homology(
                algorithm_input,
                calculate_basis=(algorithm == "recursive_basis"),
                report_sparsity=False,
            )
            betti_number = int(result[1])

        runtime_seconds = time.perf_counter() - start_time

        # Read the values while the returned result remains alive. Therefore,
        # requested output objects, including a symbolic basis, are included.
        current_bytes, peak_bytes = (tracemalloc.get_traced_memory())
        tracer_overhead_bytes = (tracemalloc.get_tracemalloc_memory())
        tracemalloc.stop()

        connection.send(
            {
                "ok": True,
                "algorithm": algorithm,
                "betti_number": betti_number,
                "tracemalloc_current_bytes": int(current_bytes),
                "tracemalloc_peak_bytes": int(peak_bytes),
                "tracemalloc_overhead_bytes": int(
                    tracer_overhead_bytes
                ),
                "runtime_seconds": float(runtime_seconds),
            }
        )

    except BaseException:
        if tracemalloc.is_tracing():
            tracemalloc.stop()

        try:
            connection.send(
                {
                    "ok": False,
                    "algorithm": algorithm,
                    "error": traceback.format_exc(),
                }
            )
        except (BrokenPipeError, EOFError, OSError):
            pass

    finally:
        result = None
        algorithm_input = None
        connection.close()

def _remaining_timeout(deadline: float | None) -> float | None:
    if deadline is None:
        return None

    return max(0.0, deadline - time.monotonic())

def _wait_for_worker_result(
    connection: Any,
    process: Any,
    deadline: float | None,
    algorithm: str,
) -> dict[str, Any]:
    """Wait for a worker result while detecting timeout or early exit."""
    while True:
        remaining = _remaining_timeout(deadline)

        if remaining is not None and remaining <= 0:
            raise TimeoutError(
                "Tracemalloc measurement timed out. "
                f"algorithm={algorithm!r}."
            )

        poll_timeout = (
            0.05
            if remaining is None
            else min(0.05, remaining)
        )

        if connection.poll(poll_timeout):
            try:
                message = connection.recv()
            except EOFError as error:
                raise RuntimeError(
                    "The tracemalloc worker closed its connection "
                    "without returning data. "
                    f"algorithm={algorithm!r}."
                ) from error

            if not isinstance(message, dict):
                raise RuntimeError(
                    "The tracemalloc worker returned an invalid "
                    f"message: {message!r}."
                )

            return message

        if not process.is_alive():
            process.join()

            if connection.poll():
                try:
                    message = connection.recv()
                except EOFError:
                    message = None

                if isinstance(message, dict):
                    return message

            raise RuntimeError(
                "The tracemalloc worker exited without returning data. "
                f"algorithm={algorithm!r}, "
                f"exit_code={process.exitcode}."
            )

def measure_algorithm_traced_memory(
    algorithm: str,
    nodes: list[str],
    edges: list[tuple[str, str]],
    max_depth: int,
    timeout: float | None = None,
) -> dict[str, Any]:
    """
    Measure peak traced allocations for one algorithm in isolation.

    A fresh spawned process is used for each run. Tracemalloc starts only
    after all imports and implementation-specific input construction have
    completed. No warm-up is performed.

    The primary value, ``tracemalloc_peak_mib``, is the maximum total size of
    memory blocks tracked by tracemalloc during the target algorithm call. It
    is not process RSS and may omit native allocations made by extensions that
    do not register their memory with Python's tracemalloc API.
    """
    if algorithm not in VALID_ALGORITHMS:
        raise ValueError(
            f"Unknown algorithm: {algorithm!r}."
        )

    if timeout is not None and timeout <= 0:
        raise ValueError(
            "timeout must be positive or None."
        )

    context = mp.get_context("spawn")
    parent_connection, child_connection = (
        context.Pipe(duplex=False)
    )

    process = context.Process(
        target=_tracemalloc_worker,
        args=(
            algorithm,
            list(nodes),
            list(edges),
            int(max_depth),
            child_connection,
        ),
    )

    deadline = (
        None
        if timeout is None
        else time.monotonic() + timeout
    )

    process.start()
    child_connection.close()

    try:
        message = _wait_for_worker_result(
            parent_connection,
            process,
            deadline,
            algorithm,
        )

        if not message.get("ok", False):
            raise RuntimeError(
                message.get(
                    "error",
                    f"Unknown tracemalloc error for {algorithm!r}.",
                )
            )

        process.join(
            timeout=_remaining_timeout(deadline)
        )

        if process.is_alive():
            process.terminate()
            process.join(timeout=5.0)
            raise TimeoutError(
                "The tracemalloc worker did not exit after returning "
                f"its result. algorithm={algorithm!r}."
            )

        if process.exitcode != 0:
            raise RuntimeError(
                "The tracemalloc worker exited abnormally. "
                f"algorithm={algorithm!r}, "
                f"exit_code={process.exitcode}."
            )

        current_bytes = int(
            message["tracemalloc_current_bytes"]
        )
        peak_bytes = int(
            message["tracemalloc_peak_bytes"]
        )
        overhead_bytes = int(
            message["tracemalloc_overhead_bytes"]
        )

        return {
            "algorithm": algorithm,
            "betti_number": int(message["betti_number"]),
            "tracemalloc_current_mib": current_bytes / MIB,
            "tracemalloc_peak_mib": peak_bytes / MIB,
            "tracemalloc_peak_minus_current_mib": (
                peak_bytes - current_bytes
            ) / MIB,
            "tracemalloc_overhead_mib": (
                overhead_bytes / MIB
            ),
            "runtime_seconds": float(
                message["runtime_seconds"]
            ),
        }

    except BaseException:
        if process.is_alive():
            process.terminate()

        process.join(timeout=5.0)
        raise

    finally:
        parent_connection.close()

        if not process.is_alive():
            try:
                process.close()
            except ValueError:
                pass

def random_stratified(nodes_per_layer: list[int], edges_between_layers: list[int], seed: int| None = None) -> nx.DiGraph:
    """
    Generate a random stratified (DAG without multi-edge) graph given the number of layers and nodes per layer.
    
    Args:
    nodes_per_layer (list of int): List of integers specifying the number of nodes in each layer.
    edges_between_layers (list of int): List of integers specifying the number of edges between each pair of adjacent layers.
    
    Returns:
    A NetworkX DiGraph representing the generated stratified DAG.
    """
    if seed is not None:
        random.seed(seed)
        np.random.seed(seed)
    L = len(nodes_per_layer)
    G = nx.DiGraph()
    node_counter = 0
    layers = []

    for i in range(L):
        layer_nodes = [f"n{node_counter + j}" for j in range(nodes_per_layer[i])]
        layers.append(layer_nodes)
        G.add_nodes_from(layer_nodes)
        node_counter += nodes_per_layer[i]

    for i in range(L - 1):
        layer_u = layers[i]
        layer_v = layers[i + 1]
        num_edges = edges_between_layers[i]
        max_possible_edges = len(layer_u) * len(layer_v)
        
        assert num_edges <= max_possible_edges, f"Too many edges between layers {i} and {i + 1}."
        
        possible_edges = [(u, v) for u in layer_u for v in layer_v]
        random.shuffle(possible_edges)
        
        selected_edges = possible_edges[:num_edges]
        G.add_edges_from(selected_edges) 

    assert nx.is_directed_acyclic_graph(G), "Graph is not acyclic."   
    
    return G

def _plot_stratified_comparison(
    densities: list[float],
    general_values: list[float],
    recursive_values: list[float],
    general_errors: list[float],
    recursive_errors: list[float],
    ylabel: str,
    title: str,
    filename: str,
    decimals: int = 3,
) -> None:
    """Plot two means with standard-deviation error bars and labels."""
    color_general = "tab:red"
    color_recursive = "tab:blue"

    figure, axis = plt.subplots()
    axis.set_xlabel("Edge Density")
    axis.set_ylabel(ylabel)

    axis.errorbar(
        densities,
        general_values,
        yerr=general_errors,
        color=color_general,
        marker="o",
        capsize=3,
        label="General Algorithm",
    )
    axis.errorbar(
        densities,
        recursive_values,
        yerr=recursive_errors,
        color=color_recursive,
        marker="s",
        capsize=3,
        label="Recursive Algorithm",
    )

    for index, value in enumerate(general_values):
        axis.annotate(
            f"{value:.{decimals}f}",
            (densities[index], value),
            textcoords="offset points",
            xytext=(0, 6),
            ha="center",
            color=color_general,
        )

    for index, value in enumerate(recursive_values):
        axis.annotate(
            f"{value:.{decimals}f}",
            (densities[index], value),
            textcoords="offset points",
            xytext=(0, -14),
            ha="center",
            color=color_recursive,
        )

    axis.legend(loc="upper left")
    axis.set_title(title, wrap=True)
    figure.savefig(filename, bbox_inches="tight")
    plt.close(figure)

def compare_stratified(
    Rep_num: int,
    nodes_per_layer: list[int],
    densities: list[float],
    seed: int | None = None,
) -> None:
    """
    Compare total runtimes only. No memory profiler or subprocess is used.

    Each realization uses the same graph for all three algorithms. Their
    execution order is cyclically rotated to avoid giving one algorithm a
    systematic order advantage.

    For each edge density, the plotted value is the sum of the runtimes over
    all Rep_num realizations, rather than the mean runtime per realization.
    """
    if Rep_num <= 0:
        raise ValueError("Rep_num must be positive.")

    depth = len(nodes_per_layer) - 1

    edge_counts = [
        [
            int(
                nodes_per_layer[layer]
                * density
                * nodes_per_layer[layer + 1]
            )
            for layer in range(depth)
        ]
        for density in densities
    ]

    algorithms = (
        "general",
        "recursive_basis",
        "recursive_no_basis",
    )

    # Each entry is the total runtime over Rep_num realizations.
    runtime_totals = {
        name: []
        for name in algorithms
    }

    raw_rows: list[list[object]] = []

    for density_index in trange(
        len(densities),
        desc="Runtime",
    ):
        runtime_runs = {
            name: []
            for name in algorithms
        }

        for repetition in range(Rep_num):
            current_seed = (
                seed + repetition
                if seed is not None
                else None
            )

            graph = random_stratified(
                nodes_per_layer,
                edge_counts[density_index],
                seed=current_seed,
            )

            nodes = list(graph.nodes())
            edgelist = list(graph.edges())
            max_L = nx.dag_longest_path_length(graph)

            shift = repetition % len(algorithms)
            execution_order = (
                algorithms[shift:]
                + algorithms[:shift]
            )

            betti_numbers: dict[str, int] = {}

            for algorithm in execution_order:
                if algorithm == "general":
                    algorithm_input = build_general_graph(
                        nodes,
                        edgelist,
                    )

                    start = time.perf_counter()

                    result = run_general_algorithm(
                        algorithm_input,
                        max_L,
                    )

                    elapsed = time.perf_counter() - start
                    betti_number = int(result[max_L])

                else:
                    algorithm_input = nx.DiGraph()
                    algorithm_input.add_nodes_from(nodes)
                    algorithm_input.add_edges_from(edgelist)

                    calculate_basis = (
                        algorithm == "recursive_basis"
                    )

                    start = time.perf_counter()

                    result = max_path_homology(
                        algorithm_input,
                        calculate_basis=calculate_basis,
                        report_sparsity=False,
                    )

                    elapsed = time.perf_counter() - start
                    betti_number = int(result[1])

                runtime_runs[algorithm].append(elapsed)
                betti_numbers[algorithm] = betti_number

                raw_rows.append(
                    [
                        densities[density_index],
                        repetition,
                        current_seed,
                        algorithm,
                        elapsed,
                        betti_number,
                    ]
                )

                del result
                del algorithm_input

            if len(set(betti_numbers.values())) != 1:
                raise RuntimeError(
                    "Homology mismatch at density "
                    f"{densities[density_index]}, "
                    f"repetition {repetition}: "
                    f"{betti_numbers}"
                )

        # Sum all Rep_num runtimes for each algorithm.
        for algorithm in algorithms:
            runtime_totals[algorithm].append(
                float(
                    np.sum(
                        runtime_runs[algorithm]
                    )
                )
            )

    # _plot_stratified_comparison expects error arrays.
    # Total runtimes have no error bars here, so use zeros.
    zero_errors = [
        0.0
        for _ in densities
    ]

    _plot_stratified_comparison(
        densities,
        runtime_totals["general"],
        runtime_totals["recursive_basis"],
        zero_errors,
        zero_errors,
        "Total Runtime (seconds)",
        (
            "Total Runtime with Basis Tracking "
            f"({Rep_num} realizations per density)"
        ),
        f"graph{depth}_comparison_overall_enable.pdf",
    )

    _plot_stratified_comparison(
        densities,
        runtime_totals["general"],
        runtime_totals["recursive_no_basis"],
        zero_errors,
        zero_errors,
        "Total Runtime (seconds)",
        (
            "Total Runtime without Basis Tracking "
            f"({Rep_num} realizations per density)"
        ),
        f"graph{depth}_comparison_overall_disable.pdf",
    )

    with open(
        f"graph{depth}_runtime_raw.csv",
        "w",
        newline="",
        encoding="utf-8",
    ) as file:
        writer = csv.writer(file)

        writer.writerow(
            [
                "density",
                "repetition",
                "seed",
                "algorithm",
                "runtime_seconds",
                "betti_number",
            ]
        )

        writer.writerows(raw_rows)

    # Also save the summed runtime for every density and algorithm.
    with open(
        f"graph{depth}_runtime_total.csv",
        "w",
        newline="",
        encoding="utf-8",
    ) as file:
        writer = csv.writer(file)

        writer.writerow(
            [
                "density",
                "algorithm",
                "total_runtime_seconds",
                "number_of_realizations",
            ]
        )

        for density_index, density in enumerate(densities):
            for algorithm in algorithms:
                writer.writerow(
                    [
                        density,
                        algorithm,
                        runtime_totals[algorithm][density_index],
                        Rep_num,
                    ]
                )

def compare_stratified_memory(
    Rep_num: int,
    nodes_per_layer: list[int],
    densities: list[float],
    seed: int | None = None,
    sample_interval: float | None = None,
    process_timeout: float | None = None,
) -> None:
    """
    Compare peak tracemalloc allocations for three implementations.

    Every algorithm run uses a fresh spawned process. In that process, module
    imports and algorithm-specific input construction finish before
    tracemalloc starts. The reported peak therefore covers allocations 
    made during the first formal algorithm call on the prepared target graph, 
    including the returned result while it remains alive.

    ``tracemalloc_peak_mib`` is not operating-system RSS. It measures the peak
    total size of allocations registered with Python's tracemalloc machinery.
    Native allocations from extensions that do not register with tracemalloc
    may not be included.

    ``sample_interval`` is accepted only for backward compatibility with
    older experiment scripts. Tracemalloc is event-based rather than sampled,
    so this argument is ignored.

    Error bars are one sample standard deviation across independently
    generated graph realizations. Raw values and summary statistics are saved
    to CSV files.
    """
    del sample_interval

    if (process_timeout is not None and process_timeout <= 0):
        raise ValueError("process_timeout must be positive or None.")

    depth = len(nodes_per_layer) - 1

    edge_counts = [
        [
            int(nodes_per_layer[layer] * density * nodes_per_layer[layer + 1]) for layer in range(depth)
        ]
        for density in densities
    ]

    algorithms = ("general", "recursive_basis", "recursive_no_basis")

    traced_peak_means = {name: [] for name in algorithms}
    traced_peak_stds = {name: [] for name in algorithms}
    traced_peak_medians = {name: [] for name in algorithms}
    traced_peak_q1 = {name: [] for name in algorithms}
    traced_peak_q3 = {name: [] for name in algorithms}
    traced_peak_maxima = {name: [] for name in algorithms}

    raw_rows: list[list[object]] = []

    for density_index in trange(len(densities), desc="Tracemalloc peak",):
        traced_peak_runs = {name: [] for name in algorithms}

        for repetition in range(Rep_num):
            current_seed = (seed + repetition if seed is not None else None)

            graph = random_stratified(
                nodes_per_layer,
                edge_counts[density_index],
                seed=current_seed,
            )

            nodes = list(graph.nodes())
            edgelist = list(graph.edges())
            max_L = nx.dag_longest_path_length(graph)

            betti_numbers: dict[str, int] = {}

            shift = repetition % len(algorithms)
            execution_order = (algorithms[shift:] + algorithms[:shift])

            for algorithm in execution_order:
                measurement = measure_algorithm_traced_memory(
                    algorithm=algorithm,
                    nodes=nodes,
                    edges=edgelist,
                    max_depth=max_L,
                    timeout=process_timeout,
                )

                traced_peak = float(measurement["tracemalloc_peak_mib"])
                traced_current = float(measurement["tracemalloc_current_mib"])
                peak_minus_current = float(measurement["tracemalloc_peak_minus_current_mib"])
                tracer_overhead = float(measurement["tracemalloc_overhead_mib"])
                runtime_seconds = float(measurement["runtime_seconds"])
                betti_number = int(measurement["betti_number"])

                traced_peak_runs[algorithm].append(traced_peak)
                betti_numbers[algorithm] = betti_number

                raw_rows.append(
                    [
                        densities[density_index],
                        repetition,
                        current_seed,
                        algorithm,
                        traced_peak,
                        traced_current,
                        peak_minus_current,
                        tracer_overhead,
                        runtime_seconds,
                        betti_number,
                    ]
                )

            if len(set(betti_numbers.values())) != 1:
                raise RuntimeError(
                    "Homology mismatch in the tracemalloc experiment "
                    f"at density {densities[density_index]}, "
                    f"repetition {repetition}: {betti_numbers}"
                )

        for algorithm in algorithms:
            values = np.asarray(traced_peak_runs[algorithm], dtype=float)

            traced_peak_means[algorithm].append(float(np.mean(values)))
            traced_peak_stds[algorithm].append(float(np.std(values, ddof=1)) if Rep_num > 1 else 0.0)
            traced_peak_medians[algorithm].append(float(np.median(values)))
            traced_peak_q1[algorithm].append(float(np.quantile(values, 0.25)))
            traced_peak_q3[algorithm].append(float(np.quantile(values, 0.75)))
            traced_peak_maxima[algorithm].append(float(np.max(values)))

    _plot_stratified_comparison(
        densities,
        traced_peak_means["general"],
        traced_peak_means["recursive_basis"],
        traced_peak_stds["general"],
        traced_peak_stds["recursive_basis"],
        "Peak Traced Memory Allocations (MiB)",
        (
            "Peak Traced Allocations with Basis Tracking "
            f"({Rep_num} isolated realizations per density)"
        ),
        f"graph{depth}_tracemalloc_enable.pdf",
        decimals=3,
    )

    _plot_stratified_comparison(
        densities,
        traced_peak_means["general"],
        traced_peak_means["recursive_no_basis"],
        traced_peak_stds["general"],
        traced_peak_stds["recursive_no_basis"],
        "Peak Traced Memory Allocations (MiB)",
        (
            "Peak Traced Allocations without Basis Tracking "
            f"({Rep_num} isolated realizations per density)"
        ),
        f"graph{depth}_tracemalloc_disable.pdf",
        decimals=3,
    )

    with open(f"graph{depth}_tracemalloc_raw.csv", "w", newline="", encoding="utf-8") as file:
        writer = csv.writer(file)
        writer.writerow(
            [
                "density",
                "repetition",
                "seed",
                "algorithm",
                "tracemalloc_peak_mib",
                "tracemalloc_current_mib",
                "tracemalloc_peak_minus_current_mib",
                "tracemalloc_overhead_mib",
                "runtime_seconds_with_tracing",
                "betti_number",
            ]
        )
        writer.writerows(raw_rows)

    with open(f"graph{depth}_tracemalloc_summary.csv", "w", newline="", encoding="utf-8") as file:
        writer = csv.writer(file)
        writer.writerow(
            [
                "density",
                "algorithm",
                "tracemalloc_peak_mean_mib",
                "tracemalloc_peak_sd_mib",
                "tracemalloc_peak_median_mib",
                "tracemalloc_peak_q1_mib",
                "tracemalloc_peak_q3_mib",
                "tracemalloc_peak_max_mib",
                "number_of_realizations",
            ]
        )

        for density_index, density in enumerate(densities):
            for algorithm in algorithms:
                writer.writerow(
                    [
                        density,
                        algorithm,
                        traced_peak_means[algorithm][density_index],
                        traced_peak_stds[algorithm][density_index],
                        traced_peak_medians[algorithm][density_index],
                        traced_peak_q1[algorithm][density_index],
                        traced_peak_q3[algorithm][density_index],
                        traced_peak_maxima[algorithm][density_index],
                        Rep_num,
                    ]
                )

def recursive_stratified(Rep_num: int, nodes_per_layer: list[int], densities: list[float], seed: int = None) -> None:
    """
    Record and plot the time spent on the Recursive Algorithm for random stratified graphs with different edge densities.

    Args:
        Rep_num (int):                 Number of repetitions per edge density.
        seed (int):                    Random seed for reproducibility.
        nodes_per_layer (list[int]):   Number of nodes in each layer.
        densities (list[float]):       Edge densities corresponding to the edge counts.

    Plot: Overall Time Spent on Recursive Algorithm at Different Edge Densities without Basis Calculation.
    """
    Total_Time_Hierarchical_withoutbasis = []
    L = len(nodes_per_layer)
    edges_between_layers = []
    for density in densities:
        N_edges = [int(nodes_per_layer[i] * density * nodes_per_layer[i + 1]) for i in range(L - 1)]
        edges_between_layers.append(N_edges)

    for j in trange(len(densities)):
        edges = edges_between_layers[j]
        Time_hierarchical_withoutbasis = 0

        for i in range(Rep_num):
            current_seed = seed + i if seed is not None else None
            G = random_stratified(nodes_per_layer, edges, seed = current_seed)
            start_withoutbasis = time.perf_counter()
            _, _, _ = max_path_homology(G, calculate_basis = False, report_sparsity = False)
            stop_withoutbasis = time.perf_counter()
            Time_hierarchical_withoutbasis += stop_withoutbasis - start_withoutbasis

        Total_Time_Hierarchical_withoutbasis.append(Time_hierarchical_withoutbasis)

    color_hierarchical = 'tab:blue'
    ############### Plots without basis calculation ################
    _, ax1 = plt.subplots()
    ax1.set_xlabel('Edge Density')
    ax1.set_ylabel('Total Time Spent (seconds)', color='black')
    ax1.plot(densities, Total_Time_Hierarchical_withoutbasis, color=color_hierarchical, marker='s', label='Recursive Algorithm')
    ax1.tick_params(axis='y')

    for i, txt in enumerate(Total_Time_Hierarchical_withoutbasis):
        ax1.annotate(f'{txt:.3f}', (densities[i], Total_Time_Hierarchical_withoutbasis[i]), textcoords="offset points", xytext=(0,-12), ha='center', color=color_hierarchical)

    ax1.legend(loc='upper left')
    plt.title(f"Overall Time Spent on Recursive Algorithm at Different Edge Densities ({Rep_num} realizations each), basis calculation disabled", wrap = True)
    plt.savefig(f'gamma{L-1}_without_basis_track.pdf')
    #plt.show()
    plt.clf()

def random_dag(n: int, p:float, seed: int = None) -> nx.DiGraph:
    """
    Generate a random DAG with n nodes.

    Each possible edge n_i -> n_j with i < j is included
    independently with probability p.
    
    Args:
    n (int): Number of nodes.
    p (float): Probability of including each possible directed edge.
    seed (int): Random seed for reproducibility.
    
    Returns:
    A NetworkX DiGraph.
    """
    if seed is not None:
        random.seed(seed)
        np.random.seed(seed)

    nodes = [f"n{i}" for i in range(n)]

    G = nx.DiGraph()
    G.add_nodes_from(nodes)

    for i in range(n - 1):
        for j in range(i + 1, n):
            if random.random() < p:
                G.add_edge(nodes[i], nodes[j])
    return G

# Example usage
if __name__ == '__main__':
    stratified_graph = random_stratified([5, 5, 5], [10, 10], seed=1)

    print("randomized stratified graph nodes: " f"{list(stratified_graph.nodes)}")
    print("randomized stratified graph edges: " f"{list(stratified_graph.edges)}\n")

    dag = random_dag(n = 15, p = 0.5, seed=1)

    print(f"randomized DAG nodes: {list(dag.nodes)}")
    print(f"randomized DAG edges: {list(dag.edges)}\n")