import networkx as nx
import random
import numpy as np
import matplotlib.pyplot as plt
import time
from tqdm import trange
from dag_preprocess import dag_process
from recursive_algorithm import max_path_homology
from general_algorithm import R_path, H_path_R, edgelist_to_graph

def random_stratified(nodes_per_layer: list[int], edges_between_layers: list[int], seed: int = None) -> list[tuple[str, str]]:
    """
    Generate a random stratified (directed and acyclic without multi-edge) graph given the number of layers and nodes per layer.
    
    Args:
    nodes_per_layer (list of int): List of integers specifying the number of nodes in each layer.
    edges_between_layers (list of int): List of integers specifying the number of edges between each pair of adjacent layers.
    
    Returns:
    edgelist (list of tuples): List of edges in the generated DAG.
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
    edgelist = list(G.edges)
    
    return edgelist

def compare_stratified(Rep_num: int, nodes_per_layer: list[int], densities: list[float], seed: int = None) -> None:
    """
    Compares the time spent on the General and Recursive Algorithms for different edge densities.
    Graphs are randomly generated stratified.

    Args:
        Rep_num (int):                 Number of repetitions per edge density.
        nodes_per_layer (list[int]):   Number of nodes in each layer.
        densities (list[float]):       Edge densities corresponding to the edge counts.
        seed (int):                    Random seed for reproducibility.

    Two plots are generated:
    Plot 1: Overall Time Spent on General and Hierarchical Algorithm at Different Edge Densities with Basis Calculation.
    Plot 2: Overall Time Spent on General and Hierarchical Algorithm at Different Edge Densities without Basis Calculation.
    """
    Total_Time_General = []
    Total_Time_Hierarchical_withbasis = []
    Total_Time_Hierarchical_withoutbasis = []
    L = len(nodes_per_layer)
    edges_between_layers = []
    for density in densities:
        N_edges = [int(nodes_per_layer[i] * density * nodes_per_layer[i + 1]) for i in range(L - 1)]
        edges_between_layers.append(N_edges)

    for j in trange(len(densities)):
        edges = edges_between_layers[j]
        Time_general = 0
        Time_hierarchical_withbasis = 0
        Time_hierarchical_withoutbasis = 0

        for i in range(Rep_num):
            current_seed = seed + i if seed is not None else None
            edgelist = random_stratified(nodes_per_layer, edges, seed = current_seed)
            general_edgelist = [f"{s}:{t}" for s, t in edgelist]
            G_com = edgelist_to_graph(general_edgelist)
            _, _, max_L, _, _ = dag_process(edgelist)

            starts = time.perf_counter()
            R = R_path(G_com, cutoff = max_L + 1)
            H, _, _, _, _ = H_path_R(R)
            ends = time.perf_counter()
            Time_general += ends - starts

            start_withbasis = time.perf_counter()
            _, Result_withbasis, _ = max_path_homology(edgelist, calculate_basis = True)
            stop_withbasis = time.perf_counter()
            Time_hierarchical_withbasis += stop_withbasis - start_withbasis

            start_withoutbasis = time.perf_counter()
            _, _, _ = max_path_homology(edgelist, calculate_basis = False)
            stop_withoutbasis = time.perf_counter()
            Time_hierarchical_withoutbasis += stop_withoutbasis - start_withoutbasis

            if H[max_L] != Result_withbasis:
                print(f"Mismatch in homology results at edge density {densities[j]} for repetition {i}")
                print(f"edgelist: {edgelist}")

        Total_Time_General.append(Time_general)
        Total_Time_Hierarchical_withbasis.append(Time_hierarchical_withbasis)
        Total_Time_Hierarchical_withoutbasis.append(Time_hierarchical_withoutbasis)

    ############### Plots with basis calculation ################
    _, ax1 = plt.subplots()
    color_general = 'tab:red'
    color_hierarchical = 'tab:blue'
    ax1.set_xlabel('Edge Density')
    ax1.set_ylabel('Total Time Spent (seconds)', color='black')
    ax1.plot(densities, Total_Time_General, color=color_general, marker='o', label='General Algorithm')
    ax1.plot(densities, Total_Time_Hierarchical_withbasis, color=color_hierarchical, marker='s', label='Recursive Algorithm')
    ax1.tick_params(axis='y')
    for i, txt in enumerate(Total_Time_General):
        ax1.annotate(f'{txt:.3f}', (densities[i], Total_Time_General[i]), textcoords="offset points", xytext=(0,4), ha='center', color=color_general)

    for i, txt in enumerate(Total_Time_Hierarchical_withbasis):
        ax1.annotate(f'{txt:.3f}', (densities[i], Total_Time_Hierarchical_withbasis[i]), textcoords="offset points", xytext=(0,-12), ha='center', color=color_hierarchical)

    ax1.legend(loc='upper left')
    plt.title(f"Overall Time Spent on General and Recursive Algorithm at Different Edge Densities ({Rep_num} realizations each), basis calculation enabled", wrap = True)
    plt.savefig(f'graph{L-1}_comparison_overall_enable.pdf')
    plt.show()

    ############### Plots without basis calculation ################
    _, ax1 = plt.subplots()
    ax1.set_xlabel('Edge Density')
    ax1.set_ylabel('Total Time Spent (seconds)', color='black')
    ax1.plot(densities, Total_Time_General, color=color_general, marker='o', label='General Algorithm')
    ax1.plot(densities, Total_Time_Hierarchical_withoutbasis, color=color_hierarchical, marker='s', label='Recursive Algorithm')
    ax1.tick_params(axis='y')
    for i, txt in enumerate(Total_Time_General):
        ax1.annotate(f'{txt:.3f}', (densities[i], Total_Time_General[i]), textcoords="offset points", xytext=(0,4), ha='center', color=color_general)

    for i, txt in enumerate(Total_Time_Hierarchical_withoutbasis):
        ax1.annotate(f'{txt:.3f}', (densities[i], Total_Time_Hierarchical_withoutbasis[i]), textcoords="offset points", xytext=(0,-12), ha='center', color=color_hierarchical)

    ax1.legend(loc='upper left')
    plt.title(f"Overall Time Spent on General and Recursive Algorithm at Different Edge Densities ({Rep_num} realizations each), basis calculation disabled", wrap = True)
    plt.savefig(f'graph{L-1}_comparison_overall_disable.pdf')
    plt.show()
    plt.clf()

def recursive_stratified(Rep_num: int, nodes_per_layer: list[int], densities: list[float], seed: int = None) -> None:
    """
    Record and plot the time spent on the Recursive Algorithm for different edge densities.
    Graphs are random stratified graphs.

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
            edgelist = random_stratified(nodes_per_layer, edges, seed = current_seed)
            start_withoutbasis = time.perf_counter()
            _, _, _ = max_path_homology(edgelist, calculate_basis = False)
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
    plt.show()
    plt.clf()

def random_dag(n: int, p:float, seed: int = None) -> list[tuple[str, str]]:
    """
    Generate a random DAG given the number of nodes and edge density.
    
    Args:
    n (int): Number of nodes.
    p (float): Probability to connect two nodes.
    
    Returns:
    edgelist (list[tuple[str, str]]): List of edges in the generated DAG.
    """
    if seed is not None:
        random.seed(seed)
        np.random.seed(seed)
    nodes = [f"n{i}" for i in range(n)]
    edgelist = []
    for i in range(n - 1):
        for j in range(i + 1, n):
            if random.random() < p:
                edgelist.append((nodes[i], nodes[j]))
    return edgelist

def compare_dag(Rep_num: int, n: int, p: list[float], seed: int = None) -> None:
    """
    Compares the time spent on the General and Recursive Algorithms for different edge densities.
    Graphs are random DAGs.

    Args:
        Rep_num (int):   Number of repetitions per edge density.
        n (int):         Number of nodes in each layer.
        p (List[float]): Probability to connect two nodes.
        seed (int):      Random seed for reproducibility.

    Two plots are generated:
    Plot 1: Overall Time Spent on General and Recursive Algorithm at Different probabilities with Basis Calculation.
    Plot 2: Overall Time Spent on General and Recursive Algorithm at Different Probabilities without Basis Calculation.
    """
    Total_Time_General = []
    Total_Time_Hierarchical_withbasis = []
    Total_Time_Hierarchical_withoutbasis = []

    for j in trange(len(p)):
        Time_general = 0
        Time_hierarchical_withbasis = 0
        Time_hierarchical_withoutbasis = 0

        for i in range(Rep_num):
            current_seed = seed + i if seed is not None else None
            edgelist = random_dag(n, p[j], current_seed)
            general_edgelist = [f"{s}:{t}" for s, t in edgelist]
            G_com = edgelist_to_graph(general_edgelist)
            _, _, max_L, _, _ = dag_process(edgelist)

            starts = time.perf_counter()
            R = R_path(G_com, cutoff = max_L + 1)
            H, _, _, _, _ = H_path_R(R)
            ends = time.perf_counter()
            Time_general += ends - starts

            start_withbasis = time.perf_counter()
            _, Result_withbasis, _ = max_path_homology(edgelist, calculate_basis = True)
            stop_withbasis = time.perf_counter()
            Time_hierarchical_withbasis += stop_withbasis - start_withbasis

            start_withoutbasis = time.perf_counter()
            _, _, _ = max_path_homology(edgelist, calculate_basis = False)
            stop_withoutbasis = time.perf_counter()
            Time_hierarchical_withoutbasis += stop_withoutbasis - start_withoutbasis

            if H[max_L] != Result_withbasis:
                print(f"Mismatch in homology results at edge density {p[j]} for repetition {i}")
                print(f"edgelist: {edgelist}, H: {H[max_L]}, Result_withbasis: {Result_withbasis}")

        Total_Time_General.append(Time_general)
        Total_Time_Hierarchical_withbasis.append(Time_hierarchical_withbasis)
        Total_Time_Hierarchical_withoutbasis.append(Time_hierarchical_withoutbasis)

    ############### Plots with basis calculation ################
    _, ax1 = plt.subplots()
    color_general = 'tab:red'
    color_hierarchical = 'tab:blue'
    ax1.set_xlabel('Edge Density')
    ax1.set_ylabel('Total Time Spent (seconds)', color='black')
    ax1.plot(p, Total_Time_General, color=color_general, marker='o', label='General Algorithm')
    ax1.plot(p, Total_Time_Hierarchical_withbasis, color=color_hierarchical, marker='s', label='Recursive Algorithm')
    ax1.tick_params(axis='y')
    for i, txt in enumerate(Total_Time_General):
        ax1.annotate(f'{txt:.3f}', (p[i], Total_Time_General[i]), textcoords="offset points", xytext=(0,4), ha='center', color=color_general)

    for i, txt in enumerate(Total_Time_Hierarchical_withbasis):
        ax1.annotate(f'{txt:.3f}', (p[i], Total_Time_Hierarchical_withbasis[i]), textcoords="offset points", xytext=(0,-12), ha='center', color=color_hierarchical)

    ax1.legend(loc='upper left')
    plt.title(f"Overall Time Spent on General and Recursive Algorithm at Different Edge Densities ({Rep_num} realizations each), basis calculation enabled", wrap = True)
    plt.savefig(f'randomdag_comparison_overall_enable_{n}nodes.pdf')
    plt.show()

    ############### Plots without basis calculation ################
    _, ax1 = plt.subplots()
    ax1.set_xlabel('Edge Density')
    ax1.set_ylabel('Total Time Spent (seconds)', color='black')
    ax1.plot(p, Total_Time_General, color=color_general, marker='o', label='General Algorithm')
    ax1.plot(p, Total_Time_Hierarchical_withoutbasis, color=color_hierarchical, marker='s', label='Recursive Algorithm')
    ax1.tick_params(axis='y')
    for i, txt in enumerate(Total_Time_General):
        ax1.annotate(f'{txt:.3f}', (p[i], Total_Time_General[i]), textcoords="offset points", xytext=(0,4), ha='center', color=color_general)

    for i, txt in enumerate(Total_Time_Hierarchical_withoutbasis):
        ax1.annotate(f'{txt:.3f}', (p[i], Total_Time_Hierarchical_withoutbasis[i]), textcoords="offset points", xytext=(0,-12), ha='center', color=color_hierarchical)

    ax1.legend(loc='upper left')
    plt.title(f"Overall Time Spent on General and Recursive Algorithm at Different Edge Densities ({Rep_num} realizations each), basis calculation disabled", wrap = True)
    plt.savefig(f'randomdag_comparison_overall_disable_{n}nodes.pdf')
    plt.show()
    plt.clf()

def recursive_dag(Rep_num: int, n: int, p: list[float], seed: int = None) -> None:
    """
    Record and plot the time spent on the Recursive Algorithm for different edge densities.

    Graphs are random DAGs.

    Args:
        Rep_num (int):   Number of repetitions per edge density.
        n (int):         Number of nodes in each layer.
        p (List[float]): Probability to connect two nodes.
        seed (int):      Random seed for reproducibility.

    Plot: Overall Time Spent on General and Recursive Algorithm at Different Probabilities without Basis Calculation.
    """
    Total_Time_Hierarchical_withoutbasis = []

    for j in trange(len(p)):
        Time_hierarchical_withoutbasis = 0

        for i in range(Rep_num):
            current_seed = seed + i if seed is not None else None
            edgelist = random_dag(n, p[j], seed = current_seed)
            start_withoutbasis = time.perf_counter()
            _, _, _ = max_path_homology(edgelist, calculate_basis = False)
            stop_withoutbasis = time.perf_counter()
            Time_hierarchical_withoutbasis += stop_withoutbasis - start_withoutbasis

        Total_Time_Hierarchical_withoutbasis.append(Time_hierarchical_withoutbasis)

    color_hierarchical = 'tab:blue'
    ############### Plots without basis calculation ################
    _, ax1 = plt.subplots()
    ax1.set_xlabel('Edge Density')
    ax1.set_ylabel('Total Time Spent (seconds)', color='black')
    ax1.plot(p, Total_Time_Hierarchical_withoutbasis, color=color_hierarchical, marker='s', label='Recursive Algorithm')
    ax1.tick_params(axis='y')

    for i, txt in enumerate(Total_Time_Hierarchical_withoutbasis):
        ax1.annotate(f'{txt:.3f}', (p[i], Total_Time_Hierarchical_withoutbasis[i]), textcoords="offset points", xytext=(0,-12), ha='center', color=color_hierarchical)

    ax1.legend(loc='upper left')
    plt.title(f"Overall Time Spent on Recursive Algorithm at Different Edge Densities ({Rep_num} realizations each), basis calculation disabled", wrap = True)
    plt.savefig(f'randomdag_recursive_overall_disable_{n}nodes.pdf')
    plt.show()
    plt.clf()

# Example of generating random graph
if __name__ == '__main__':
    graph = random_stratified([5, 5, 5], [10, 10], seed = 1)
    print(f"randomized stratified graph: {graph}\n")

    dag_edgelist = random_dag(15, 0.5, seed = 1)
    print(f"randomized dag: {dag_edgelist}\n")