import random

import numpy as np
import networkx as nx
import matplotlib.pyplot as plt

from tqdm import tqdm
from experiment_func import random_stratified
from recursive_algorithm import max_path_homology

def persistent_homology(nodes_per_layer: list[int], edges_between_layers: list[int], seed: int = None) -> None:
    """
    Performs decreasing persistent homology plot on a random stratified graph given edge weights sampled from a uniform distribution.

    Args:
        nodes_per_layer (list[int]):   Number of nodes in each layer.
        edges_between_layers (list[int]): Number of edges between layers.
    """
    if seed is not None:
        random.seed(seed)
        np.random.seed(seed)
    
    G = random_stratified(nodes_per_layer, edges_between_layers, seed = seed)
    lp = nx.dag_longest_path_length(G)

    nodes = list(G.nodes())

    weighted_edgelist = [(u, v, random.uniform(0, 1)) for u, v in G.edges()]
    weighted_edgelist.sort(key=lambda edge: edge[2])
    edges = [(u, v) for u, v, _ in weighted_edgelist]
    weights = [w for _, _, w in weighted_edgelist]

    betti_numbers = []
    thresholds = []
    threshold = 0

    with tqdm(total = len(edges), desc = "Removing edges and computing Betti numbers", unit = "edge") as progress_bar:
        while edges:
            G_threshold = nx.DiGraph()
            G_threshold.add_nodes_from(nodes)
            G_threshold.add_edges_from(edges)

            current_lp, betti, _ = max_path_homology(G_threshold, calculate_basis = False)

            thresholds.append(threshold)
            if current_lp == lp:
                betti_numbers.append(betti)
            else:
                betti_numbers.append(0)
                break   

            if not edges:
                break

            edges.pop(0)
            threshold = weights.pop(0)
            progress_bar.update(1)

    thresholds.append(1)
    betti_numbers.append(0)
    
    plt.plot(thresholds, betti_numbers)
    plt.xlabel('t')
    plt.ylabel('Betti Number')
    #plt.savefig('4_10_10_10 persistent_homology standard uniform weights.pdf', bbox_inches='tight')
    plt.show()
    plt.clf()

if __name__ == '__main__':
    persistent_homology([4, 10, 10, 10], [40, 100, 100], seed = 1)