import numpy as np
import random
import matplotlib.pyplot as plt
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
        edgelist = random_stratified(nodes_per_layer, edges_between_layers, seed = seed)
    else:
        edgelist = random_stratified(nodes_per_layer, edges_between_layers)

    weighted_edgelist = edgelist.copy()
    for idx, edge in enumerate(weighted_edgelist):
        weight = random.uniform(0, 1)
        u, v = edge
        weighted_edgelist[idx] = (u, v, weight)
    
    weighted_edgelist = sorted(weighted_edgelist, key=lambda x: x[2])
    edges = [(u, v) for u, v, _ in weighted_edgelist]
    weights = [w for _, _, w in weighted_edgelist]
    betti_number = []
    thresholds = []
    threshold = 0
    while edges:
        _, results, _ = max_path_homology(edges, calculate_basis = False)
        betti_number.append(results)
        thresholds.append(threshold)
        edge = edges.pop(0)
        threshold = weights.pop(0)
    betti_number.append(0)
    thresholds.append(1)
    plt.plot(thresholds, betti_number)
    plt.xlabel('t')
    plt.ylabel('Betti Number')
    plt.savefig('4_10_10_10 persistent_homology standard uniform weights.pdf')
    plt.show()
    plt.clf()

if __name__ == '__main__':
    persistent_homology([4, 10, 10, 10], [40, 100, 100], seed = 1)