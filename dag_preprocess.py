import networkx as nx

from collections import deque
from collections.abc import Hashable

def layer_dict(G: nx.DiGraph) -> dict[int, set[Hashable]]:
    """
    Label the nodes of a stratified graph. 
    The key is the depth of the layer, and the corresponding value is the set of nodes in that layer.

    Args:
        G (nx.DiGraph): A directed graph.

    Returns:
        layer_dict (dict[int, set[str]]): A dictionary that maps each layer to a set of nodes in that layer.
    """
    roots = [node for node in G if G.in_degree(node) == 0]
    visited = set()
    layer_dict = {}
    stack = [(root, 0) for root in roots]

    while stack:
        node, depth = stack.pop()
        if node in visited:
            continue
        visited.add(node)
        layer_dict.setdefault(depth, set()).add(node)
        for succ in G.successors(node):
            if succ not in visited:
                stack.append((succ, depth + 1))

    return layer_dict

def graph_decomposition(G: nx.DiGraph) -> nx.DiGraph:
    """
    The algorithm 2. Given the DAG G, construct a stratified graph of G that contains all the longest paths of G.

    Args:
        G (nx.DiGraph): A DAG.

    Returns:
        nx.DiGraph: A stratified graph that contains all the longest paths of G.
    """

    if not nx.is_directed_acyclic_graph(G):
        raise ValueError("The input graph must be a DAG.")
    
    G = G.copy()
    lp = nx.dag_longest_path_length(G)

    if lp == 0:
        return G
    
    topological_order = list(nx.topological_sort(G))
    from_top = {v: 0 for v in G.nodes}
    from_bottom = {v: 0 for v in G.nodes}

    for u in topological_order:
        for v in G.successors(u):
            from_top[v] = max(from_top[v], from_top[u] + 1)
    
    for v in reversed(topological_order):
        for u in G.predecessors(v):
            from_bottom[u] = max(from_bottom[u], from_bottom[v] + 1)
    
    edges_to_keep = [
        (u, v)
        for u, v in G.edges
        if from_top[u] + 1 + from_bottom[v] == lp
    ]

    G_star = nx.DiGraph(edges_to_keep)
    return G_star

def prune(G: nx.DiGraph) -> nx.DiGraph:
    """
    Recursively prune the graph from graph_decomposition.
    Pruning rule: a node and its connected edges are removed if it has only one successor or one predecessor.

    Args:
        G (nx.DiGraph):   G* from graph_decomposition with positive length.

    Returns:
        nx.DiGraph: A directed subgraph pruned from G*.
    """
    G_prune = G.copy()

    def is_prunable(node):
        return (G_prune.in_degree(node) == 1 or G_prune.out_degree(node) == 1)
    
    queue = deque(node for node in G_prune.nodes if is_prunable(node))
    plan = set(queue)

    while queue:
        node = queue.popleft()
        if node not in G_prune:
            continue
        predecessors = list(G_prune.predecessors(node))
        successors = list(G_prune.successors(node))
        G_prune.remove_node(node)

        for neighbor in predecessors + successors:
            if (neighbor in G_prune and neighbor not in plan and is_prunable(neighbor)):
                queue.append(neighbor)
                plan.add(neighbor)

    return G_prune

def dag_process(G: nx.DiGraph) -> tuple[list[dict[int, set[Hashable]]], list[list[int]], int, int, list[nx.DiGraph],]:
    """
    Preprocess a finite DAG for maximal path-homology computation.

    For positive maximum path length, construct G_*, recursively
    prune it, and separate the remaining weakly connected components.
    """

    if not isinstance(G, nx.DiGraph) or isinstance(G, nx.MultiDiGraph):
        raise TypeError("G must be a NetworkX DiGraph.")

    if G.number_of_nodes() == 0:
        raise ValueError("G must contain at least one node.")

    if not nx.is_directed_acyclic_graph(G):
        raise ValueError("G must be a DAG.")

    G_copy = G.copy()
    lp = nx.dag_longest_path_length(G_copy)

    if lp == 0:
        layers = {0: set(G_copy.nodes())}
        return [layers], [[G_copy.number_of_nodes()]], lp, 1, [G_copy]
    else:
        G_decomposed = graph_decomposition(G_copy)
        G_pruned = prune(G_decomposed)
        lp_pruned = nx.dag_longest_path_length(G_pruned)

        if lp != lp_pruned:
            return [], [], lp, 0, []
        else:
            components = [G_pruned.subgraph(c).copy() for c in nx.weakly_connected_components(G_pruned)]
            subgraph_dict, node_counts, graph_list = [], [], []
            num_graph = 0
            for subgraph in components:
                layers = layer_dict(subgraph)
                num_graph += 1
                s_dict = {K: layers[K] for K in sorted(layers)}
                node_count = [len(layers[i]) for i in range(len(layers))]
                graph_list.append(subgraph)
                subgraph_dict.append(s_dict)
                node_counts.append(node_count)
        
            return (subgraph_dict, node_counts, lp, num_graph, graph_list)

# An example usage
if __name__ == "__main__":
    edgelist = [('a0', 'b2'), ('a0', 'b3'), ('a1', 'b2'), ('a1', 'b3'), 
                ('a1', 'c1'), ('b0', 'c0'), ('b0', 'c1'), ('b1', 'c0'), ('b1', 'c1'), 
                ('b2', 'c2'), ('b2', 'c3'), ('b3', 'c2'), ('b3', 'c3'), ('b4', 'c4'), 
                ('b4', 'c5'), ('b5', 'c4'), ('b5', 'c5'), ('b0', 'c2'), ('b0', 'c3'), 
                ('b1', 'c2'), ('b1', 'c3'), ('b4', 'c2'), ('b4', 'c3'), ('b5', 'c2'), 
                ('b5', 'c3'), ('c0', 'd0'), ('c0', 'd1'), ('c1', 'd0'), ('c1', 'd1'), 
                ('c4', 'd2'), ('c4', 'd3'), ('c5', 'd2'), ('c5', 'd3')]

    G = nx.DiGraph(edgelist)
    lp = nx.dag_longest_path_length(G)

    print(f"original graph:{G.edges()},\n number of edges = {G.number_of_edges()},\n max path length:{lp}\n")
    (subgrph_dict, node_counts, lp, num_graph, graph_list) = dag_process(G)
    print(f"subgrph_dict:{subgrph_dict},\n node_counts:{node_counts},\n lp:{lp},\n num_graph:{num_graph},\n graph_list:{graph_list}\n")

    iso_nodes = ['a0', 'a1', 'a2', 'b0', 'b1', 'b2', 'b3', 'b4', 'b5', 'c0', 'c1', 'c2', 'c3', 'c4', 'c5', 'd0', 'd1', 'd2', 'd3']
    G1 = nx.DiGraph()
    G1.add_nodes_from(iso_nodes)
    lp1 = nx.dag_longest_path_length(G1)
    print(f"original graph:{G1.edges() if G1.edges() else G1.nodes()},\n number of edges = {G1.number_of_edges()},\n max path length:{lp1}\n")
    (subgrph_dict1, node_counts1, lp1, num_graph1, graph_list1) = dag_process(G1)
    print(f"subgrph_dict:{subgrph_dict1},\n node_counts:{node_counts1},\n lp:{lp1},\n num_graph:{num_graph1},\n graph_list:{graph_list1}\n")