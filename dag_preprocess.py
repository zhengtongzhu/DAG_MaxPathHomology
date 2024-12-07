import networkx as nx

def edge_density(G: nx.DiGraph, layer0: set[str], layer1: set[str]) -> tuple[float, int]:
    """
    Calculates the density of the edge between 0th layer and 1st layer.

    Args:
        G (nx.DiGraph): A directed graph.
        layer0 (set[str]): Set of nodes in layer 0.
        layer1 (set[str]): Set of nodes in layer 1.

    Returns:
        density (float): The density of the edge between 0th layer and 1st layer.
        num_edges (int): The number of edges between 0th layer and 1st layer.
    """
    num_edges = sum(1 for u, v in G.edges() if u in layer0 and v in layer1)
    density = num_edges / (len(layer0) * len(layer1))
    return density, num_edges

def layer_dict(G: nx.DiGraph) -> dict[int, set[str]]:
    """
    Label the nodes in the stratified graph. Key is the depth of the layer, value is the set of nodes in that layer.

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

def prune(G: nx.DiGraph) -> nx.DiGraph:
    """
    Recursively prune the graph.
    Pruning rule: a node and edges connected to it will be removed if the node has only one successor or one predecessor.

    Args:
        G (nx.DiGraph):   A directed graph.

    Returns:
        nx.DiGraph: A directed subgraphs that were pruned from the original graph.
    """
    while True:
        to_remove = {node for node in G if G.in_degree(node) == 1 or G.out_degree(node) == 1}
        if not to_remove:
            break
        G.remove_nodes_from(to_remove)
    return G

def lp_edgelist(G: nx.DiGraph) -> nx.DiGraph:
    """
    The G* algorithm. Given the DiGraph G, return the DiGraph that contains all the edges in longest paths.
    Isolated nodes are removed.

    Args:
        G (nx.DiGraph): A directed graph.

    Returns:
        nx.DiGraph: A DiGraph that contains all the edges in longest paths.
    """
    global lp
    if lp == 0:
        return G
    topological_sort = list(nx.topological_sort(G))
    from_top = {node: 0 for node in G}
    from_bottom = from_top.copy()

    for node in topological_sort:
        for succ in G.successors(node):
            if from_top[succ] < from_top[node] + 1:
                from_top[succ] = from_top[node] + 1
    
    for node in reversed(topological_sort):
        for pred in G.predecessors(node):
            if from_bottom[pred] < from_bottom[node] + 1:
                from_bottom[pred] = from_bottom[node] + 1
    
    edge_to_remove = [(u, v) for u, v in G.edges() if from_top[u] + 1 + from_bottom[v] != lp]
    G.remove_edges_from(edge_to_remove)
    isolated_nodes = list(nx.isolates(G))
    G.remove_nodes_from(isolated_nodes)
    return G

def dag_process(edgelist: list[tuple[str, str]], iso_nodes: list[str] = None) -> tuple[list[dict[int, set[str]]], list[int], int, int, list[nx.DiGraph]]:
    """
    This function takes an edgelist representing a directed acyclic graph (DAG) without multi-edges, dividing it into subgraphs 
    that are equivalent for calculating maximal path homology. If the edge density between the top two layers exceeds that between 
    the bottom two layers, the graph is reversed for faster calculations.

    Args:
        edgelist (list[tuple[str, str]]): A list of tuples where each tuple represents a directed edge in the graph.
        iso_nodes (list[str]):            A list of nodes that are isolated in the graph.

    Returns:
        Tuple:
            subgraph_dict (list[dict[int, set[str]]]): A list of dictionaries where each dictionary represents the nodes grouped by 
                                                       their layer number within a subgraph. The keys are the layer numbers and 
                                                       the values are sets of nodes in the corresponding layers.
            node_counts (list[list[int]]):             A list of lists, where each inner list contains the number of nodes at each 
                                                       layer in the corresponding subgraph.
            lp (int):                                  The length of the longest path identified in the original graph.
            num_graph (int):                           The total number of subgraphs generated after processing.
            graph_list (list[nx.DiGraph]):             A list of directed subgraphs, each representing weakly connected components 
                                                       of the pruned original graph.
    """
    global lp
    if len(edgelist) != len(set(edgelist)):
        raise ValueError("Error: The graph has duplicate edges.")
    G = nx.DiGraph(edgelist)
    if not nx.is_directed_acyclic_graph(G):
        raise ValueError("Error: The graph must be a DAG.")

    lp = nx.dag_longest_path_length(G)
    subgraph_dict, node_counts, graph_list, num_graph = [], [], [], 0

    if lp == 0:
        if iso_nodes is not None:
            G.add_nodes_from(iso_nodes)
            num_graph = G.number_of_nodes() - 1
        return subgraph_dict, node_counts, lp, num_graph, graph_list
    else:
        while True:
            G = prune(G)
            n1 = G.number_of_edges()
            G = lp_edgelist(G)
            n2 = G.number_of_edges()
            if n1 == n2:
                break
        
        wcc = [G.subgraph(c) for c in nx.weakly_connected_components(G)]
        for subgraph in wcc:
            layers = layer_dict(subgraph)
            num_graph += 1
            s_dict = {K: layers[K] for K in sorted(layers)}
            node_count = [len(layers[i]) for i in range(len(layers))]
               
            density_f, n_f = edge_density(subgraph, layers[0], layers[1])
            density_b, n_b = edge_density(subgraph, layers[len(layers) - 2], layers[len(layers) - 1])   
            if (density_f > density_b) or (abs(density_f - density_b) < 1e-6 and n_f > n_b):
                subgraph = nx.reverse(subgraph, copy = True)
                s_dict = {new_k: s_dict[old_k] for new_k, old_k in enumerate(sorted(s_dict.keys(), reverse = True))}
                node_count.reverse()

            graph_list.append(subgraph)
            subgraph_dict.append(s_dict)
            node_counts.append(node_count)
        
        return subgraph_dict, node_counts, lp, num_graph, graph_list

# An example usage
if __name__ == "__main__":
    edgelist = [('a0', 'b2'), ('a0', 'b3'), ('a1', 'b2'), ('a1', 'b3'), ('b4', 'd1'), ('a1', 'c1'), ('b0', 'c0'), ('b0', 'c1'), 
                ('b1', 'c0'), ('b1', 'c1'), ('b2', 'c2'), ('b2', 'c3'), ('b3', 'c2'), ('b3', 'c3'), ('b4', 'c4'), ('b4', 'c5'), 
                ('b5', 'c4'), ('b5', 'c5'), ('b0', 'c2'), ('b0', 'c3'), ('b1', 'c2'), ('b1', 'c3'), ('b4', 'c2'), ('b4', 'c3'), 
                ('b5', 'c2'), ('b5', 'c3'), ('c0', 'd0'), ('c0', 'd1'), ('c1', 'd0'), ('c1', 'd1'), ('c4', 'd2'), ('c4', 'd3'), 
                ('c5', 'd2'), ('c5', 'd3'), ('c6', 'd2'), ('c6', 'd3'), ('c7', 'd2'), ('c7', 'd3')]

    G = nx.DiGraph(edgelist)
    lp = nx.dag_longest_path_length(G)
    print(f"lp:{lp}\n")
    after_prune = prune(G)
    print(f"after prune:{after_prune.edges()},\n number of edges:{after_prune.number_of_edges()},\n new_lp:{nx.dag_longest_path_length(after_prune)}\n")
    after_lp = lp_edgelist(after_prune)
    print(f"lp_edgelist:{after_lp.edges()},\n number of edges:{after_lp.number_of_edges()}\n")
    wcc = [G.subgraph(c) for c in nx.weakly_connected_components(after_lp)]
    print(f"wcc:{wcc}\n")
    layer_dicts = [layer_dict(subgraph) for subgraph in wcc]
    print(f"layer_dicts:{layer_dicts}\n")
    print('--'*50)
    subgrph_dict, node_counts, lp, num_graph, graph_list = dag_process(edgelist)
    print(f"subgrph_dict:{subgrph_dict},\n node_counts:{node_counts},\n lp:{lp},\n num_graph:{num_graph},\n graph_list:{graph_list}\n")