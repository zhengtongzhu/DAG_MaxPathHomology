import networkx as nx
from collections import deque

def layer_dict(G: nx.DiGraph) -> dict[int, set[str]]:
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

def dag_process(edgelist: list[tuple[str, str]], iso_nodes: list[str] | None = None) -> tuple[list[dict[int, set[str]]], list[list[int]], int, int, list[nx.DiGraph]]:
    """
    This function takes an edgelist representing a DAG without multi-edges, dividing it into subgraphs 
    that are equivalent for calculating maximal path homology.

    Args:
        edgelist (list[tuple[str, str]]): A list of tuples where each tuple represents a directed edge in the graph.
        iso_nodes (list[str]):            A list of isolated nodes of the graph.

    Returns:
        Tuple:
            subgraph_dict (list[dict[int, set[str]]]): A list of dictionaries where each dictionary represents the nodes grouped by 
                                                       their layer number within a subgraph. The keys are the layer numbers and 
                                                       the values are sets of nodes in the corresponding layers.
            node_counts (list[list[int]]):             A list of lists, where each inner list contains the number of nodes at each 
                                                       layer in the corresponding subgraph.
            lp (int):                                  The length of the longest path identified in the original graph.
            num_graph (int):                           The total number of stratified weakly connected subgraphs generated after processing.
            graph_list (list[nx.DiGraph]):             A list of directed subgraphs, each representing weakly connected components 
                                                       of the pruned original graph.
    """
    if len(edgelist) != len(set(edgelist)):
        raise ValueError("Error: The graph has multi-edges.")
    
    G = nx.DiGraph(edgelist)

    if iso_nodes is not None:
        G.add_nodes_from(iso_nodes)
    
    if not nx.is_directed_acyclic_graph(G):
        raise ValueError("Error: The graph must be a DAG.")

    lp = nx.dag_longest_path_length(G)
    subgraph_dict, node_counts, graph_list, num_graph = [], [], [], 0

    if lp == 0:
        if iso_nodes is not None:
            num_graph = G.number_of_nodes() - 1
        return subgraph_dict, node_counts, lp, num_graph, graph_list
    else:
        G_decompose = graph_decomposition(G)
        G_prune = prune(G_decompose)

        if G_prune.number_of_nodes() == 0:
            lp_new = 0
        else:
            lp_new = nx.dag_longest_path_length(G_prune)

        if lp != lp_new:
            print(f"The {lp}-th order path homology of the original graph is trivial.")
            return subgraph_dict, node_counts, lp, num_graph, graph_list
        else:
            wcc = [G_prune.subgraph(c).copy() for c in nx.weakly_connected_components(G_prune)]
            for subgraph in wcc:
                layers = layer_dict(subgraph)
                num_graph += 1
                s_dict = {K: layers[K] for K in sorted(layers)}
                node_count = [len(layers[i]) for i in range(len(layers))]
                graph_list.append(subgraph)
                subgraph_dict.append(s_dict)
                node_counts.append(node_count)
        
            return subgraph_dict, node_counts, lp, num_graph, graph_list

# An example usage
if __name__ == "__main__":
    edgelist = [('a0', 'b2'), ('a0', 'b3'), ('a1', 'b2'), ('a1', 'b3'), ('b4', 'd1'), 
                ('a1', 'c1'), ('b0', 'c0'), ('b0', 'c1'), ('b1', 'c0'), ('b1', 'c1'), 
                ('b2', 'c2'), ('b2', 'c3'), ('b3', 'c2'), ('b3', 'c3'), ('b4', 'c4'), 
                ('b4', 'c5'), ('b5', 'c4'), ('b5', 'c5'), ('b0', 'c2'), ('b0', 'c3'), 
                ('b1', 'c2'), ('b1', 'c3'), ('b4', 'c2'), ('b4', 'c3'), ('b5', 'c2'), 
                ('b5', 'c3'), ('c0', 'd0'), ('c0', 'd1'), ('c1', 'd0'), ('c1', 'd1'), 
                ('c4', 'd2'), ('c4', 'd3'), ('c5', 'd2'), ('c5', 'd3'), ('a2', 'b4')]

    G = nx.DiGraph(edgelist)
    lp = nx.dag_longest_path_length(G)

    print(f"original graph:{G.edges()},\n number of edges = {G.number_of_edges()},\n lp:{lp}\n")
    subgrph_dict, node_counts, lp, num_graph, graph_list = dag_process(edgelist)
    print(f"subgrph_dict:{subgrph_dict},\n node_counts:{node_counts},\n lp:{lp},\n num_graph:{num_graph},\n graph_list:{graph_list}\n")