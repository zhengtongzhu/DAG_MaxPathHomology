# MaxPathHomology
Zhu, Z. and Chi, Z. (2024). "*Recursive Computation of Path Homology for Stratified Digraphs*".

This Python script is an implementation of the above paper that recursively computes the maximal (reduced) path homology of an unweighted stratified graph or unweighted directed acyclic graph (DAG).

## Installation
To get started with this project, follow the steps below to install the necessary dependencies:

- **Python**: Ensure that Python 3.8 or higher is installed. You can check your Python version by running:
    ```bash
    python --version
    ```
- **Clone the repository**:
    ```bash
    git clone https://github.com/zhengtongzhu/DAG_MaxPathHomology.git
    cd DAG_MaxPathHomology
    ```
- **Install dependencies**:
    ```bash
    pip install -r requirements.txt
    ```
- **Run the project**:
    ```bash
    python recursive_algorithm.py
    ```

    You will see the following output from the test example:
    ```python
    lp: 2,
    dim: 3,
    basis: [
        [(-b0 + b1)*((-c0 + c1)*(-c3 + d1) + (c0 - c1)*(-c3 + d0)), 
        (-a0 + a1)*(-b2 + b3)*(c2 - c3)], 
        [(-b4 + b5)*(-c4 + c5)*(d2 - d3)]
        ]
    ```

## dag_preprocess.py
The `graph_preprocess.py` module provides functions to compute `l(G)` and the stratified decomposition $G_*$ for a given unweighted DAG without multiple edges. These methods are described in Section 5.2 of the paper.

### Example Usage
Consider a DAG `G` with the following `edgelist`:

```python
edgelist = [('a0', 'b2'), ('a0', 'b3'), ('a1', 'b2'), ('a1', 'b3'), ('b4', 'd1'), 
            ('a1', 'c1'), ('b0', 'c0'), ('b0', 'c1'), ('b1', 'c0'), ('b1', 'c1'), 
            ('b2', 'c2'), ('b2', 'c3'), ('b3', 'c2'), ('b3', 'c3'), ('b4', 'c4'), 
            ('b4', 'c5'), ('b5', 'c4'), ('b5', 'c5'), ('b0', 'c2'), ('b0', 'c3'), 
            ('b1', 'c2'), ('b1', 'c3'), ('b4', 'c2'), ('b4', 'c3'), ('b5', 'c2'), 
            ('b5', 'c3'), ('c0', 'd0'), ('c0', 'd1'), ('c1', 'd0'), ('c1', 'd1'), 
            ('c4', 'd2'), ('c4', 'd3'), ('c5', 'd2'), ('c5', 'd3'), ('a2', 'b4')]
```

each tuple `(a, b)` in this edgelist represents a directed unweighted edge in `G` from node `a` to node `b`. Here's a visualization of `G`:

<img src="figures\example_G.pdf" width="70%" style="display: block; margin: auto;" />

The `dag_process` function first check if the DAG `G` contains multi-edges or has a loop (based on [NetworkX](https://networkx.org/)):

```python
if len(edgelist) != len(set(edgelist)):
    raise ValueError("Error: The graph has duplicate edges.")
G = nx.DiGraph(edgelist)
if not nx.is_directed_acyclic_graph(G):
    raise ValueError("Error: The graph must be a DAG.")
```

then compute the longest path length of `G` (by [dag_longest_path_length](https://networkx.org/documentation/stable/reference/algorithms/generated/networkx.algorithms.dag.dag_longest_path_length.html#networkx.algorithms.dag.dag_longest_path_length)):

```python
lp = nx.dag_longest_path_length(G)
```

We repeatedly prune the graph using `prune` and $G_*$-algorithm (`lp_edgelist`) until the graph structure can no longer be simplified. The repetition of the `prune` and `lp_edgelist` provides a set of weakly connected components `{G_i}`. These `{G_i}` are stratified graphs, and computing the direct sum of the maximal path homology of all `G_i` is equivalent to computing the maximal path homology of `G`, which are guaranteed by Corollary 3.5 and Proposition 3.7 of the paper. For the `edgelist` above, `10` edges are removed after pruning:

- `('a2', 'b4')`: removed by `prune`.
- `('b0', 'c2')`, `('b0', 'c3')`, `('b1', 'c2')`, `('b1', 'c3')`, `('b4', 'c2')`, `('b4', 'c3')`, `('b5', 'c2')`, `('b5', 'c3')`, `('b4', 'd1')`: removed by `lp_edgelist`.

The original graph `G` is splitted into two weakly connected components:
```python
G_0: [('a0', 'b2'), ('a0', 'b3'), ('a1', 'b2'), ('a1', 'b3'), ('a1', 'c1'),
      ('b2', 'c2'), ('b2', 'c3'), ('b3', 'c2'), ('b3', 'c3'),
      ('b0', 'c0'), ('b0', 'c1'), ('b1', 'c0'), ('b1', 'c1'),
      ('c0', 'd0'), ('c0', 'd1'), ('c1', 'd0'), ('c1', 'd1')]

G_1: [('b4', 'c4'), ('b4', 'c5'), ('b5', 'c4'), ('b5', 'c5'),
      ('c4', 'd2'), ('c4', 'd3'), ('c5', 'd2'), ('c5', 'd3')]
```

The `dag_process` returns 5 components of all `G_i`:

```python
subgraph_dict, node_counts, lp, num_graph, graph_list = dag_process(edgelist)
```

- The $i^{th}$ element of `subgraph_dict` is a dictionary represents `G_i`, where each key is a layer index (start from `0`) and the value corresponding to the key is the set of nodes in that layer. For example, in the `edgelist` above, the `subgraph_dict` is:

```python
subgraph_dict = [
    {0: {'b1', 'a1', 'b0', 'a0'}, 1: {'c1', 'b2', 'b3', 'c0'}, 2: {'c2', 'c3', 'd1', 'd0'}},
    {0: {'b5', 'b4'}, 1: {'c5', 'c4'}, 2: {'d3', 'd2'}}
    ]
```
`c2` is a node in `G_0`'s `2`nd layer and `c4` is a node in `G_1`'s `1`st layer.

- The `i`$^{th}$ element of `node_counts` represents the number of nodes in different layer of `G_i`. From the `subgraph_dict` above we know that

```python
node_counts = [
    [4, 4, 4], 
    [2, 2, 2]
    ]
```
where `[4, 4, 4]` means `G_0` has `4` nodes in layer `0`, `4` nodes in layer `1`, and `4` nodes in layer `2`.

- `lp` is the longest path length of `G` as described above (`lp = l(G)`). Here `lp = 2`.
- `num_graph` is the cardinality of `{G_i}`. Here `num_graph = 2`.
- `graph_list` is a list where each element `G_i` is stored as an instance of the `nx.DiGraph` class.

## recursive_algrithm.py
The `recursive_algrithm.py` includes the main algorithm `max_path_homology` to compute the maximal path homology for an unweighted stratified graph. Along with `graph_preprocess.py`, it also works for unweighted DAGs.

The `max_path_homology` returns `lp`, `sum_dim` and `basis`.
```python
lp, dim, basis = max_path_homology(edgelist, calculate_basis)
```
Here:
- `lp=l(G)` is the longest path length of `G`.
- `dim` is the sum of the Betti Numbers of the `lp`-dimentional (maximal) path homologies of all `G_i`.
- If the input `calculate_basis == True`, then `basis` returns a basis of the `lp`-dimensional (maximal) reduced path homology of all `G_i`, otherwise it returns `None`.

For the `edgelist` above, if `calculate_basis == True`:
```python
lp = 2
sum_dim = 3
basis = [
    [(-b0 + b1)*((-c0 + c1)*(-c3 + d1) + (c0 - c1)*(-c3 + d0)), 
    (-a0 + a1)*(-b2 + b3)*(c2 - c3)], 
    [(-b4 + b5)*(-c4 + c5)*(d2 - d3)]
    ]
```
## general_algorithm.py
The `general_algorithm.py` is based on an implementation from the following project:

Carranza, D., Doherty, B., Kapulkin, K., Opie, M., Sarazola, M., & Wong, L. Z. (2022). *Python script for computing path homology of digraphs* (Version 1.0.0) [Computer software]. https://github.com/sheaves/path_homology.

This script implements a general algorithm for computing reduced path homology. If you use this project or its comparative implementation in your work, please also cite the above project to acknowledge their contributions.

## Other Modules
`maxpph.py`: Produces a decreasing persistence path homology plot (Section 6.2 of the paper).

`experiment_func.py`, `stratified_gamma_1_2_3.py` and `stratified_gamma_4_5.py`: Implement additional functions and simulations discussed in Section 6.1 of the paper.

`maxph_matrix.py`: Contains functions for matrix operations.

## Citation
If you find this code useful, please cite it using the following BibTeX entry:

```bibtex
@software{Zhu_Computing_the_maximal_2024,
author = {Zhu, Zhengtong and Chi, Zhiyi},
month = dec,
title = {{Computing the maximal path homology of directed acyclic graph}},
url = {https://github.com/zhengtongzhu/DAG_MaxPathHomology},
version = {1.0.0},
year = {2024}
}
```