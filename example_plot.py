import networkx as nx
import matplotlib.pyplot as plt

# create graph
G = nx.DiGraph()
edgelist = [
    ('a0', 'b2', 'blue'), ('a0', 'b3', 'blue'), ('a1', 'b2', 'blue'),
    ('a1', 'b3', 'blue'), ('b4', 'd1', 'red'), ('a1', 'c1', 'blue'),
    ('b0', 'c0', 'blue'), ('b0', 'c1', 'blue'), ('b1', 'c0', 'blue'),
    ('b1', 'c1', 'blue'), ('b2', 'c2', 'green'), ('b2', 'c3', 'green'),
    ('b3', 'c2', 'green'), ('b3', 'c3', 'green'), ('b4', 'c4', 'green'),
    ('b4', 'c5', 'green'), ('b5', 'c4', 'blue'), ('b5', 'c5', 'blue'),
    ('b0', 'c2', 'red'), ('b0', 'c3', 'red'), ('b1', 'c2', 'red'),
    ('b1', 'c3', 'red'), ('b4', 'c2', 'red'), ('b4', 'c3', 'red'),
    ('b5', 'c2', 'red'), ('b5', 'c3', 'red'), ('c0', 'd0', 'green'),
    ('c0', 'd1', 'green'), ('c1', 'd0', 'green'), ('c1', 'd1', 'green'),
    ('c4', 'd2', 'green'), ('c4', 'd3', 'green'), ('c5', 'd2', 'green'),
    ('c5', 'd3', 'green'), ('a2', 'b4', 'red')
]

for src, dst, color in edgelist:
    G.add_edge(src, dst, color = color)

# label nodes
G_first_row = ['a0', 'a1', 'a2']
G_second_row = ['b0', 'b1', 'b2', 'b3', 'b4', 'b5']
G_third_row = ['c0', 'c1', 'c2', 'c3', 'c4', 'c5']
G_fourth_row = ['d0', 'd1', 'd2', 'd3']

# set position
G_pos = {}
for i, node in enumerate(G_first_row):
    G_pos[node] = (i * 2, 4)
for i, node in enumerate(G_second_row):
    G_pos[node] = (i * 2, 3)
for i, node in enumerate(G_third_row):
    G_pos[node] = (i * 2, 2)
for i, node in enumerate(G_fourth_row):
    G_pos[node] = (i * 2, 1)

# edge color
G_edge_colors = [G[src][dst]['color'] for src, dst in G.edges()]

# plot
plt.figure(figsize=(18, 12))
nx.draw(
    G,
    G_pos,
    with_labels = True,
    edge_color = G_edge_colors,
    arrows = True,
    node_size = 800,
    font_size = 12,
    connectionstyle = "arc3,rad=0.2",
    arrowstyle = "->",
)
plt.title("Original Graph G", fontsize = 16, pad = 20)
plt.tight_layout()
#plt.savefig("example_G.pdf")
plt.show()
plt.close()

##G_0
G_0 = nx.DiGraph()
edgelist_0 = [('a0', 'b2', 'blue'), ('a0', 'b3', 'blue'), ('a1', 'b2', 'blue'),
    ('a1', 'b3', 'blue'), ('a1', 'c1', 'blue'),
    ('b0', 'c0', 'blue'), ('b0', 'c1', 'blue'), ('b1', 'c0', 'blue'),
    ('b1', 'c1', 'blue'), ('b2', 'c2', 'green'), ('b2', 'c3', 'green'),
    ('b3', 'c2', 'green'), ('b3', 'c3', 'green'), 
    ('c0', 'd0', 'green'),
    ('c0', 'd1', 'green'), ('c1', 'd0', 'green'), ('c1', 'd1', 'green')
]

for src, dst, color in edgelist_0:
    G_0.add_edge(src, dst, color = color)

# label nodes
G_0_first_row = ['a0', 'a1', 'b0', 'b1']
G_0_second_row = ['c0', 'c1', 'b2', 'b3']
G_0_third_row = ['d0', 'd1', 'c2', 'c3']

# set position
G_0_pos = {}
for i, node in enumerate(G_0_first_row):
    G_0_pos[node] = (i * 2, 3)
for i, node in enumerate(G_0_second_row):
    G_0_pos[node] = (i * 2, 2)
for i, node in enumerate(G_0_third_row):
    G_0_pos[node] = (i * 2, 1)

# edge color
G_0_edge_colors = [G_0[src][dst]['color'] for src, dst in G_0.edges()]

# plot
plt.figure(figsize=(18, 12))
nx.draw(
    G_0,
    G_0_pos,
    with_labels = True,
    edge_color = G_0_edge_colors,
    arrows = True,
    node_size = 800,
    font_size = 12,
    arrowstyle = "->",
)
plt.title("G_0", fontsize = 16, pad = 20)
plt.tight_layout()
#plt.savefig("G_0.pdf")
plt.show()
plt.close()

##G_1
G_1 = nx.DiGraph()
edgelist_1 = [
    ('c4', 'd2', 'green'), ('c4', 'd3', 'green'), ('c5', 'd2', 'green'),
    ('c5', 'd3', 'green'), ('b4', 'c4', 'blue'),
    ('b4', 'c5', 'blue'), ('b5', 'c4', 'blue'), ('b5', 'c5', 'blue')
]

for src, dst, color in edgelist_1:
    G_1.add_edge(src, dst, color = color)

# label nodes
G_1_first_row = ['b4', 'b5']
G_1_second_row = ['c4', 'c5']
G_1_third_row = ['d2', 'd3']

# set position
G_1_pos = {}
for i, node in enumerate(G_1_first_row):
    G_1_pos[node] = (i * 2, 3)
for i, node in enumerate(G_1_second_row):
    G_1_pos[node] = (i * 2, 2)
for i, node in enumerate(G_1_third_row):
    G_1_pos[node] = (i * 2, 1)

# edge color
G_1_edge_colors = [G_1[src][dst]['color'] for src, dst in G_1.edges()]

# plot
plt.figure(figsize=(18, 12))
nx.draw(
    G_1,
    G_1_pos,
    with_labels = True,
    edge_color = G_1_edge_colors,
    arrows = True,
    node_size = 800,
    font_size = 12,
    arrowstyle = "->",
)
plt.title("G_1", fontsize = 16, pad = 20)
plt.tight_layout()
#plt.savefig("G_1.pdf")
plt.show()
plt.close()