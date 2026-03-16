#!/usr/bin/env python3

import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from pathlib import Path

# File paths
network_file = 'LMI1-WEE1-filtered-network.tsv'
nodes_file = 'LMI1-WEE1-nodes.txt'
output_file = 'LMI1-WEE1-network.png'

print("Loading nodes of interest...")
nodes_of_interest = set()
with open(nodes_file, 'r') as f:
    for line in f:
        node = line.strip()
        if node and not node.startswith('#'):
            nodes_of_interest.add(node)

print(f"Loaded {len(nodes_of_interest)} nodes of interest")

print("Building network from filtered edges...")
G = nx.Graph()

with open(network_file, 'r') as f:
    header = f.readline()
    line_count = 0
    for line in f:
        parts = line.strip().split()
        if len(parts) >= 3:
            protein1 = parts[0]
            protein2 = parts[1]
            try:
                score = int(parts[-1]) if parts[-1].isdigit() else 500
            except:
                score = 500
            
            G.add_edge(protein1, protein2, weight=score)
            line_count += 1

print(f"Loaded {G.number_of_nodes()} nodes and {G.number_of_edges()} edges")

# Create layout
print("Computing network layout (this may take a moment)...")
pos = nx.spring_layout(G, k=2, iterations=50, seed=42)

# Create figure
fig, ax = plt.subplots(figsize=(20, 20), dpi=150)

# Separate nodes into categories
nodes_of_interest_in_graph = [node for node in G.nodes() if node in nodes_of_interest]
other_nodes = [node for node in G.nodes() if node not in nodes_of_interest]

# Draw edges
nx.draw_networkx_edges(G, pos, ax=ax, edge_color='gray', width=0.5, alpha=0.3)

# Draw nodes
nx.draw_networkx_nodes(G, pos, nodelist=other_nodes, node_color='lightblue', 
                       node_size=300, ax=ax, alpha=0.7, label='Other proteins')

nx.draw_networkx_nodes(G, pos, nodelist=nodes_of_interest_in_graph, node_color='red', 
                       node_size=800, ax=ax, label='Nodes of interest')

# Draw labels
nx.draw_networkx_labels(G, pos, font_size=8, font_weight='bold', ax=ax)

# Create legend
legend_elements = [
    mpatches.Patch(facecolor='red', edgecolor='black', label='Nodes of interest'),
    mpatches.Patch(facecolor='lightblue', edgecolor='black', label='Other proteins')
]
ax.legend(handles=legend_elements, loc='upper left', fontsize=12)

ax.set_title(f'LMI1-WEE1 Filtered Network\n({G.number_of_nodes()} nodes, {G.number_of_edges()} edges)', 
             fontsize=16, fontweight='bold')
ax.axis('off')

print(f"Saving visualization to {output_file}...")
plt.tight_layout()
plt.savefig(output_file, dpi=150, bbox_inches='tight', facecolor='white')
print(f"✓ Visualization saved to {output_file}")

# Print summary
print(f"\nNetwork Summary:")
print(f"  Total nodes: {G.number_of_nodes()}")
print(f"  Total edges: {G.number_of_edges()}")
print(f"  Nodes of interest in network: {len(nodes_of_interest_in_graph)}")
print(f"  Intermediate nodes: {len(other_nodes)}")
print(f"  Network density: {nx.density(G):.4f}")

if nx.is_connected(G):
    print(f"  ✓ Network is fully connected")
else:
    components = list(nx.connected_components(G))
    print(f"  ✗ Network has {len(components)} connected components")
