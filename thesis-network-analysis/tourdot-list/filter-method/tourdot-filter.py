#!/usr/bin/env python3

import networkx as nx
from collections import defaultdict, deque
import sys

# File paths
links_file = '3702-full-UniProt-GN-Name.tsv'
nodes_file = 'tourdot-notes-filter-revised.txt'
output_file = 'tourdot-filtered-network-revised.tsv'

print("Loading nodes of interest...")
# Read the nodes of interest
nodes_of_interest = set()
with open(nodes_file, 'r') as f:
    for line in f:
        node = line.strip()
        if node and not node.startswith('#'):  # Skip comments and empty lines
            nodes_of_interest.add(node)

print(f"Loaded {len(nodes_of_interest)} nodes of interest")

print("Building network graph...")
# Build network graph
G = nx.Graph()
edges_loaded = 0
edge_lines = {}  # Store line data for each edge

with open(links_file, 'r') as f:
    header = f.readline()  # Skip header
    for line in f:
        parts = line.strip().split()
        if len(parts) >= 3:
            protein1 = parts[0]
            protein2 = parts[1]
            edge_key = tuple(sorted([protein1, protein2]))
            
            # Add edge with score (combined_score is last column)
            try:
                score = int(parts[-1]) if parts[-1].isdigit() else 0
                weight = 1000 - score  # Invert score so higher score = lower weight
            except:
                weight = 500
            
            if not G.has_edge(protein1, protein2):
                G.add_edge(protein1, protein2, weight=weight)
                edge_lines[edge_key] = line.strip()
            
            edges_loaded += 1
            if edges_loaded % 1000000 == 0:
                print(f"Loaded {edges_loaded} edges...")

print(f"Total edges loaded: {edges_loaded}")
print(f"Total nodes in graph: {G.number_of_nodes()}")

# Find which nodes of interest are in the graph
nodes_in_graph = nodes_of_interest & set(G.nodes())
nodes_missing = nodes_of_interest - nodes_in_graph

print(f"Nodes of interest in graph: {len(nodes_in_graph)}")
print(f"Nodes of interest NOT in graph: {len(nodes_missing)}")

if nodes_missing:
    print(f"Missing nodes: {sorted(nodes_missing)}")

# Build minimal network connecting all nodes of interest using shortest paths
print("\nFinding shortest paths to connect nodes...")

# Start with the first node, then iteratively add closest nodes
selected_edges = set()
selected_nodes = set()

if nodes_in_graph:
    # Start with first node
    start_node = min(sorted(nodes_in_graph))  # Use sorted for consistency
    selected_nodes.add(start_node)
    remaining = nodes_in_graph - {start_node}
    
    iteration = 0
    max_iterations = len(nodes_in_graph)
    
    while remaining and iteration < max_iterations:
        iteration += 1
        
        # Find closest node to the current set
        closest_node = None
        shortest_path_length = float('inf')
        shortest_path = None
        
        # Use BFS from the entire selected set for efficiency
        for target in remaining:
            # Find shortest distance from any selected node to target
            best_source = None
            best_path = None
            best_length = float('inf')
            
            for source in list(selected_nodes)[:5]:  # Check top 5 selected nodes for speed
                try:
                    path = nx.shortest_path(G, source, target, weight='weight')
                    path_length = len(path) - 1
                    
                    if path_length < best_length:
                        best_length = path_length
                        best_source = source
                        best_path = path
                except (nx.NetworkXNoPath, nx.NetworkXError):
                    continue
            
            if best_path and best_length < shortest_path_length:
                shortest_path_length = best_length
                closest_node = target
                shortest_path = best_path
        
        if closest_node and shortest_path:
            # Add all edges from the shortest path
            for i in range(len(shortest_path) - 1):
                edge = tuple(sorted([shortest_path[i], shortest_path[i+1]]))
                selected_edges.add(edge)
            
            selected_nodes.add(closest_node)
            remaining.remove(closest_node)
            
            print(f"Iteration {iteration}: Connected {closest_node} with path length {shortest_path_length}")
        else:
            if remaining:
                print(f"Could not connect {len(remaining)} remaining nodes: {sorted(remaining)[:5]}...")
            break

print(f"\nSelected {len(selected_edges)} edges")
print(f"Selected {len(selected_nodes)} nodes from nodes of interest")

# Create output network with selected edges and their data
print(f"Writing output to {output_file}...")

# Write output
with open(output_file, 'w') as f:
    f.write(header)
    for edge in sorted(selected_edges):
        if edge in edge_lines:
            f.write(edge_lines[edge] + '\n')

print(f"Output written to {output_file}")

# Verify connectivity
print("\nVerifying network connectivity...")
subgraph_nodes = set()
for edge in selected_edges:
    subgraph_nodes.add(edge[0])
    subgraph_nodes.add(edge[1])

subgraph = G.subgraph(subgraph_nodes)

# Check if the subgraph is connected
if nx.is_connected(subgraph):
    print("✓ Network is connected (no islands)")
else:
    components = list(nx.connected_components(subgraph))
    print(f"✗ Network has {len(components)} connected components (islands detected)")
    for i, comp in enumerate(components):
        # Show nodes of interest in each component
        comp_nodes_of_interest = nodes_of_interest & comp
        print(f"  Component {i+1}: {len(comp)} total nodes, {len(comp_nodes_of_interest)} nodes of interest")
        if len(comp_nodes_of_interest) <= 5:
            print(f"    Nodes of interest: {sorted(comp_nodes_of_interest)}")

# Check nodes of interest connectivity
nodes_of_interest_in_subgraph = nodes_of_interest & subgraph_nodes
print(f"\nNodes of interest in output network: {len(nodes_of_interest_in_subgraph)}/{len(nodes_in_graph)}")

if len(nodes_of_interest_in_subgraph) < len(nodes_in_graph):
    missing_from_output = nodes_in_graph - nodes_of_interest_in_subgraph
    print(f"Nodes of interest NOT in output: {sorted(missing_from_output)}")