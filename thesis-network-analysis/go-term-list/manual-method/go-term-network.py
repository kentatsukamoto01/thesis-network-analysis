import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np

# List of TSV files to include in the network
# Add or remove files from this list to include/exclude them in the visualization
tsv_files = ['DPB.tsv', 'CDKA-1.tsv', 'E2FC.tsv', 'CDC6.tsv', 'CYCD5-1.tsv', 'CYCA2-1.tsv', 'CDKB1-1.tsv', 'UPL3.tsv', 'KRP2.tsv', 'CKS1.tsv', 'TCP15.tsv', 'E2FE.tsv', 'PYM.tsv', 'BIN4.tsv', 'FZR2.tsv', 'CYCA2-3.tsv', 'SIM.tsv']

# Column names for the TSV files
column_names = ['node1', 'node2', 'node1_string_id', 'node2_string_id', 'neighborhood_on_chromosome', 'gene_fusion', 'phylogenetic_cooccurrence', 'homology', 'coexpression', 'experimentally_determined_interaction', 'database_annotated', 'automated_textmining', 'combined_score']

# Load and combine all TSV files
all_data = []
for tsv_file in tsv_files:
    df = pd.read_csv(tsv_file, sep='\t', comment='#', header=None, names=column_names)
    all_data.append(df)

# Combine all dataframes
combined_df = pd.concat(all_data, ignore_index=True)

# Ensure unique undirected edges (display only one edge per pair, keeping data intact)
combined_df['min_node'] = combined_df[['node1', 'node2']].min(axis=1)
combined_df['max_node'] = combined_df[['node1', 'node2']].max(axis=1)
df_unique = combined_df.drop_duplicates(subset=['min_node', 'max_node'])

# Create a NetworkX graph
G = nx.Graph()

# Add edges with weights from unique pairs
for index, row in df_unique.iterrows():
    G.add_edge(row['min_node'], row['max_node'], weight=row['combined_score'])

# Draw the network
plt.figure(figsize=(12, 12))
pos = nx.spring_layout(G, k=0.2)  # Layout for positioning nodes
edges = G.edges()
weights = [G[u][v]['weight'] for u, v in edges]

# Draw nodes
nx.draw_networkx_nodes(G, pos, node_size=200, node_color='lightblue')

# Draw edges with thinner varying thickness based on weight
nx.draw_networkx_edges(G, pos, width=[w*3 for w in weights], edge_color='gray')

# Draw labels
nx.draw_networkx_labels(G, pos, font_size=6)

# Add title
plt.title('GO Term Labeled DNA Endoreduplication Network')

# Save as PNG
plt.savefig('go-term-network.png', dpi=300, bbox_inches='tight')

print("Network visualization saved as 'go-term-network.png'")

# Create a TSV file with node and edge count
node_degrees = dict(G.degree())
nodes_df = pd.DataFrame(list(node_degrees.items()), columns=['node', 'edge_count'])
nodes_df = nodes_df.sort_values('edge_count', ascending=False)

# Save to TSV
nodes_df.to_csv('go-term-network-list.tsv', sep='\t', index=False)

print("Node edge count list saved as 'go-term-network-list.tsv'")