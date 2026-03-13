import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np

# List of TSV files to include in the network
# Add or remove files from this list to include/exclude them in the visualization
tsv_files = ['BASL.tsv', 'CCS52A2.tsv', 'CDC20-1.tsv', 'CDC20-2.tsv', 'CDC20-3.tsv', 'CDC20-4.tsv', 'CDC20-5.tsv', 'CDC20-6.tsv', 'CDKA-1.tsv', 'CDKB1-2.tsv', 'CDKB1;1.tsv', 'CAS1.tsv', 'CYCA2;3.tsv', 'DEL1.tsv', 'E2Fa.tsv', 'E2Fb.tsv', 'E2Fc.tsv', 'EBP1.tsv', 'FAMA.tsv', 'FBL17.tsv', 'FZR2.tsv', 'GIG1.tsv', 'GTL1.tsv', 'KRP2.tsv', 'KRP3.tsv', 'KRP4.tsv', 'KRP7.tsv', 'MED16.tsv', 'MMS21.tsv', 'MYB124.tsv', 'MYB88.tsv', 'PBP1-2.tsv', 'PL1.tsv', 'RBR1.tsv', 'RECQSIM.tsv', 'RKP.tsv', 'SCL28-2.tsv', 'SCL28.tsv', 'SIM.tsv', 'SKIP.tsv', 'SKP2A.tsv', 'SMR1.tsv', 'SMR2.tsv', 'SOG1.tsv', 'SPCH.tsv', 'UBC19.tsv', 'UBP14.tsv', 'UVI4.tsv', 'WEE1.tsv', 'ATPK1.tsv', 'E2Fd.tsv', 'PYM.tsv', 'ATXR6.tsv', 'CYCL1-1.tsv', 'CDKG1.tsv', 'RS2Z33.tsv', 'SMR5.tsv', 'KRP1.tsv', 'KRP1-2.tsv', 'SKP2B.tsv']

# Column names for the TSV files
column_names = ['node1', 'node2', 'node1_string_id', 'node2_string_id', 'neighborhood_on_chromosome', 'gene_fusion', 'phylogenetic_cooccurrence', 'homology', 'coexpression', 'experimentally_determined_interaction', 'database_annotated', 'automated_textmining', 'combined_score']

# Tsv files whose nodes should be colored red
red_node_files = ['SKP2A.tsv', 'MMS21.tsv', 'EBP1.tsv', 'SCL28-2.tsv', 'RECQSIM.tsv', 'SMR2.tsv', 'PBP1-2.tsv', 'RKP.tsv', 'MYB88.tsv', 'CDKB1-2.tsv', 'FAMA.tsv', 'MYB124.tsv', 'KRP3.tsv', 'KRP4.tsv', 'KRP7.tsv', 'GIG1.tsv', 'ATPK1.tsv', 'E2Fd.tsv', 'ATXR6.tsv', 'RS2Z33.tsv', 'SKP2B.tsv']

# Load and combine all TSV files, tracking nodes from red files
all_data = []
red_nodes = set()
for tsv_file in tsv_files:
    df = pd.read_csv(tsv_file, sep='\t', comment='#', header=None, names=column_names)
    all_data.append(df)
    
    # If this file is in the red_node_files list, add its nodes to red_nodes set
    if tsv_file in red_node_files:
        red_nodes.update(df['node1'].values)
        red_nodes.update(df['node2'].values)

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
pos = nx.spring_layout(G, k=0.1)  # Layout for positioning nodes
edges = G.edges()
weights = [G[u][v]['weight'] for u, v in edges]

# Create node colors based on red_nodes set
node_colors = ['red' if node in red_nodes else 'lightblue' for node in G.nodes()]

# Draw nodes
nx.draw_networkx_nodes(G, pos, node_size=200, node_color=node_colors)

# Draw edges with thinner varying thickness based on weight
nx.draw_networkx_edges(G, pos, width=[w*3 for w in weights], edge_color='gray')

# Draw labels
nx.draw_networkx_labels(G, pos, font_size=6)

# Add title
plt.title('Tourdot Endoreduplication Network')

# Save as PNG
plt.savefig('tourdot-network.png', dpi=300, bbox_inches='tight')

print("Network visualization saved as 'tourdot-network.png'")

# Create a TSV file with node and edge count
node_degrees = dict(G.degree())
nodes_df = pd.DataFrame(list(node_degrees.items()), columns=['node', 'edge_count'])
nodes_df = nodes_df.sort_values('edge_count', ascending=False)

# Save to TSV
nodes_df.to_csv('tourdot-network-list.tsv', sep='\t', index=False)

print("Node edge count list saved as 'tourdot-network-list.tsv'")