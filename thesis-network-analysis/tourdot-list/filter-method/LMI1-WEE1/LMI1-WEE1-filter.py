#!/usr/bin/env python3

# File paths
links_file = '../3702-full-UniProt-GN-Name.tsv'
nodes_file = 'LMI1-WEE1-nodes.txt'
output_file = 'LMI1-WEE1-filtered-network.tsv'

print("Loading nodes of interest from LMI1-WEE1-nodes.txt...")
# Read the nodes of interest
nodes_of_interest = set()
with open(nodes_file, 'r') as f:
    for line in f:
        node = line.strip()
        if node and not node.startswith('#'):  # Skip comments and empty lines
            nodes_of_interest.add(node)

print(f"Loaded {len(nodes_of_interest)} nodes of interest")

print(f"Filtering edges from {links_file}...")
# Filter the network to include only rows with nodes of interest
filtered_edges = []
header = None
rows_checked = 0

with open(links_file, 'r') as f:
    # Skip the header (which may span multiple lines)
    header = f.readline().strip()
    # Continue reading until we find a line that looks like data (starts with a protein name)
    while header and not header[0].isupper():
        header = f.readline().strip()
    
    for line in f:
        parts = line.strip().split()
        if len(parts) >= 2:
            protein1 = parts[0]
            protein2 = parts[1]
            
            # Keep row if either protein is in nodes of interest
            if protein1 in nodes_of_interest or protein2 in nodes_of_interest:
                filtered_edges.append(line.strip())
            
            rows_checked += 1
            if rows_checked % 1000000 == 0:
                print(f"  Checked {rows_checked} rows, found {len(filtered_edges)} matching edges...")

print(f"Total rows checked: {rows_checked}")
print(f"Total matching edges found: {len(filtered_edges)}")

# Write filtered edges to output file
print(f"Writing filtered network to {output_file}...")
with open(output_file, 'w') as f:
    f.write(header + '\n')
    for edge in filtered_edges:
        f.write(edge + '\n')

print(f"✓ Filtered network saved to {output_file}")
print(f"\nNext step: Run LMI1-WEE1-visualize.py to visualize the filtered network")
