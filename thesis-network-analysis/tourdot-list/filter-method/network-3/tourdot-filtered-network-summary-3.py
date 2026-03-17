#!/usr/bin/env python3

from collections import defaultdict
from pathlib import Path

# File paths
network_file = 'tourdot-filtered-network-3.tsv'
output_file = 'tourdot-filtered-network-summary.tsv'

print("Counting edges for each protein...")
edge_count = defaultdict(int)

with open(network_file, 'r') as f:
    header = f.readline()
    for line in f:
        parts = line.strip().split()
        if len(parts) >= 2:
            protein1 = parts[0]
            protein2 = parts[1]
            edge_count[protein1] += 1
            edge_count[protein2] += 1

print(f"Found {len(edge_count)} unique proteins")

# Sort by edge count (highest to lowest)
sorted_proteins = sorted(edge_count.items(), key=lambda x: x[1], reverse=True)

# Write to output file
with open(output_file, 'w') as f:
    f.write("protein\tedges\n")
    for protein, count in sorted_proteins:
        f.write(f"{protein}\t{count}\n")

print(f"✓ Summary saved to {output_file}")

# Print summary statistics
print(f"\nSummary Statistics:")
print(f"  Total unique proteins: {len(edge_count)}")
print(f"  Max edges: {max(edge_count.values())}")
print(f"  Min edges: {min(edge_count.values())}")
print(f"  Average edges: {sum(edge_count.values()) / len(edge_count):.2f}")
