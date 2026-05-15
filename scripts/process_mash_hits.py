#!/usr/bin/env python3
"""
Process mash distance results to extract the 10 best hits (lowest distances).
Sorts on column 3 (ascending) and outputs assembly names in GCF_XXXXXXX.X format.
"""

import sys
import re
from pathlib import Path


def extract_gcf_id(filename):
    """Extract GCF_XXXXXXX.X format from the beginning of the filename."""
    match = re.match(r'(GCF_\d+\.\d+)', filename)
    if match:
        return match.group(1)
    return filename


def process_mash_file(input_file, output_file, top_n=10):
    """
    Sort mash distances by column 3 (ascending) and extract top N hits.
    
    Args:
        input_file: Path to mash distances file
        output_file: Path to output file with assembly names
        top_n: Number of top hits to extract (default: 10)
    """
    # Read and parse the mash file
    hits = []
    with open(input_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            
            fields = line.split('\t')
            if len(fields) >= 3:
                try:
                    distance = float(fields[2])
                    hits.append((distance, fields[0]))
                except ValueError:
                    print(f"Warning: Could not parse distance from line: {line}", file=sys.stderr)
                    continue
    
    # Sort by distance (column 3, ascending)
    hits.sort(key=lambda x: x[0])
    
    # Extract top N hits
    top_hits = hits[:top_n]
    
    # Extract GCF IDs and write to output file
    with open(output_file, 'w') as f:
        for distance, filename in top_hits:
            gcf_id = extract_gcf_id(filename)
            f.write(f"{gcf_id}\n")
    
    return top_hits


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python process_mash_hits.py <mash_input_file> [output_file] [top_n]")
        print("Example: python process_mash_hits.py results/3-Analysis/mash/SRR5006289_mash_distances.txt results/3-Analysis/mash/SRR5006289_best_hits.txt 10")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2] if len(sys.argv) > 2 else input_file.replace('.txt', '_best_hits.txt')
    top_n = int(sys.argv[3]) if len(sys.argv) > 3 else 10
    
    input_path = Path(input_file)
    if not input_path.exists():
        print(f"Error: Input file not found: {input_file}", file=sys.stderr)
        sys.exit(1)
    
    process_mash_file(input_file, output_file, top_n)
