#!/usr/bin/env python3
"""
Takes a chain file (gzipped or not) where Pan troglodytes is the TARGET
(e.g. Pan_troglodytes_vs_Homo_sapiens.chain.gz) and outputs three BED files
in Pan troglodytes coordinates:
  - <prefix>.unique.bed      : positions covered by exactly 1 chain
  - <prefix>.nonunique.bed   : positions covered by >1 chains
  - <prefix>.noalign.bed     : positions not covered by any chain

Usage: python chain_to_beds.py Pan_troglodytes_vs_Homo_sapiens.chain.gz Pan_troglodytes
"""

import sys
import gzip
import argparse
from collections import defaultdict

def open_chain(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")

def parse_chain_target_intervals(path):
    """
    Parse chain file and return TARGET-side covered intervals.
    Chain format:
      chain score tName tSize tStrand tStart tEnd qName qSize qStrand qStart qEnd id
      followed by alignment blocks: size [dt dq]

    We read tName/tSize/tStart (fields 2-6) and advance position by dt (fields[1]).
    Target coordinates are always on the + strand.
    """
    intervals = defaultdict(list)  # chrom -> [(start, end), ...]

    with open_chain(path) as f:
        t_name = None
        t_size = None
        t_start = None
        pos = None  # current position in target

        for line in f:
            line = line.strip()
            if not line:
                continue

            if line.startswith("chain"):
                fields = line.split()
                # chain score tName tSize tStrand tStart tEnd qName qSize qStrand qStart qEnd id
                t_name  = fields[2]       # e.g. NC_072398.2 (chimp)
                t_size  = int(fields[3])  # chimp chrom size
                t_start = int(fields[5])  # tStart
                pos = t_start             # walk along target from tStart

            else:
                fields = line.split()
                if not fields:
                    continue

                # Skip comment/header lines
                if fields[0].startswith("#"):
                    continue

                size = int(fields[0])

                # This alignment block covers [pos, pos+size) in target
                intervals[t_name].append((pos, pos + size))

                if len(fields) == 3:
                    # dt = gap in target, dq = gap in query
                    dt = int(fields[1])   # advance by target gap
                    pos += size + dt
                # if len(fields) == 1, it's the last block in the chain

    return intervals

def merge_adjacent(intervals):
    """
    Merge adjacent or overlapping (start, end) intervals.
    Input must be sorted.
    """
    if not intervals:
        return []
    merged = [list(intervals[0])]
    for start, end in intervals[1:]:
        if start <= merged[-1][1]:  # adjacent (start == prev end) or overlapping
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([start, end])
    return [tuple(x) for x in merged]

def compute_coverage_depth(intervals):
    """
    Given a list of (start, end) intervals (possibly overlapping),
    return a list of (start, end, depth) non-overlapping segments.
    Uses a sweep line approach.
    """
    if not intervals:
        return []

    events = []
    for s, e in intervals:
        events.append((s, +1))
        events.append((e, -1))
    events.sort()

    segments = []
    depth = 0
    prev_pos = None

    for pos, delta in events:
        if prev_pos is not None and pos > prev_pos and depth > 0:
            segments.append((prev_pos, pos, depth))
        depth += delta
        prev_pos = pos

    return segments  # [(start, end, depth), ...]

def get_chrom_sizes_from_chain(path):
    """Extract TARGET chrom sizes from chain header lines."""
    sizes = {}
    with open_chain(path) as f:
        for line in f:
            if line.startswith("chain"):
                fields = line.split()
                t_name = fields[2]
                t_size = int(fields[3])
                sizes[t_name] = t_size
    return sizes

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("chain", help="Input chain file (.chain or .chain.gz), species-of-interest as TARGET")
    parser.add_argument("prefix", help="Output prefix (e.g. Pan_troglodytes)")
    args = parser.parse_args()

    print("Parsing chain file...", file=sys.stderr)
    intervals_by_chrom = parse_chain_target_intervals(args.chain)
    chrom_sizes = get_chrom_sizes_from_chain(args.chain)

    unique_path    = f"{args.prefix}.unique.bed"
    nonunique_path = f"{args.prefix}.nonunique.bed"
    noalign_path   = f"{args.prefix}.noalign.bed"

    with open(unique_path, "w") as fu, \
         open(nonunique_path, "w") as fn, \
         open(noalign_path, "w") as fna:

        for chrom, size in sorted(chrom_sizes.items()):
            raw_intervals = intervals_by_chrom.get(chrom, [])

            # Compute depth per segment
            depth_segments = compute_coverage_depth(raw_intervals)

            # Collect intervals per class before merging
            unique_ivs    = []
            nonunique_ivs = []
            noalign_ivs   = []

            covered_end = 0
            for start, end, depth in depth_segments:
                if start > covered_end:
                    noalign_ivs.append((covered_end, start))
                if depth == 1:
                    unique_ivs.append((start, end))
                else:
                    nonunique_ivs.append((start, end))
                covered_end = end

            if covered_end < size:
                noalign_ivs.append((covered_end, size))

            # Merge adjacent/overlapping intervals within each class
            for start, end in merge_adjacent(sorted(unique_ivs)):
                fu.write(f"{chrom}\t{start}\t{end}\n")
            for start, end in merge_adjacent(sorted(nonunique_ivs)):
                fn.write(f"{chrom}\t{start}\t{end}\n")
            for start, end in merge_adjacent(sorted(noalign_ivs)):
                fna.write(f"{chrom}\t{start}\t{end}\n")

    print(f"Written: {unique_path}", file=sys.stderr)
    print(f"Written: {nonunique_path}", file=sys.stderr)
    print(f"Written: {noalign_path}", file=sys.stderr)

if __name__ == "__main__":
    main()
