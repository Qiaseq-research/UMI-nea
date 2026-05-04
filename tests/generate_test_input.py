#!/usr/bin/env python3
"""
Generate synthetic UMI clustering input files for reproducibility testing.

Output format (UMI-nea input): <group_id> TAB <UMI_seq> TAB <read_count>
Sorted by read_count descending within each group.

Usage:
    python3 generate_test_input.py <n_founders> <umi_len> <seed> <output_file>
    python3 generate_test_input.py 1000 18 42 sim2.input
"""

import sys
import random
from collections import defaultdict

BASES = "ACGT"


def random_umi(length, rng):
    return "".join(rng.choice(BASES) for _ in range(length))


def mutate(umi, n_mutations, rng):
    """Apply up to n_mutations substitutions or indels to umi."""
    seq = list(umi)
    for _ in range(n_mutations):
        if not seq:
            break
        pos = rng.randint(0, len(seq) - 1)
        op = rng.randint(0, 2)
        if op == 0:  # substitution
            seq[pos] = rng.choice([b for b in BASES if b != seq[pos]])
        elif op == 1 and len(seq) > 1:  # deletion
            seq.pop(pos)
        else:  # insertion
            seq.insert(pos, rng.choice(BASES))
    return "".join(seq[:len(umi) + 2])  # cap length growth


def generate(n_founders, umi_len, seed, out_path,
             mean_founder_reads=15, children_per_founder=8,
             max_child_reads=4, max_dist=2):
    rng = random.Random(seed)

    # Generate distinct founder sequences
    founder_seqs = set()
    while len(founder_seqs) < n_founders:
        founder_seqs.add(random_umi(umi_len, rng))
    founders = list(founder_seqs)

    counts = defaultdict(int)

    for f in founders:
        # Assign read count from a geometric-like distribution
        rc = max(1, int(rng.expovariate(1.0 / mean_founder_reads)))
        counts[f] += rc

        # Generate children (sequencing errors off the founder)
        n_children = rng.randint(0, children_per_founder * 2)
        for _ in range(n_children):
            dist = rng.choices([1, 2], weights=[3, 1])[0]
            child = mutate(f, dist, rng)
            if child != f:
                child_rc = rng.randint(1, max_child_reads)
                counts[child] += child_rc

    # Sort by count descending (UMI-nea expects this order)
    sorted_umis = sorted(counts.items(), key=lambda x: -x[1])

    with open(out_path, "w") as fh:
        for umi, rc in sorted_umis:
            fh.write(f"1\t{umi}\t{rc}\n")

    total = len(sorted_umis)
    print(f"Generated {total} unique UMIs from {n_founders} founders → {out_path}")
    return total


if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: generate_test_input.py <n_founders> <umi_len> <seed> <output>")
        sys.exit(1)
    n_founders = int(sys.argv[1])
    umi_len    = int(sys.argv[2])
    seed       = int(sys.argv[3])
    out_path   = sys.argv[4]
    generate(n_founders, umi_len, seed, out_path)
