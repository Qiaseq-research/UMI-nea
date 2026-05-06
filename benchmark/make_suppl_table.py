#!/usr/bin/env python3
"""
Merge all per-config performance.txt files into a single performance.csv
in the UMI-nea publication supplemental table format.

Usage: python make_suppl_table.py [benchmark_dir] [output_csv]
Defaults: benchmark_dir=. output_csv=performance.csv
"""
import sys
import os
import glob
import pandas as pd

benchmark_dir = sys.argv[1] if len(sys.argv) > 1 else "."
output_csv    = sys.argv[2] if len(sys.argv) > 2 else os.path.join(benchmark_dir, "performance.csv")

# Collect all per-config performance.txt files
pattern = os.path.join(benchmark_dir, "sim_*", "performance.txt")
files = sorted(glob.glob(pattern))
if not files:
    sys.exit(f"No performance.txt files found under {benchmark_dir}")

frames = []
for f in files:
    try:
        df = pd.read_csv(f, sep=r"\s+")
        frames.append(df)
        print(f"  loaded {f}  ({len(df)} rows)")
    except Exception as e:
        print(f"  SKIP {f}: {e}")

if not frames:
    sys.exit("No data loaded.")

df = pd.concat(frames, ignore_index=True)

# ── tool name normalisation to match publication labels ─────────────────────
tool_map = {
    "UMI-nea":  "UMI-nea",
    "umi-tools": "UMI-tools",
    "calib":    "Calib",
    "UMIC-seq": "UMIc-seq",
}
df["tool"] = df["tool"].map(tool_map).fillna(df["tool"])

# ── read_type column ────────────────────────────────────────────────────────
read_type_map = {12: "short_reads", 18: "short_reads",
                 25: "long_reads",  50: "long_reads"}
df["read_type"] = df["umi_len"].map(read_type_map).fillna("unknown")

# ── column order matching plot scripts ──────────────────────────────────────
front_cols = [
    "num_founder", "mean_children_num", "variance_children_num",
    "umi_len", "read_type", "err_rate",
    "insertion-deletion-substitution", "replicate",
    "total_umi", "simulated_mean_children_num", "simulated_variance_children_num",
    "substitution_base", "indel_base",
    "simulated_insertion-deletion-substitution",
    "substitution_only_umi", "indel_umi", "uniq_umi",
    "tool", "clustering_threshold", "thread",
    "runtime_in_sec", "dedup_umi_cluster",
    "V-measure", "homogeneity_score", "completeness_score",
    "RPU_cutoff", "RPU_cutoff_model", "estimated_molecule",
]
# keep any extra columns at the end
extra = [c for c in df.columns if c not in front_cols]
df = df[[c for c in front_cols if c in df.columns] + extra]

# ── sort: founder, umi_len, err_rate, tool, replicate ───────────────────────
sort_cols = ["num_founder", "umi_len", "err_rate", "tool", "replicate"]
sort_cols = [c for c in sort_cols if c in df.columns]
df = df.sort_values(sort_cols).reset_index(drop=True)

df.to_csv(output_csv, index=False)
print(f"\nWrote {len(df)} rows → {output_csv}")

# ── quick summary ────────────────────────────────────────────────────────────
print("\nMean V-measure by tool × read_type × num_founder:")
summary = (df.groupby(["tool", "read_type", "num_founder"])["V-measure"]
             .agg(mean="mean", std="std", n="count")
             .round(4))
print(summary.to_string())
