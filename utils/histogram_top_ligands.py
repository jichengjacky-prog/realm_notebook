#!/usr/bin/env python3
"""
histogram_top_ligands.py — histogram of the top-N round-1 ligands by ddg

Reads a top-ligands CSV produced by utils/manual_filter_top_round1.py
(header: ligand,file,ddg,total_motifs,significant_motifs,real_motif_ratio,
hbond_motif_count,hbond_motif_energy_sum,total — rows already ranked by
ddg, most negative first) and plots two histograms:

  1. ddg distribution of the top-N ligands
  2. real_motif_ratio distribution of the top-N ligands

Saves a PNG and prints percentile summary statistics.  No filesystem scan
is performed — the top list is read directly.

Usage:
    python3 utils/histogram_top_ligands.py \
        [--input output/snakemake_state4/top_ligands_round1.txt] \
        [--top-n 100000] [--out path/to/hist.png]
"""

import argparse
import csv
import os
import sys

import matplotlib
matplotlib.use("Agg")  # headless
import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REALM = os.path.dirname(HERE)  # utils/ is directly under the realm root


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input",
        default=os.path.join(REALM, "output", "snakemake_state4", "top_ligands_round1.txt"),
        help="top-ligands CSV from manual_filter_top_round1.py",
    )
    parser.add_argument("--top-n", type=int, default=100000)
    parser.add_argument(
        "--out",
        default=None,
        help="output PNG path (default: <input dir>/top<top-n>_histograms.png)",
    )
    args = parser.parse_args()

    in_path = os.path.abspath(args.input)
    if not os.path.isfile(in_path):
        print(f"ERROR: input file not found: {in_path}", file=sys.stderr)
        sys.exit(1)
    out_png = args.out or os.path.join(
        os.path.dirname(in_path), f"top{args.top_n//1000}k_histograms.png"
    )

    print(f"input file : {in_path}")
    print(f"top n      : {args.top_n}")
    print(f"output png : {out_png}")

    # ── Read the top list ────────────────────────────────────────────────
    ddgs = []
    mrs = []
    with open(in_path, "r", newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            if len(ddgs) >= args.top_n:
                break
            try:
                ddgs.append(float(row["ddg"]))
                mrs.append(float(row["real_motif_ratio"]))
            except (KeyError, ValueError, TypeError):
                continue
    ddgs = np.array(ddgs, dtype=float)
    mrs = np.array(mrs, dtype=float)
    n = len(ddgs)
    if n == 0:
        print("ERROR: no rows read", file=sys.stderr)
        sys.exit(1)
    print(f"read {n} ligands (requested top {args.top_n})")

    def pct(a, q):
        return float(np.percentile(a, q))

    print("\n── ddg percentiles (more negative = better) ──")
    for q in (1, 5, 25, 50, 75, 95, 99):
        print(f"  P{q:02d}: {pct(ddgs, q):8.3f}")
    print(f"  min: {ddgs.min():.3f}   max: {ddgs.max():.3f}")

    print("\n── real_motif_ratio percentiles ──")
    for q in (1, 5, 25, 50, 75, 95, 99):
        print(f"  P{q:02d}: {pct(mrs, q):8.3f}")
    print(f"  min: {mrs.min():.3f}   max: {mrs.max():.3f}")

    # ── Histograms ──────────────────────────────────────────────────────
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

    ax1.hist(ddgs, bins=100, color="#4C72B0", edgecolor="white", linewidth=0.3)
    ax1.axvline(pct(ddgs, 50), color="red", linestyle="--", linewidth=1,
                label=f"median {pct(ddgs, 50):.2f}")
    ax1.axvline(pct(ddgs, 25), color="orange", linestyle=":", linewidth=1,
                label=f"P25 {pct(ddgs, 25):.2f}")
    ax1.axvline(pct(ddgs, 75), color="orange", linestyle=":", linewidth=1,
                label=f"P75 {pct(ddgs, 75):.2f}")
    ax1.set_xlabel("ddg (kcal/mol)")
    ax1.set_ylabel("count")
    ax1.set_title(f"Top {n:,} ligands — ddg distribution")
    ax1.legend(fontsize=8)

    ax2.hist(mrs, bins=50, color="#55A868", edgecolor="white", linewidth=0.3)
    ax2.axvline(pct(mrs, 50), color="red", linestyle="--", linewidth=1,
                label=f"median {pct(mrs, 50):.3f}")
    ax2.set_xlabel("real_motif_ratio")
    ax2.set_ylabel("count")
    ax2.set_title(f"Top {n:,} ligands — real_motif_ratio distribution")
    ax2.legend(fontsize=8)

    fig.suptitle(f"Round-1 top {n:,} ligands by ddg", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(out_png, dpi=150)
    print(f"\nhistograms saved to {out_png}")


if __name__ == "__main__":
    main()
