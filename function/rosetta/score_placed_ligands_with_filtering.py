#!/usr/bin/env python3
"""
Score placed ligand PDB files from Rosetta's ligand_discovery_search.

Walks the current directory looking for placement PDB files, parses scoring
metrics from their headers, filters by real_motif_ratio, computes weighted
totals, and writes raw_scores.csv and weighted_scores.csv.

Usage (called from within the placements/ directory):
    python score_placed_ligands_with_filtering.py
"""

import os
import re
import sys

# ── Configuration ──────────────────────────────────────────────────────────
MINIMUM_REAL_MOTIF_RATIO = 0.25  # placements below this ratio are ignored

# Score: ddg only  ─  score = -ddg  (higher = better)
#   ddg is negative (more negative = better binding), so -ddg makes it positive.
#   Filtered by real_motif_ratio > {MINIMUM_REAL_MOTIF_RATIO} below.
SCORE_WEIGHTS = {
    "ddg": -1.0,              # -ddg: better binders get higher positive score
    "total_motifs": 0.0,
    "significant_motifs": 0.0,
    "real_motif_ratio": 0.0,
    "hbond_motif_count": 0.0,
    "hbond_motif_energy_sum": 0.0,
}


# ── PDB parsing ────────────────────────────────────────────────────────────

def parse_placement_scores(pdb_path):
    """Parse a Rosetta placement PDB and extract scoring metrics from headers.

    Returns dict with keys: ddg, total_motifs, significant_motifs,
    real_motif_ratio, hbond_motif_count, hbond_motif_energy_sum,
    found_motif_residues (list).  Returns None on failure.
    """
    if not os.path.isfile(pdb_path):
        print(f"ERROR: PDB file not found: {pdb_path}")
        return None

    scores = {
        "ddg": 0.0,
        "total_motifs": 0.0,
        "significant_motifs": 0.0,
        "real_motif_ratio": 0.0,
        "hbond_motif_count": 0,
        "hbond_motif_energy_sum": 0.0,
        "found_motif_residues": [],
    }

    try:
        with open(pdb_path, "r") as fh:
            for line in fh:
                # ddG
                if line.startswith("Scoring: Post-HighResDock system ddG:"):
                    scores["ddg"] = float(line.split()[-1].strip())

                # Total motifs
                elif line.startswith("Placement motifs: Total motifs made:"):
                    scores["total_motifs"] = float(line.split()[-1].strip())

                # Significant motifs
                elif line.startswith(
                    "Placement motifs: Motifs made against significant residues count:"
                ):
                    scores["significant_motifs"] = float(line.split()[-1].strip())

                # Real motif ratio
                elif line.startswith("Placement motifs: Real motif ratio:"):
                    scores["real_motif_ratio"] = float(line.split()[-1].strip())

                # Per-motif hbond lines
                elif ": Placement motif " in line:
                    index_match = re.search(r"Hbond_score.*?_(\d{3})[A-Z]?", line)
                    if index_match:
                        scores["found_motif_residues"].append(index_match.group(1))

                    hbond_score_str = line.split(":")[-1].strip()
                    try:
                        hbond_score = float(hbond_score_str)
                    except ValueError:
                        continue

                    scores["hbond_motif_energy_sum"] += hbond_score
                    if hbond_score != 0.0:
                        scores["hbond_motif_count"] += 1

    except Exception as e:
        print(f"ERROR reading {pdb_path}: {e}")
        return None

    return scores


def compute_weighted_total(scores_dict, weights=None):
    """Apply score weights and return (raw_total, weighted_total, breakdown_dict)."""
    if weights is None:
        weights = {}

    defaults = {
        "ddg": 1.0,
        "total_motifs": 1.0,
        "significant_motifs": 1.0,
        "real_motif_ratio": 1.0,
        "hbond_motif_count": 1.0,
        "hbond_motif_energy_sum": 1.0,
    }
    defaults.update(weights)

    raw_total = 0.0
    weighted_total = 0.0
    breakdown = {}

    for term, weight in defaults.items():
        raw = scores_dict.get(term, 0.0)
        raw_total += raw
        weighted_total += raw * weight
        breakdown[term] = raw

    return raw_total, weighted_total, breakdown


# ── Main ───────────────────────────────────────────────────────────────────

def main():
    location = os.getcwd()
    print(f"Scoring placements in: {location}")

    # Collect PDB files (exclude minipose files; all others are placement PDBs)
    pdb_files = []
    for root, dirs, files in os.walk(location):
        for f in files:
            if f.endswith(".pdb") and "minipose" not in f:
                pdb_files.append(os.path.join(root, f))

    if not pdb_files:
        print("No placement PDB files found.")
        # Write empty CSVs so the pipeline doesn't break
        _write_empty_csvs()
        return

    print(f"Found {len(pdb_files)} placement PDB(s)")

    # Open output files
    raw_fh = open("raw_scores.csv", "w")
    weighted_fh = open("weighted_scores.csv", "w")

    # Write headers
    header = (
        "file,ddg,total_motifs,significant_motifs,real_motif_ratio,"
        "hbond_motif_count,hbond_motif_energy_sum,total\n"
    )
    raw_fh.write(header)
    weighted_fh.write(header)

    scored_count = 0
    for pdb_path in sorted(pdb_files):
        parsed = parse_placement_scores(pdb_path)
        if parsed is None:
            continue

        # Filter by real motif ratio
        if parsed["real_motif_ratio"] < MINIMUM_REAL_MOTIF_RATIO:
            print(
                f"  SKIP {os.path.basename(pdb_path)}: "
                f"real_motif_ratio={parsed['real_motif_ratio']:.4f} < {MINIMUM_REAL_MOTIF_RATIO}"
            )
            continue

        raw_total, weighted_total, breakdown = compute_weighted_total(
            parsed, SCORE_WEIGHTS
        )

        # Write raw scores
        raw_fh.write(
            f"{pdb_path},{breakdown['ddg']:.6f},{breakdown['total_motifs']:.6f},"
            f"{breakdown['significant_motifs']:.6f},{breakdown['real_motif_ratio']:.6f},"
            f"{breakdown['hbond_motif_count']},{breakdown['hbond_motif_energy_sum']:.6f},"
            f"{raw_total:.6f}\n"
        )

        # Write weighted scores — individual columns store RAW (unweighted) values
        # so downstream steps can read the actual metrics (e.g. real_motif_ratio).
        # Only the 'total' column reflects the weighted combination.
        weighted_fh.write(
            f"{pdb_path},"
            f"{breakdown['ddg']:.6f},"
            f"{breakdown['total_motifs']:.6f},"
            f"{breakdown['significant_motifs']:.6f},"
            f"{breakdown['real_motif_ratio']:.6f},"
            f"{breakdown['hbond_motif_count']:.6f},"
            f"{breakdown['hbond_motif_energy_sum']:.6f},"
            f"{weighted_total:.6f}\n"
        )

        scored_count += 1
        print(
            f"  {os.path.basename(pdb_path)}: "
            f"raw={raw_total:.4f}  weighted={weighted_total:.4f}  "
            f"ddg={parsed['ddg']:.4f}  motifs={parsed['real_motif_ratio']:.4f}"
        )

    raw_fh.close()
    weighted_fh.close()
    print(f"Scored {scored_count} placements (filtered from {len(pdb_files)} total)")


def _write_empty_csvs():
    """Write empty score CSVs with just headers so downstream steps don't break."""
    header = (
        "file,ddg,total_motifs,significant_motifs,real_motif_ratio,"
        "hbond_motif_count,hbond_motif_energy_sum,total\n"
    )
    for name in ("raw_scores.csv", "weighted_scores.csv"):
        with open(name, "w") as fh:
            fh.write(header)


if __name__ == "__main__":
    main()
