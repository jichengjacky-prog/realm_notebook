#!/usr/bin/env python3
"""
Regenerate weighted_scores.csv from raw_scores.csv files.

Reads each raw_scores.csv, applies SCORE_WEIGHTS to compute the weighted
total, and writes weighted_scores.csv alongside it.  Individual metric
columns are stored as raw (unweighted) values — only the 'total' column
reflects the weighted combination.

Usage:
    python utils/regenerate_weighted_scores.py output/snakemake_state4/tmp/
    python utils/regenerate_weighted_scores.py output/ --dry-run
"""

import argparse
import csv
import os
import sys
from pathlib import Path

# ── Score weights (must match score_placed_ligands_with_filtering.py) ─────
SCORE_WEIGHTS = {
    "ddg": -1.0,
    "total_motifs": 0.0,
    "significant_motifs": 0.0,
    "real_motif_ratio": 0.0,
    "hbond_motif_count": 0.0,
    "hbond_motif_energy_sum": 0.0,
}

# Columns that appear in both CSV files (order must match the CSV header)
METRIC_COLUMNS = [
    "ddg",
    "total_motifs",
    "significant_motifs",
    "real_motif_ratio",
    "hbond_motif_count",
    "hbond_motif_energy_sum",
]

CSV_HEADER = (
    "file,ddg,total_motifs,significant_motifs,real_motif_ratio,"
    "hbond_motif_count,hbond_motif_energy_sum,total\n"
)


def compute_weighted_total(row: dict) -> float:
    """Apply SCORE_WEIGHTS to a row dict and return the weighted total."""
    total = 0.0
    for col in METRIC_COLUMNS:
        try:
            raw = float(row.get(col, 0))
        except (ValueError, TypeError):
            raw = 0.0
        total += raw * SCORE_WEIGHTS.get(col, 0.0)
    return total


def regenerate_one(raw_csv_path: Path, dry_run: bool = False) -> bool:
    """Read raw_scores.csv and write weighted_scores.csv in the same directory."""
    weighted_csv_path = raw_csv_path.parent / "weighted_scores.csv"

    rows = []
    try:
        with open(raw_csv_path, "r", newline="") as fh:
            reader = csv.DictReader(fh)
            for row in reader:
                rows.append(row)
    except Exception as e:
        print(f"  [ERROR] Failed to read {raw_csv_path}: {e}", file=sys.stderr)
        return False

    if not rows:
        print(f"  [SKIP] {raw_csv_path} is empty")
        return False

    if dry_run:
        print(f"  [DRY-RUN] Would write {len(rows)} row(s) to {weighted_csv_path}")
        return True

    try:
        with open(weighted_csv_path, "w", newline="") as fh:
            fh.write(CSV_HEADER)
            for row in rows:
                file_path = row.get("file", "")
                weighted_total = compute_weighted_total(row)
                # Write raw metric values + weighted total
                metric_str = ",".join(
                    f"{float(row.get(col, 0)):.6f}" for col in METRIC_COLUMNS
                )
                fh.write(f"{file_path},{metric_str},{weighted_total:.6f}\n")
    except Exception as e:
        print(f"  [ERROR] Failed to write {weighted_csv_path}: {e}", file=sys.stderr)
        return False

    return True


def main():
    parser = argparse.ArgumentParser(
        description="Regenerate weighted_scores.csv from raw_scores.csv files."
    )
    parser.add_argument(
        "root",
        type=Path,
        help="Root directory (tmp/) containing batch_* subdirectories.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Show what would be done without writing any files.",
    )
    parser.add_argument(
        "--batch",
        type=str,
        default=None,
        help="Process only a specific batch (e.g. 'batch_3570').",
    )
    args = parser.parse_args()

    root = args.root.resolve()
    if not root.exists():
        print(f"Error: directory not found: {root}", file=sys.stderr)
        sys.exit(1)

    # Discover batch subdirectories to avoid resolving the entire tree at once
    if args.batch:
        batch_dirs = [root / args.batch]
        if not batch_dirs[0].is_dir():
            print(f"Error: batch directory not found: {batch_dirs[0]}", file=sys.stderr)
            sys.exit(1)
    else:
        batch_dirs = sorted(
            p for p in root.iterdir() if p.is_dir() and p.name.startswith("batch_")
        )

    if not batch_dirs:
        print(f"No batch_* directories found under {root}")
        sys.exit(0)

    print(f"Processing {len(batch_dirs)} batch director{'y' if len(batch_dirs) == 1 else 'ies'} under {root}")
    if args.dry_run:
        print("DRY RUN — no files will be written.\n")
    else:
        print()

    total_regenerated = 0
    total_found = 0
    for batch_dir in batch_dirs:
        raw_files = sorted(batch_dir.rglob("raw_scores.csv"))
        if not raw_files:
            continue
        total_found += len(raw_files)
        batch_regenerated = 0
        for fpath in raw_files:
            if regenerate_one(fpath, dry_run=args.dry_run):
                batch_regenerated += 1
        total_regenerated += batch_regenerated
        if batch_regenerated > 0:
            action = "Would regenerate" if args.dry_run else "Regenerated"
            print(f"  {batch_dir.name}: {action} {batch_regenerated} file(s)")

    print(f"\n{'Would regenerate' if args.dry_run else 'Regenerated'} {total_regenerated} of {total_found} file(s) across {len(batch_dirs)} batch(es).")


if __name__ == "__main__":
    main()
