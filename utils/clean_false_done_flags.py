#!/usr/bin/env python3
"""
Find and remove false per-ligand .done flags where weighted_scores.csv count
is below the expected anchor residue count.

Also removes the batch-level rosetta_round1.done if any ligand in that batch
had a false .done flag, so Snakemake will re-process the batch on the next run.

Usage:
  python utils/clean_false_done_flags.py --dry-run    # count only
  python utils/clean_false_done_flags.py --delete      # actually remove
"""

import os
import sys
import argparse
import glob

BASE_DIR = "/pi/summer.thyme-umw/Ji_rosetta_discovery"

# Both config_state3.yaml and config_state4.yaml have 10 anchor residues.
# If this changes, update accordingly.
EXPECTED_CSV_COUNT = 10


def clean_false_done_flags(dry_run=True):
    """Remove <ligand>.done files where CSVs < EXPECTED_CSV_COUNT.
    If any ligand in a batch was corrected, also remove rosetta_round1.done."""

    total_lig_done_removed = 0
    total_batch_done_removed = 0
    batches_affected = 0

    for state in ["snakemake_state3", "snakemake_state4"]:
        state_path = os.path.join(BASE_DIR, "output", state, "tmp")
        if not os.path.isdir(state_path):
            continue

        for batch_name in sorted(os.listdir(state_path)):
            batch_dir = os.path.join(state_path, batch_name)
            if not batch_name.startswith("batch_") or not os.path.isdir(batch_dir):
                continue

            batch_done_flag = os.path.join(batch_dir, "rosetta_round1.done")
            if not os.path.exists(batch_done_flag):
                continue

            round1_dir = os.path.join(batch_dir, "round1")
            if not os.path.isdir(round1_dir):
                continue

            batch_lig_removed = 0

            for entry in sorted(os.listdir(round1_dir)):
                lig_path = os.path.join(round1_dir, entry)
                if not os.path.isdir(lig_path):
                    continue

                lig_done_file = os.path.join(round1_dir, f"{entry}.done")
                if not os.path.exists(lig_done_file):
                    continue

                # Count weighted_scores.csv files under this ligand's anchor dirs
                csv_files = glob.glob(os.path.join(lig_path, "*", "weighted_scores.csv"))
                csv_count = len(csv_files)

                if csv_count < EXPECTED_CSV_COUNT:
                    batch_lig_removed += 1
                    if dry_run:
                        print(f"  [WOULD REMOVE] {lig_done_file}  ({csv_count}/{EXPECTED_CSV_COUNT} CSVs)")
                    else:
                        try:
                            os.remove(lig_done_file)
                            print(f"  [REMOVED] {lig_done_file}  ({csv_count}/{EXPECTED_CSV_COUNT} CSVs)")
                        except OSError as e:
                            print(f"  [ERROR] {lig_done_file}: {e}")

            if batch_lig_removed > 0:
                total_lig_done_removed += batch_lig_removed
                batches_affected += 1
                if dry_run:
                    print(f"  [WOULD REMOVE] {batch_done_flag}  ({batch_lig_removed} false ligand .done flags in batch)")
                else:
                    try:
                        os.remove(batch_done_flag)
                        total_batch_done_removed += 1
                        print(f"  [REMOVED] {batch_done_flag}  ({batch_lig_removed} false ligand .done flags in batch)")
                    except OSError as e:
                        print(f"  [ERROR] {batch_done_flag}: {e}")

                print()

    return total_lig_done_removed, total_batch_done_removed, batches_affected


def main():
    parser = argparse.ArgumentParser(
        description="Clean up false .done flags (CSVs < expected anchor count)"
    )
    parser.add_argument("--dry-run", action="store_true", default=True,
                        help="Only report what would be removed (default)")
    parser.add_argument("--delete", action="store_true",
                        help="Actually delete the false done flags")
    args = parser.parse_args()

    dry_run = not args.delete
    if args.delete:
        print("=== DELETING false done flags ===\n")
    else:
        print("=== DRY RUN: Finding false done flags ===\n")

    n_lig, n_batch, n_affected = clean_false_done_flags(dry_run=dry_run)

    print(f"{'='*60}")
    if dry_run:
        print(f"DRY RUN SUMMARY:")
        print(f"  False ligand .done flags:  {n_lig:,}")
        print(f"  Batches to re-run:         {n_batch:,} (of {n_affected:,} affected)")
        print(f"\nRun with --delete to actually remove these flags.")
    else:
        print(f"REMOVAL SUMMARY:")
        print(f"  Ligand .done flags removed:  {n_lig:,}")
        print(f"  Batch rosetta_round1.done:   {n_batch:,}")
        print(f"  Batches affected:            {n_affected:,}")


if __name__ == "__main__":
    main()
