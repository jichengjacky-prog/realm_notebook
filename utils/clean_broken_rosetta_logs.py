#!/usr/bin/env python3
"""
Find and remove broken/stale Rosetta log files.

Cleans up three types of broken logs:
  1. Orphaned ros1_*.out/.err/.py  — where <LIGAND>.done exists (job completed, logs stale)
  2. Empty rosetta_round1.log      — Snakemake rule logs with 0 bytes (failed batch jobs)
  3. ROSETTA_CRASH.log             — individual Rosetta anchor crash logs

Usage:
  python clean_broken_rosetta_logs.py --dry-run    # count only
  python clean_broken_rosetta_logs.py --delete      # actually remove
"""

import os
import sys
import argparse
import glob

BASE_DIR = "/pi/summer.thyme-umw/Ji_rosetta_discovery"


def clean_orphaned_ros1_files(dry_run=True):
    """Remove orphaned ros1_*.out, ros1_*.err, ros1_*.py where .done exists."""
    total_removed = 0
    total_size_bytes = 0
    batches_affected = 0

    for state in ["snakemake_state3", "snakemake_state4"]:
        state_path = os.path.join(BASE_DIR, "output", state, "tmp")
        if not os.path.isdir(state_path):
            continue

        for batch_name in sorted(os.listdir(state_path)):
            batch_dir = os.path.join(state_path, batch_name)
            if not batch_name.startswith("batch_") or not os.path.isdir(batch_dir):
                continue

            round1_dir = os.path.join(batch_dir, "round1")
            if not os.path.isdir(round1_dir):
                continue

            batch_removed = 0
            batch_size = 0

            for pattern in ["ros1_*.out", "ros1_*.err", "ros1_*.py"]:
                for filepath in glob.glob(os.path.join(round1_dir, pattern)):
                    basename = os.path.basename(filepath)
                    if not basename.startswith("ros1_"):
                        continue

                    inner = basename[5:]
                    if inner.endswith(".out") or inner.endswith(".err"):
                        lig_name = inner[:-4]
                    elif inner.endswith(".py"):
                        lig_name = inner[:-3]
                    else:
                        continue

                    done_file = os.path.join(round1_dir, f"{lig_name}.done")
                    if os.path.exists(done_file):
                        file_size = os.path.getsize(filepath)
                        batch_removed += 1
                        batch_size += file_size
                        if not dry_run:
                            try:
                                os.remove(filepath)
                            except OSError as e:
                                print(f"  [ERROR] {filepath}: {e}")

            if batch_removed > 0:
                batches_affected += 1
                total_removed += batch_removed
                total_size_bytes += batch_size
                action = "would be removed" if dry_run else "removed"
                print(f"  Batch {batch_name}: {batch_removed} orphaned ros1_* files ({batch_size:,} bytes) {action}")

    return total_removed, total_size_bytes, batches_affected


def clean_empty_rosetta_logs(dry_run=True):
    """Remove empty (0-byte) rosetta_round1.log files."""
    total_removed = 0

    for state in ["snakemake_state3", "snakemake_state4"]:
        state_path = os.path.join(BASE_DIR, "output", state, "tmp")
        if not os.path.isdir(state_path):
            continue

        for batch_name in sorted(os.listdir(state_path)):
            batch_dir = os.path.join(state_path, batch_name)
            if not batch_name.startswith("batch_") or not os.path.isdir(batch_dir):
                continue

            log_file = os.path.join(batch_dir, "rosetta_round1.log")
            if os.path.isfile(log_file) and os.path.getsize(log_file) == 0:
                total_removed += 1
                if dry_run:
                    print(f"  [WOULD REMOVE] {log_file} (0 bytes)")
                else:
                    try:
                        os.remove(log_file)
                        print(f"  [REMOVED] {log_file} (0 bytes)")
                    except OSError as e:
                        print(f"  [ERROR] {log_file}: {e}")

    return total_removed


def clean_rosetta_crash_logs(dry_run=True):
    """Remove ROSETTA_CRASH.log files."""
    total_removed = 0
    total_size_bytes = 0

    for state in ["snakemake_state3", "snakemake_state4"]:
        state_path = os.path.join(BASE_DIR, "output", state, "tmp")
        if not os.path.isdir(state_path):
            continue

        for batch_name in sorted(os.listdir(state_path)):
            batch_dir = os.path.join(state_path, batch_name)
            if not batch_name.startswith("batch_") or not os.path.isdir(batch_dir):
                continue

            round1_dir = os.path.join(batch_dir, "round1")
            if not os.path.isdir(round1_dir):
                continue

            # Walk the ligand/anchor directory tree looking for ROSETTA_CRASH.log
            for root, dirs, files in os.walk(round1_dir):
                for fname in files:
                    if fname == "ROSETTA_CRASH.log":
                        filepath = os.path.join(root, fname)
                        file_size = os.path.getsize(filepath)
                        total_removed += 1
                        total_size_bytes += file_size
                        if dry_run:
                            if total_removed <= 20:
                                print(f"  [WOULD REMOVE] {filepath} ({file_size} bytes)")
                        else:
                            try:
                                os.remove(filepath)
                            except OSError as e:
                                print(f"  [ERROR] {filepath}: {e}")

    return total_removed, total_size_bytes


def main():
    parser = argparse.ArgumentParser(description="Clean up broken Rosetta log files")
    parser.add_argument("--dry-run", action="store_true", default=True,
                        help="Only count files, do not delete (default)")
    parser.add_argument("--delete", action="store_true",
                        help="Actually delete the files")
    args = parser.parse_args()

    dry_run = not args.delete
    if args.delete:
        print("=== DELETING broken Rosetta log files ===\n")
    else:
        print("=== DRY RUN: Finding broken Rosetta log files ===\n")

    # 1. Clean orphaned ros1_* files
    print("--- Orphaned ros1_*.out/.err/.py (where .done exists) ---")
    n1, s1, b1 = clean_orphaned_ros1_files(dry_run=dry_run)
    print(f"  Total: {n1:,} files, {s1/(1024*1024):.2f} MB, {b1} batches\n")

    # 2. Clean empty rosetta_round1.log
    print("--- Empty rosetta_round1.log (0 bytes) ---")
    n2 = clean_empty_rosetta_logs(dry_run=dry_run)
    print(f"  Total: {n2:,} files\n")

    # 3. Clean ROSETTA_CRASH.log
    print("--- ROSETTA_CRASH.log files ---")
    n3, s3 = clean_rosetta_crash_logs(dry_run=dry_run)
    print(f"  Total: {n3:,} files, {s3/(1024*1024):.2f} MB\n")

    print(f"{'='*60}")
    total_files = n1 + n2 + n3
    total_size = s1 + s3
    if dry_run:
        print(f"GRAND TOTAL (dry run): {total_files:,} files, {total_size/(1024*1024):.2f} MB would be removed")
        print(f"Run with --delete to actually remove these files.")
    else:
        print(f"GRAND TOTAL: {total_files:,} files removed, {total_size/(1024*1024):.2f} MB freed")


if __name__ == "__main__":
    main()
