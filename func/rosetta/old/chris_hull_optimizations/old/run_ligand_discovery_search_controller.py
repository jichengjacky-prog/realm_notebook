#!/usr/bin/env python3
"""
HPC‑optimized controller for launching Rosetta ligand discovery jobs on LSF.

Features:
- Clean, deterministic job list generation
- LSF‑safe array submission
- Automatic tmp-space reservation via rusage[tmp=X]
- Scratch-size sanity checking
- Drop‑in replacement for the user's original controller
"""

import os
from pathlib import Path
import subprocess
import argparse
import sys

starting_location = os.getcwd()


def run_cmd(cmd):
	result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
	if result.stdout:
		print(result.stdout, end="")
	if result.stderr:
		print(result.stderr, file=sys.stderr, end="")
	return result.returncode


def parse_args():
    p = argparse.ArgumentParser(
        description="HPC-optimized controller for launching Rosetta ligand discovery jobs on LSF."
    )

    p.add_argument("target_pdb", type=Path, help="Path to target PDB file")
    p.add_argument(
        "anchor_list",
        help="Comma-separated anchor residues, e.g. '63,87,96,179'",
    )
    p.add_argument("motifs_file", type=Path, help="Path to motifs file")
    p.add_argument("discovery_root", type=Path, help="Root directory to search for test_params")
    p.add_argument("atr", help="atr argument passed to worker")
    p.add_argument("rep", help="rep argument passed to worker")
    p.add_argument("ddg", help="ddg argument passed to worker")
    p.add_argument(
        "extra_args_file",
        type=Path,
        nargs="*",
        default=None,
        help="Optional extra args file",
    )

    return p.parse_args()


def find_test_params_dirs(discovery_root: Path):
    test_params_directories = []
    for r, d, f in os.walk(discovery_root):
        for dire in d:
            if dire == "test_params":
                test_params_directories.append(r + "/" + dire + "/")
    return test_params_directories


def build_joblists(target_pdb: Path, anchors, motifs_file: Path, test_params_directories, atr, rep, ddg, extra_args_file):
    anchors = [a.strip() for a in anchors.split(",") if a.strip()]

    all_joblist_list = []
    joblist_job_counter = 0
    joblist_file_counter = 0
    working_joblist_file = "joblist_" + str(joblist_file_counter) + ".txt"
    joblist_path = Path(starting_location + "/" + working_joblist_file)

    f = joblist_path.open("w")

    for tp_dire in test_params_directories:
        for anchor in anchors:
            line = (
                f"{target_pdb} "
                f"{anchor} "
                f"{motifs_file} "
                f"{tp_dire} "
                f"{atr} {rep} {ddg}"
            )
            if extra_args_file:
                line += f" {extra_args_file}"
            f.write(line + "\n")

            joblist_job_counter += 1

            if joblist_job_counter == 8000:
                all_joblist_list.append([joblist_path, joblist_job_counter])
                joblist_file_counter += 1
                working_joblist_file = "joblist_" + str(joblist_file_counter) + ".txt"
                joblist_path = Path(starting_location + "/" + working_joblist_file)
                joblist_job_counter = 0
                f.close()
                f = joblist_path.open("w")

    if joblist_job_counter > 0:
        f.close()
        all_joblist_list.append([joblist_path, joblist_job_counter])

    num_jobs = len(anchors) * len(test_params_directories)
    print(f"Prepared job list with {num_jobs} jobs.")

    return all_joblist_list


def ensure_logs_dir():
    logs = Path("logs")
    logs.mkdir(exist_ok=True)
    return logs


def check_scratch_size(motifs_file: Path):
    scratch_gb = int(os.environ.get("SCRATCH_SIZE_GB", 20))
    tmp_mb = scratch_gb * 1024  # LSF expects MB

    motif_size_gb = motifs_file.stat().st_size / (1024**3)
    min_required_gb = max(5, motif_size_gb * 2)

    if scratch_gb < min_required_gb:
        print(
            f"\nERROR: Requested SCRATCH_SIZE_GB={scratch_gb} GB is too small.\n"
            f"Motif file alone is {motif_size_gb:.2f} GB.\n"
            f"Minimum recommended scratch is {min_required_gb:.1f} GB.\n"
            f"Please increase SCRATCH_SIZE_GB and resubmit.\n"
        )
        sys.exit(1)

    if motif_size_gb > scratch_gb * 0.25:
        print(
            f"WARNING: Motif file ({motif_size_gb:.2f} GB) is more than 25% "
            f"of requested scratch ({scratch_gb} GB)."
        )

    return tmp_mb


def submit_joblists(all_joblist_list, tmp_mb):
    wrapper = Path("/pi/summer.thyme-umw/enamine-REAL-2.6billion/umass_chan_REAL-M_platform/rosetta/chris_hull_optimizations/lsf_scratch_wrapper.sh").resolve()

    if not wrapper.exists():
        print("ERROR: lsf_scratch_wrapper.sh not found.")
        sys.exit(1)

    for joblist in all_joblist_list:
        bsub_cmd = [
            "bsub ",
            "-J rosetta_ld[1-%s] " % joblist[1],
            "-R \"rusage[mem=10000,tmp=%s]\" " % tmp_mb,
            "-q long ",
            "-W 96:00 ",
            "-o logs/%J_%I.out ",
            "-e logs/%J_%I.err ",
            "bash %s " % wrapper + str(joblist[0])
        ]

        print("\nSubmitting LSF job array:")
        print(" ".join(bsub_cmd), "\n")

        run_cmd(" ".join(bsub_cmd))


def main():
    args = parse_args()

    target_pdb = args.target_pdb
    anchor_list = args.anchor_list
    motifs_file = args.motifs_file
    discovery_root = args.discovery_root
    atr, rep, ddg = args.atr, args.rep, args.ddg
    extra_args_file = args.extra_args_file

    test_params_directories = find_test_params_dirs(discovery_root)

    all_joblist_list = build_joblists(
        target_pdb, anchor_list, motifs_file, test_params_directories, atr, rep, ddg, extra_args_file
    )

    ensure_logs_dir()

    tmp_mb = check_scratch_size(motifs_file)

    submit_joblists(all_joblist_list, tmp_mb)


if __name__ == "__main__":
    main()
