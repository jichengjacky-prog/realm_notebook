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
import sys
from pathlib import Path
import subprocess

starting_location = os.getcwd()

def run(cmd: str):
    """Run a shell command safely."""
    subprocess.run(cmd, shell=True, check=True)


# ---------------------------------------------------------------------
# Parse controller arguments
# ---------------------------------------------------------------------
if len(sys.argv) < 8:
    print(
        "Usage: run_ligand_discovery_search_controller.py "
        "<target_pdb> <anchor_list> <motifs_file> <discovery_root> "
        "<atr> <rep> <ddg> [extra_args_file]"
    )
    sys.exit(1)

target_pdb = Path(sys.argv[1])
anchor_list = sys.argv[2]        # e.g. "63,87,96,179"
motifs_file = Path(sys.argv[3])
discovery_root = Path(sys.argv[4])
atr, rep, ddg = sys.argv[5:8]
extra_args_file = Path(sys.argv[8]) if len(sys.argv) >= 9 else None


#look down the discovery root for all test_params directories

test_params_directories = []

for r,d,f in os.walk(discovery_root):
    for dire in d:
        if dire == "test_params":
            test_params_directories.append(r + "/" + dire + "/")

# ---------------------------------------------------------------------
# Build job list (one job per anchor residue)
# ---------------------------------------------------------------------
anchors = [a.strip() for a in anchor_list.split(",") if a.strip()]

joblist_path = Path(starting_location + "/joblist.txt")
with joblist_path.open("w") as f:
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

num_jobs = len(anchors) * len(test_params_directories)
print(f"Prepared job list with {num_jobs} jobs.")


# ---------------------------------------------------------------------
# Ensure logs directory exists
# ---------------------------------------------------------------------
logs = Path("logs")
logs.mkdir(exist_ok=True)


# ---------------------------------------------------------------------
# SCRATCH SIZE COMMENTARY + SANITY CHECK
#
# The worker script performs all heavy I/O inside a loopback-mounted
# filesystem located in /tmp. The wrapper script creates this filesystem
# with a size determined by SCRATCH_SIZE_GB.
#
# To prevent oversubscription of /tmp on compute nodes, we request the
# same amount of space from LSF using:
#
#     -R "rusage[tmp=X]"
#
# where X is in megabytes.
#
# This ensures:
#   • LSF will only dispatch the job to nodes with enough free /tmp space.
#   • Multiple jobs cannot overrun a node's local disk.
#   • The wrapper and scheduler always agree on the required scratch size.
#
# Users may override the default (10 GB) by setting:
#
#     export SCRATCH_SIZE_GB=20
#
# before running this controller.
#
# SANITY CHECK:
#   We ensure the scratch size is large enough to hold:
#     • The motif file
#     • Temporary Rosetta outputs
#     • Placement PDBs
#   Minimum recommended scratch = max(5 GB, 2 × motif file size)
# ---------------------------------------------------------------------

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


# ---------------------------------------------------------------------
# Submit LSF job array
# ---------------------------------------------------------------------
wrapper = Path("/pi/summer.thyme-umw/enamine-REAL-2.6billion/umass_chan_REAL-M_platform/rosetta/chris_hull_optimizations/lsf_scratch_wrapper.sh").resolve()

if not wrapper.exists():
    print("ERROR: lsf_scratch_wrapper.sh not found.")
    sys.exit(1)

bsub_cmd = (
    f"bsub "
    f"-J rosetta_ld[1-{num_jobs}] "
    f"-R \"rusage[mem=10000,tmp={tmp_mb}]\" "
    f"-q long "
    f"-W 96:00 "
#    f"-o logs/%J_%I.out "
#    f"-e logs/%J_%I.err "
#    f"bash {wrapper} $(sed -n \"\\$LSB_JOBINDEX\"p " + starting_location + "/joblist.txt)"
    f"bash {wrapper} " + str(joblist_path)
)

print("\nSubmitting LSF job array:")
print(bsub_cmd, "\n")

run(bsub_cmd)

