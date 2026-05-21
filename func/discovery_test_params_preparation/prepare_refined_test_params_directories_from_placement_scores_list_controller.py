#!/usr/bin/env python3
"""
Multithreaded controller to prepare refined test params directories from a placements summary file.
This script extracts ligand SMILES from placement files and submits per-ligand jobs (via bsub) to
run the helper script that generates conformers and params.

Usage:
    python prepare_refined_test_params_directories_from_placement_scores_list_controller.py SUMMARY_FILE WORKING_LOCATION [--realm-location PATH] [--license-key KEY] [--workers N]
"""

import os
import sys
import argparse
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from threading import Lock


def parse_args():
    parser = argparse.ArgumentParser(
        description='Prepare conformers for multiple placements from a summary file'
    )
    parser.add_argument(
        'summary_file',
        type=argparse.FileType('r'),
        help='Path to placement summary file (CSV format)'
    )
    parser.add_argument(
        'working_location',
        help='Working directory to use for processing'
    )
    parser.add_argument(
        '--license-key',
        dest='license_key',
        default='',
        help='Optional Conformator license key (if needed for container)'
    )
    parser.add_argument(
        '--realm-location',
        dest='realm_location',
        default=None,
        help='Optional repository root (defaults to current working directory)'
    )
    parser.add_argument(
        '--workers',
        dest='workers',
        type=int,
        default=4,
        help='Number of worker threads to use (default: 4)'
    )
    return parser.parse_args()


def run_cmd(cmd, cwd=None):
    """Run a shell command, capture and print stdout/stderr, return exit code."""
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True, cwd=cwd)
    if result.stdout:
        print(result.stdout, end="")
    if result.stderr:
        print(result.stderr, file=sys.stderr, end="")
    return result.returncode


args = parse_args()
summary_file = args.summary_file
working_location = os.path.abspath(args.working_location)
license_key = args.license_key
workers = args.workers
realm_location = args.realm_location or os.getcwd()

# Ensure working location exists
os.makedirs(working_location, exist_ok=True)

# thread-safety primitives and state
dir_lock = Lock()
ligand_name_set = set()
# numeric directory counter (top-level directories named 0,1,2,... each contain up to 100 ligand dirs)
_dir_counter = 0


def _get_dir_counter():
    global _dir_counter
    return _dir_counter


def _increment_dir_counter():
    global _dir_counter
    _dir_counter += 1
    os.makedirs(os.path.join(working_location, str(_dir_counter)), exist_ok=True)


def allocate_directory_for_ligand(name):
    """Allocate (and create) a per-ligand directory under a numeric top-level folder.
    Ensures each top-level folder holds up to 100 ligand subdirectories.
    Returns the absolute path to the ligand directory.
    """
    global _dir_counter
    with dir_lock:
        top_dir = os.path.join(working_location, str(_dir_counter))
        os.makedirs(top_dir, exist_ok=True)
        target = os.path.join(top_dir, name)
        if os.path.exists(target):
            return target
        os.makedirs(target, exist_ok=True)
        # if the top_dir has reached 100 subdirs, advance counter
        try:
            n = len([n for n in os.listdir(top_dir) if os.path.isdir(os.path.join(top_dir, n))])
        except Exception:
            n = 0
        if n >= 100:
            _increment_dir_counter()
        return target


def process_placement_line(line):
    """Process a single summary-file line: extract placement, get SMILES, and submit bsub job."""
    if not line or line.strip() == "":
        return None
    if line.startswith("file,ddg,total_motifs,significant_motifs"):
        return None

    full_file = line.split(",")[0].strip()
    if not full_file:
        return None

    compressed_placements_path = full_file.split("/placements/")[0]
    placement_file = full_file.split("/")[-1]

    # ligand name is in the placement file name as the 3rd-from-last underscore-separated field
    parts = placement_file.split("_")
    if len(parts) < 3:
        return None
    ligand_name = parts[-3]

    # deduplicate so we only submit one job per ligand
    with dir_lock:
        if ligand_name in ligand_name_set:
            return None
        ligand_name_set.add(ligand_name)

    # create ligand directory
    target_dir = allocate_directory_for_ligand(ligand_name)

    # extract the placement into the ligand directory
    tar_cmd = f"tar -xzf {compressed_placements_path}/placements.tar.gz placements/{placement_file} --strip-components=1 -C ."
    run_cmd(tar_cmd, cwd=target_dir)

    # extract ligand pdb
    run_cmd(f"grep HETATM {placement_file} > ligand.pdb", cwd=target_dir)

    # convert to SMILES
    run_cmd("obabel -ipdb ligand.pdb -osmi -O ligand.smi -xn", cwd=target_dir)

    # read smiles
    smiles_string = ""
    try:
        with open(os.path.join(target_dir, "ligand.smi"), 'r') as sf:
            for l in sf:
                smiles_string = l.strip()
                break
    except Exception:
        smiles_string = ""

    # cleanup local files
    run_cmd("rm -f *pdb *smi", cwd=target_dir)

    # prepare and submit the job using bsub
    helper = os.path.join(realm_location, 'func', 'discovery_test_params_preparation', 'prepare_refined_test_params_directories_from_placement_scores_list.py')
    cmd_list = ['bsub -q short -n 1', '-W', '1:00', '-u', '-R', 'rusage[mem=4000]',
        '-o', 'log.out', '-e', 'log.err',
        'python', helper,
        str(smiles_string),
        str(ligand_name),
        str(license_key),
        realm_location
    ]
    cmd = ' '.join(cmd_list)
    print(cmd)
    run_cmd(cmd, cwd=target_dir)

    # small throttle
    run_cmd('sleep 0.1', cwd=target_dir)
    return ligand_name


def main():
    # dispatch work across worker threads
    futures = []
    with ThreadPoolExecutor(max_workers=workers) as exe:
        for line in summary_file:
            futures.append(exe.submit(process_placement_line, line))

        for fut in as_completed(futures):
            try:
                _ = fut.result()
            except Exception as e:
                print(f"Error processing placement: {e}", file=sys.stderr)


if __name__ == '__main__':
    main()
