#!/usr/bin/env python3
"""Run Rosetta ligand discovery search.
CLI wrapper — import run_rosetta_discovery_search from utils.py for programmatic use."""

import os, sys

# Allow running directly or importing from utils
try:
	from utils import run_rosetta_discovery_search
except ImportError:
	sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "discovery"))
	from utils import run_rosetta_discovery_search


if __name__ == "__main__":
	import argparse
	parser = argparse.ArgumentParser(
		description="Run Rosetta ligand discovery search for a single test_params directory."
	)
	parser.add_argument("target_pdb", help="Target PDB file path")
	parser.add_argument("anchor_residue_string", help="Comma-separated Rosetta-indexed anchor residue(s)")
	parser.add_argument("motifs_file", help="Motifs file path")
	parser.add_argument("test_params_dir", help="Path to the test_params directory")
	parser.add_argument("discovery_directory_root", help="Root directory of the discovery pipeline")
	parser.add_argument("atr", help="fa_atr cutoff")
	parser.add_argument("rep", help="fa_rep cutoff")
	parser.add_argument("ddg", help="ddg cutoff")
	parser.add_argument("--extra-args-file", default="", help="Optional file with additional Rosetta args")
	args = parser.parse_args()

	success = run_rosetta_discovery_search(
		args.target_pdb, args.anchor_residue_string, args.motifs_file,
		args.test_params_dir, args.discovery_directory_root,
		args.atr, args.rep, args.ddg, args.extra_args_file
	)
	sys.exit(0 if success else 1)