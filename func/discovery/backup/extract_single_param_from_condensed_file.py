#!/usr/bin/env python3
"""Extract a single conformer's params from a condensed params file.
CLI wrapper — import extract_single_param from utils.py for programmatic use."""

import os, sys

# Allow running directly or importing from utils
try:
	from utils import extract_single_param
except ImportError:
	sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
	from utils import extract_single_param


if __name__ == "__main__":
	if len(sys.argv) < 4:
		print("Usage: python extract_single_param_from_condensed_file.py <condensed_file> <conf_identifier> <ligand_name> [output_dir]")
		sys.exit(1)

	condensed_file_path = sys.argv[1]
	conf_identifier = int(sys.argv[2])
	ligand_name = sys.argv[3]
	output_dir = sys.argv[4] if len(sys.argv) > 4 else None

	result = extract_single_param(condensed_file_path, conf_identifier, ligand_name, output_dir)
	if result:
		print(f"Params file written to: {result}")
	else:
		print("Failed to extract params file.")
		sys.exit(1)

