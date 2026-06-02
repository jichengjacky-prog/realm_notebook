#!/usr/bin/env python3
"""Fix column spacing in Rosetta .params files produced from condensed format.
CLI wrapper — import fix_params_spacing from utils.py for programmatic use."""

import os, sys

# Allow running directly or importing from utils
try:
	from utils import fix_params_spacing
except ImportError:
	# Fallback: if run standalone, add script dir to path
	sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
	from utils import fix_params_spacing


if __name__ == "__main__":
	if len(sys.argv) < 2:
		print("Usage: python fix_condensed_param_file_spacing.py <params_file>")
		sys.exit(1)

	result = fix_params_spacing(sys.argv[1])
	if result:
		print(f"Fixed params file: {result}")
	else:
		print("Failed to fix params file.")
		sys.exit(1)

