#!/usr/bin/env python3
"""
Shared utility functions for the Rosetta discovery pipeline.

Consolidates reusable logic from:
  - fix_condensed_param_file_spacing.py
  - extract_single_param_from_condensed_file.py
  - score_placed_ligands_with_filtering.py
  - unified_discovery_controller.py
"""

import os
import re
import shutil
import subprocess
import sys
import tempfile


# ---------------------------------------------------------------------------
# Shell helpers
# ---------------------------------------------------------------------------

def run_cmd(cmd, cwd=None, description=""):
	"""Run a shell command. Returns returncode."""
	if description:
		print(f"  {description}")
	result = subprocess.run(cmd, shell=True, capture_output=True, text=True, cwd=cwd)
	if result.stdout:
		print(result.stdout, end="")
	if result.stderr:
		print(result.stderr, file=sys.stderr, end="")
	return result.returncode


# ---------------------------------------------------------------------------
# Enamine library path helpers
# ---------------------------------------------------------------------------

def chunk_to_path(chunk_str):
	"""Convert a 5-digit zero-padded chunk string to (superchunk_str, chunk_str).
	Example: '00123' -> ('1', '00123')"""
	superchunk_str = str(int(chunk_str[:3]))
	return superchunk_str, chunk_str


# ---------------------------------------------------------------------------
# Condensed params file helpers
# ---------------------------------------------------------------------------

def fix_params_spacing_text(params_text):
	"""Fix column spacing in Rosetta .params content (in-memory).
	Takes a string of params file content, returns fixed string."""
	lines = []
	for line in params_text.splitlines(True):
		if line.startswith("ATOM"):
			parts = line.strip().split()
			if len(parts) >= 5:
				lines.append("ATOM %-4s %-4s %-4s %.2f\n" % (
					parts[1], parts[2], parts[3], float(parts[4])))
		elif line.startswith("BOND_TYPE"):
			parts = line.strip().split()
			if len(parts) >= 4:
				lines.append("BOND_TYPE %-4s %-4s %-4s\n" % (
					parts[1], parts[2], parts[3]))
		elif line.startswith("CHI"):
			parts = line.strip().split()
			if len(parts) >= 6:
				lines.append("CHI %i %-4s %-4s %-4s %-4s\n" % (
					int(parts[1]), parts[2], parts[3], parts[4], parts[5]))
		else:
			lines.append(line)
	return "".join(lines)


def fix_params_spacing(params_file, output_dir=None):
	"""Fix column spacing in a Rosetta .params file (file-based wrapper).
	Reads *params_file*, writes a fixed version, and replaces the original (or
	writes to *output_dir*).

	Returns the path to the fixed file, or None on failure.
	"""
	if not os.path.isfile(params_file):
		print(f"ERROR: params file not found: {params_file}")
		return None

	try:
		with open(params_file, "r") as fh:
			content = fh.read()
		fixed = fix_params_spacing_text(content)

		if output_dir:
			os.makedirs(output_dir, exist_ok=True)
			dest = os.path.join(output_dir, os.path.basename(params_file))
			with open(dest, "w") as fh:
				fh.write(fixed)
			return dest
		else:
			with open(params_file, "w") as fh:
				fh.write(fixed)
			return params_file
	except Exception as e:
		print(f"ERROR fixing params spacing in {params_file}: {e}")
		return None


# ---------------------------------------------------------------------------
# Condensed params extraction
# ---------------------------------------------------------------------------

def extract_single_param_text(condensed_text, conf_identifier):
	"""Extract a single conformer's params from condensed text content (in-memory).
	Takes the full text of a condensed params file, returns the .params content
	as a string for the chosen conformer, or None on failure."""
	condense_dict = {}
	in_keys = False
	in_params = False
	output_lines = []

	for line in condensed_text.splitlines(True):
		line_no_newline = line.rstrip("\n")

		# Section headers
		if len(line_no_newline.split()) == 1 and line_no_newline == "KEYS":
			in_keys = True
			continue
		if len(line_no_newline.split()) == 1 and line_no_newline == "PARAMS":
			in_params = True
			in_keys = False
			continue

		# Read key definitions
		if in_keys:
			if not line_no_newline.startswith("_"):
				print(f"Bad key line (missing underscore): {line_no_newline}")
				continue
			if len(line_no_newline.split()) != 2:
				print(f"Bad key line (expected 2 parts): {line_no_newline}")
				continue

			key, entry = line_no_newline.split()
			if key in condense_dict:
				print(f"Duplicate key {key} in dictionary!")
				continue
			condense_dict[key] = entry

		# Read data lines
		if in_params:
			if line_no_newline.count(":") > 1:
				print(f"Line has too many colons: {line_no_newline}")
				continue

			header_side, data_side = line_no_newline.split(":")

			# Check for conserved marker (*)
			conserved = header_side.endswith("*")
			if conserved:
				header_side = header_side[:-1]

			# Translate header
			if header_side in condense_dict:
				header_side = condense_dict[header_side]

			# Pick the correct conformer's data
			split_data = data_side.split(",")
			if conserved:
				my_data = split_data[0]
			else:
				my_data = split_data[conf_identifier - 1]

			# Translate each token
			my_data_split = my_data.split()
			for i in range(len(my_data_split)):
				if my_data_split[i] in condense_dict:
					my_data_split[i] = condense_dict[my_data_split[i]]

			my_data = " " + " ".join(my_data_split)
			output_lines.append(header_side + my_data + "\n")

	return "".join(output_lines)


def extract_single_param(condensed_file_path, conf_identifier, ligand_name, output_dir=None):
	"""Extract a single conformer's params from a condensed params file and write
	it as a standalone .params file readable by Rosetta (file-based wrapper).

	Args:
		condensed_file_path: path to the condensed params text file
		conf_identifier:    1-based conformer index within the condensed file (1-15)
		ligand_name:        base name for the output .params file
		output_dir:         directory to write the .params file (default: cwd)

	Returns:
		path to the written .params file, or None on failure
	"""
	try:
		with open(condensed_file_path, "r") as fh:
			content = fh.read()
	except IOError as e:
		print(f"ERROR: Cannot open file: {e}")
		return None

	params_text = extract_single_param_text(content, conf_identifier)
	if params_text is None:
		return None

	if output_dir:
		os.makedirs(output_dir, exist_ok=True)
		output_path = os.path.join(output_dir, f"{ligand_name}.params")
	else:
		output_path = f"{ligand_name}.params"

	with open(output_path, "w") as fh:
		fh.write(params_text)
	return output_path


# ---------------------------------------------------------------------------
# Placement PDB scoring
# ---------------------------------------------------------------------------

def parse_placement_scores(pdb_path, residue_index_dict=None):
	"""Parse a Rosetta ligand-discovery placement PDB file and extract scoring
	metrics from its header comments.

	Args:
		pdb_path:            path to the placed ligand PDB file
		residue_index_dict:  optional dict mapping original->translated residue indices

	Returns:
		dict with keys: ddg, total_motifs, significant_motifs, real_motif_ratio,
		hbond_motif_count, hbond_motif_energy_sum, found_motif_residues (list),
		or None if the file cannot be read.
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
				elif line.startswith("Placement motifs: Motifs made against significant residues count:"):
					scores["significant_motifs"] = float(line.split()[-1].strip())

				# Real motif ratio
				elif line.startswith("Placement motifs: Real motif ratio:"):
					scores["real_motif_ratio"] = float(line.split()[-1].strip())

				# Per-motif hbond lines
				elif ": Placement motif " in line:
					# Extract residue index from the motif description
					index_match = re.search(r"Hbond_score.*?_(\d{3})[A-Z]?", line)
					if index_match:
						index = index_match.group(1)
						# Translate if a key dict is provided
						if residue_index_dict and index in residue_index_dict:
							index = residue_index_dict[index]
						scores["found_motif_residues"].append(index)

					# Extract hbond score (last colon-separated field)
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


def compute_weighted_total(scores_dict, score_weights=None):
	"""Apply score weights to a scores dict (from parse_placement_scores) and
	return (weighted_total, weighted_breakdown_dict).

	Default weights are 1.0 for every term.
	"""
	defaults = {
		"ddg": 1.0,
		"total_motifs": 1.0,
		"significant_motifs": 1.0,
		"real_motif_ratio": 1.0,
		"hbond_motif_count": 1.0,
		"hbond_motif_energy_sum": 1.0,
		"closest_autodock_recovery_rmsd": 1.0,
		"closest_autodock_recovery_ddg": 1.0,
		"strain_energy": 1.0,
	}
	if score_weights:
		defaults.update(score_weights)

	weighted = {}
	total = 0.0
	for term, weight in defaults.items():
		raw = scores_dict.get(term, 0.0)
		w = raw * weight
		weighted[term] = w
		total += w

	return total, weighted


def load_score_weights(weights_csv_path):
	"""Load score weights from a CSV file (term,weight per line, no header).
	Returns a dict of {term: weight}."""
	weights = {}
	if not weights_csv_path or not os.path.isfile(weights_csv_path):
		return weights

	with open(weights_csv_path, "r") as fh:
		for line in fh:
			line = line.strip()
			if not line:
				continue
			parts = line.split(",")
			if len(parts) >= 2:
				try:
					weights[parts[0].strip()] = float(parts[1].strip())
				except ValueError:
					print(f"WARNING: bad weight line: {line}")
	return weights


def load_residue_index_key(key_file_path):
	"""Load a residue-index translation key CSV.
	Format: res_type,original_index,translated_index[,difference]
	Returns dict mapping original_index -> translated_index."""
	key_dict = {}
	if not key_file_path or not os.path.isfile(key_file_path):
		return key_dict

	with open(key_file_path, "r") as fh:
		for line in fh:
			if line.startswith("res_type"):
				continue
			parts = line.strip().split(",")
			if len(parts) >= 3:
				key_dict[parts[1].strip()] = parts[2].strip()
	return key_dict


# ---------------------------------------------------------------------------
# Rosetta ligand discovery search
# ---------------------------------------------------------------------------

def run_rosetta_discovery_search(target_pdb, anchor_residue_string, motifs_file,
								  test_params_dir, discovery_root, atr, rep, ddg,
								  extra_args_file=None, work_dir=None):
	"""Run a single Rosetta ligand discovery search.

	Writes a Rosetta args file, executes the search via singularity, and
	post-processes the output (placements directory, scoring, compression).

	Args:
		target_pdb:             path to target PDB file
		anchor_residue_string:  comma-separated Rosetta-indexed anchor residues
		motifs_file:            path to motifs file
		test_params_dir:        path to test_params directory (must be named "test_params")
		discovery_root:         root directory of the discovery pipeline
		atr, rep, ddg:          fa_atr, fa_rep, ddg score cutoffs
		extra_args_file:        optional file with additional Rosetta args
		work_dir:               working directory (default: current directory)

	Returns:
		True on success, False on failure.
	"""
	if work_dir:
		os.makedirs(work_dir, exist_ok=True)
		orig_cwd = os.getcwd()
		os.chdir(work_dir)
	else:
		orig_cwd = None

	try:
		# Ensure test_params_dir ends with /
		if not test_params_dir.endswith("/"):
			test_params_dir = test_params_dir + "/"

		# Write Rosetta args file
		with open("args", "w") as args_fh:
			args_fh.write("#keep seed constant\n")
			args_fh.write("-constant_seed 1\n")
			args_fh.write("#ignore unrecognized residues to help mitigate crashes\n")
			args_fh.write("-ignore_unrecognized_res\n")
			args_fh.write("#handle ligand repeats if using multiple anchor residues\n")
			args_fh.write("-in::file::override_database_params true\n")
			args_fh.write("#constrain coordinates\n")
			args_fh.write("-constrain_relax_to_start_coords\n")
			args_fh.write("#keep all placements\n")
			args_fh.write("-best_pdbs_to_keep 0\n")

			# User-input dependent
			args_fh.write("#mapped protein system\n")
			args_fh.write("-s /input/" + os.path.basename(target_pdb) + "\n")
			args_fh.write("#mapped motifs file\n")
			args_fh.write("-motif_filename /input/" + os.path.basename(motifs_file) + "\n")
			args_fh.write("#mapped test_params directory\n")
			args_fh.write("-params_directory_path /input/" + os.path.basename(test_params_dir.rstrip("/")) + "/\n")
			args_fh.write("#rosetta-indexed anchor residue index/indices\n")
			args_fh.write("-protein_discovery_locus " + anchor_residue_string + "\n")
			args_fh.write("#fa_atr cutoff\n")
			args_fh.write("-fa_atr_cutoff = " + str(atr) + "\n")
			args_fh.write("#fa_rep cutoff\n")
			args_fh.write("-fa_rep_cutoff = " + str(rep) + "\n")
			args_fh.write("#ddg cutoff\n")
			args_fh.write("-ddg_cutoff = " + str(ddg) + "\n")

			# Extra user args
			if extra_args_file and os.path.isfile(extra_args_file):
				args_fh.write("###################################################\n")
				args_fh.write("#extra user args from: " + extra_args_file + "\n")
				with open(extra_args_file, "r") as extra_fh:
					args_fh.write(extra_fh.read())

		# Build singularity command
		rosetta_sif = os.path.join(discovery_root, "sif", "rosetta_condensed_6_25_2024.sif")
		rosetta_cmd = " ".join([
			"singularity exec",
			"--bind " + test_params_dir + ":/input/test_params/",
			"--bind " + os.getcwd() + "/args:/input/args",
			"--bind " + target_pdb + ":/input/" + os.path.basename(target_pdb),
			"--bind " + motifs_file + ":/input/" + os.path.basename(motifs_file),
			rosetta_sif,
			"/rosetta/source/bin/ligand_discovery_search_protocol.linuxgccrelease @/input/args",
		])
		##debug print the command
		print(rosetta_cmd)

		print("Running Rosetta discovery search...")
		if run_cmd(rosetta_cmd) != 0:
			print("ERROR: Rosetta discovery search failed")
			return False

		# Post-processing: organize placements
		run_cmd("mkdir -p placements")
		run_cmd("mv *pdb placements 2>/dev/null; true")

		os.chdir("placements")
		# Rename each PDB with the anchor residue prefix
		for r, d, f in os.walk(os.getcwd()):
			for file in f:
				if file.endswith(".pdb") and r == os.getcwd():
					os.rename(file, "res" + anchor_residue_string + "_" + file)

		# Run placement scoring
		score_script = os.path.join(discovery_root, "func", "rosetta",
									"score_placed_ligands_with_filtering.py")
		run_cmd("python " + score_script)

		# Copy score CSVs up one level
		run_cmd("cp *csv .. 2>/dev/null; true")

		os.chdir("..")

		# Compress placements
		run_cmd("tar -czf placements.tar.gz placements")
		run_cmd("rm -drf placements")

		# Dehydrate to minimize overhead
		dehydrate_script = os.path.join(discovery_root, "func", "tidying",
										"shrink_placement_pdbs_to_placement_and_surrounding_residues.py")
		run_cmd("python " + dehydrate_script + " " + target_pdb)

		print("Rosetta discovery search complete.")
		return True

	except Exception as e:
		print(f"ERROR in Rosetta discovery search: {e}")
		return False

	finally:
		if orig_cwd:
			os.chdir(orig_cwd)


# ---------------------------------------------------------------------------
# CLI wrappers (kept for backward compatibility)
# ---------------------------------------------------------------------------

if __name__ == "__main__":
	import argparse

	parser = argparse.ArgumentParser(description="Discovery pipeline utilities")
	sub = parser.add_subparsers(dest="command")

	# fix-params
	p_fix = sub.add_parser("fix-params", help="Fix spacing in a condensed params file")
	p_fix.add_argument("params_file", help="Path to the .params file")

	# extract-param
	p_ext = sub.add_parser("extract-param", help="Extract single conformer from condensed file")
	p_ext.add_argument("condensed_file", help="Path to condensed params file")
	p_ext.add_argument("conf_id", type=int, help="Conformer identifier (1-based)")
	p_ext.add_argument("ligand_name", help="Base ligand name")
	p_ext.add_argument("--output-dir", "-o", default=None, help="Output directory")

	# parse-placement
	p_place = sub.add_parser("parse-placement", help="Parse scores from a placement PDB")
	p_place.add_argument("pdb_path", help="Path to the placed ligand PDB")
	p_place.add_argument("--residue-key", "-k", default=None, help="Residue index translation key CSV")

	args = parser.parse_args()

	if args.command == "fix-params":
		result = fix_params_spacing(args.params_file)
		print(f"Fixed: {result}" if result else "Failed.")

	elif args.command == "extract-param":
		result = extract_single_param(args.condensed_file, args.conf_id,
									  args.ligand_name, args.output_dir)
		print(f"Written: {result}" if result else "Failed.")

	elif args.command == "parse-placement":
		key_dict = load_residue_index_key(args.residue_key) if args.residue_key else None
		scores = parse_placement_scores(args.pdb_path, key_dict)
		if scores:
			total, weighted = compute_weighted_total(scores)
			print(f"ddg={scores['ddg']} motifs={scores['total_motifs']} "
				  f"sig_motifs={scores['significant_motifs']} ratio={scores['real_motif_ratio']} "
				  f"hbond_count={scores['hbond_motif_count']} hbond_sum={scores['hbond_motif_energy_sum']}")
			print(f"Found motif residues: {scores['found_motif_residues']}")
		else:
			print("Failed to parse placement PDB.")
			sys.exit(1)

	else:
		parser.print_help()
