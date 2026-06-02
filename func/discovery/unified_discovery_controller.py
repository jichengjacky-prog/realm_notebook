#!/usr/bin/env python3
"""
Unified Discovery Controller

Combines the entire discovery pipeline into a single controller:
1. Takes a shapedb results list (conformer info + chunk/subchunk)
2. Extracts conformer params from the Enamine library
3. Runs Rosetta ligand discovery search
4. Scores placements with filtering
5. Maintains a top-1000 deduplicated ligand heap

Usage:
    python unified_discovery_controller.py \
        --shapedb-list <file> \
        --target-pdb <file> \
        --anchor-residues <res1,res2,...> \
        --motifs-file <file> \
        --output-dir <dir> \
        [--realm-location <path>] \
        [--workers <n>] \
        [--max-ligands <n>] \
        [--top-hits <n>] \
        [--conformator-license <file>] \
        [--atr <val>] [--rep <val>] [--ddg <val>]
"""

import argparse
import concurrent.futures
import heapq
import os
import shutil
import subprocess
import sys
import tempfile
from collections import defaultdict
from pathlib import Path
from threading import Lock
from utils import *


# Thread-safe heap and state management
heap_lock = Lock()
ligand_seen = set()
ligand_scores_heap = []  # min-heap of (-score, ligand_name_conf, data_dict)
ligand_seen_lock = Lock()

ERROR_LOG = None


def _write_error(msg):
	"""Write an error message to both stderr and the error log file."""
	print(msg, file=sys.stderr)
	if ERROR_LOG:
		ERROR_LOG.write(msg + "\n")
		ERROR_LOG.flush()


def extract_conformer_params(ligand_name, conf_num, chunk, subchunk, 
							  enamine_path, tmp_root=None, output_dir=None):
	"""Extract and fix a single conformer's params from the Enamine library.

	Args:
		output_dir: if given, writes the fixed .params file there and returns its path.
		            If None, returns the params content as a string (no file written).

	Returns:
		str (params text) if output_dir is None, or str (file path) if output_dir is given.
		Returns None on failure.
	"""
	superchunk_str, chunk_str = chunk_to_path(chunk)
	
	# Pipe tar output directly to memory — no temp files on disk
	tar_cmd = (
		f"tar -xzf {enamine_path}/{superchunk_str}/{chunk_str}/"
		f"condensed_params_and_db_{subchunk}.tar.gz "
		f"condensed_params_and_db_{subchunk}/single_conf_params/"
		f"{ligand_name}_shorthand_params.txt "
		f"--strip-components=2 -O"
	)
	result = subprocess.run(tar_cmd, shell=True, capture_output=True, text=True)
	if result.returncode != 0 or not result.stdout:
		print(f"ERROR: Failed to extract params for {ligand_name}_{conf_num}")
		if result.stderr:
			print(result.stderr, file=sys.stderr, end="")
		return None
	
	# All processing is in-memory from here — no intermediate files
	params_text = extract_single_param_text(result.stdout, conf_num)
	if not params_text:
		print(f"ERROR: Failed to extract conformer {conf_num} for {ligand_name}")
		return None

	params_text = fix_params_spacing_text(params_text)
	
	# Write to output_dir if requested, otherwise return text directly
	if output_dir:
		os.makedirs(output_dir, exist_ok=True)
		dest = os.path.join(output_dir, f"{ligand_name}_{conf_num}.params")
		with open(dest, "w") as fh:
			fh.write(params_text)
		return dest
	
	return params_text


def extract_smiles_from_sdf(sdf_file):
	"""Extract SMILES string from first line of SDF file (assumes obabel format)."""
	try:
		# Use obabel to convert SDF to SMILES
		result = subprocess.run(
			f"obabel {sdf_file} -O - -osmi 2>/dev/null | head -1",
			shell=True, capture_output=True, text=True
		)
		if result.stdout:
			smiles = result.stdout.strip().split()[0]  # Get first field
			return smiles
	except Exception as e:
		print(f"ERROR extracting SMILES: {e}")
	return None


def generate_conformers_and_params(ligand_name, smiles, output_dir, 
									realm_location, conformator_container_sif,
									num_conformers=150, license_key="", tmp_root=None, work_dir=None):
	"""Generate 150 conformers from SMILES and create Rosetta params files.
	If work_dir is not given, a temp directory is created under tmp_root."""
	
	# Use caller-provided work_dir, or create one
	_own_work_dir = None
	if work_dir is None:
		work_dir = tempfile.mkdtemp(dir=tmp_root)
		_own_work_dir = work_dir
	
	try:
		# Write SMILES to file
		smiles_file = os.path.join(work_dir, f"{ligand_name}.smi")
		with open(smiles_file, "w") as fh:
			fh.write(smiles)
		
		# Run conformator in singularity container
		if license_key:
			conf_cmd = (
				f"singularity exec {conformator_container_sif} bash -lc "
				f"'/conformator_for_container/conformator_1.2.1/conformator --license \"{license_key}\" && "
				f"/conformator_for_container/conformator_1.2.1/conformator -i {ligand_name}.smi "
				f"-o {ligand_name}_confs.sdf --keep3d --hydrogens -n {num_conformers} -v 0'"
			)
		else:
			conf_cmd = (
				f"singularity exec {conformator_container_sif} bash -lc "
				f"'/conformator_for_container/conformator_1.2.1/conformator -i {ligand_name}.smi "
				f"-o {ligand_name}_confs.sdf --keep3d --hydrogens -n {num_conformers} -v 0'"
			)
		
		if run_cmd(conf_cmd, cwd=work_dir) != 0:
			print(f"ERROR: Conformator failed for {ligand_name}")
			return []
		
		# Split conformers into individual SDFs
		split_cmd = f"obabel -isdf {ligand_name}_confs.sdf -O {ligand_name}_.sdf -m"
		if run_cmd(split_cmd, cwd=work_dir) != 0:
			print(f"ERROR: obabel split failed for {ligand_name}")
			return []
		
		# Generate params for each conformer
		params_files = []
		for i in range(1, num_conformers + 1):
			sdf_file = os.path.join(work_dir, f"{ligand_name}_{i}.sdf")
			if not os.path.exists(sdf_file):
				continue
			
			# Run molfile_to_params
			params_cmd = (
				f"singularity exec {conformator_container_sif} python "
				f"/conformator_for_container/molfile_to_params.py {sdf_file} "
				f"-n {ligand_name}_{i} --keep-names --long-names --clobber --no-pdb"
			)
			
			if run_cmd(params_cmd, cwd=work_dir) != 0:
				print(f"WARNING: molfile_to_params failed for {ligand_name}_{i}")
				continue
			
			params_file = os.path.join(work_dir, f"{ligand_name}_{i}.params")
			if os.path.exists(params_file):
				# Copy to output directory
				dest = os.path.join(output_dir, f"{ligand_name}_{i}.params")
				shutil.copy2(params_file, dest)
				params_files.append(f"{ligand_name}_{i}.params")
		
		print(f"Generated {len(params_files)} conformer params for {ligand_name}")
		return params_files
	
	except Exception as e:
		print(f"ERROR generating conformers for {ligand_name}: {e}")
		return []
	
	finally:
		if _own_work_dir:
			shutil.rmtree(work_dir, ignore_errors=True)


def create_test_params_dir(ligands_list, batch_dir, realm_location, enamine_path, tmp_root=None):
	"""Create a test_params directory from a list of (ligand_name, conf_num, chunk, subchunk) tuples."""
	test_params_dir = os.path.join(batch_dir, "test_params")
	os.makedirs(test_params_dir, exist_ok=True)
	
	# Create required Rosetta files
	run_cmd("touch exclude_pdb_component_list.txt patches.txt", cwd=test_params_dir)
	
	# Create residue_types.txt header
	res_types_file = os.path.join(test_params_dir, "residue_types.txt")
	with open(res_types_file, "w") as fh:
		fh.write("## atom_type_set and mm-atom_type_set for Rosetta\n")
		fh.write("TYPE_SET_MODE full_atom\n")
		fh.write("ATOM_TYPE_SET fa_standard\n")
		fh.write("ELEMENT_SET default\n")
		fh.write("MM_ATOM_TYPE_SET fa_standard\n")
		fh.write("ORBITAL_TYPE_SET fa_standard\n")
		fh.write("## Params files\n")
	
	# Extract and fix each conformer
	failed_ligands = []
	for ligand_name, conf_num, chunk, subchunk in ligands_list:
		try:
			params_file = extract_conformer_params(
				ligand_name, conf_num, chunk, subchunk, enamine_path,
				batch_dir, output_dir=test_params_dir
			)
			
			if params_file:
				# Add to residue_types.txt
				with open(res_types_file, "a") as fh:
					fh.write(f"{ligand_name}_{conf_num}.params\n")
			else:
				failed_ligands.append(f"{ligand_name}_{conf_num}")
		
		except Exception as e:
			print(f"ERROR processing {ligand_name}_{conf_num}: {e}")
			failed_ligands.append(f"{ligand_name}_{conf_num}")
	
	if failed_ligands:
		print(f"WARNING: Failed to process {len(failed_ligands)} ligands")
	
	return test_params_dir


def generate_and_add_conformers_to_test_params(test_params_dir, ligands_list, 
											   realm_location, enamine_path,
											   num_conformers=150, license_key="", tmp_root=None):
	"""Generate conformers from SMILES and add to existing test_params directory (Step 1.5).
	Returns True on success, False on failure."""
	
	# Path to conformator container
	conformator_sif = os.path.join(realm_location, "sif", "conformator_container.sif")
	if not os.path.exists(conformator_sif):
		print(f"ERROR: Conformator container not found at {conformator_sif}")
		return False
	
	# Get residue_types.txt path
	res_types_file = os.path.join(test_params_dir, "residue_types.txt")
	
	# For each unique ligand, extract SMILES and generate 150 conformers
	processed_ligands = set()
	shared_work_dir = tempfile.mkdtemp(dir=tmp_root)
	
	any_success = False
	try:
		for ligand_name, conf_num, chunk, subchunk in ligands_list:
			# Only process each ligand once (skip if already generated conformers)
			if ligand_name in processed_ligands:
				continue
			processed_ligands.add(ligand_name)
			
			print(f"Generating 150 conformers for {ligand_name}...")
			
			# Extract original conformer params (returns text, no file written)
			params_text = extract_conformer_params(
				ligand_name, conf_num, chunk, subchunk, enamine_path, tmp_root=tmp_root
			)
			
			if not params_text:
				print(f"WARNING: Could not extract params for {ligand_name}, skipping conformer generation")
				continue
			
			# The params file doesn't have SMILES, so we need to extract from original SDF
			# For now, use a placeholder approach: we'll create a minimal SMILES representation
			# In production, you'd extract actual SMILES from the Enamine database
			
			# Generate 150 conformers and params
			new_params = generate_conformers_and_params(
				ligand_name, "C",  # Placeholder SMILES - should be extracted from library
				test_params_dir, realm_location, conformator_sif,
				num_conformers=num_conformers, license_key=license_key,
				tmp_root=tmp_root, work_dir=shared_work_dir
			)
			
			if new_params:
				# Add generated params to residue_types.txt
				with open(res_types_file, "a") as fh:
					for params_name in new_params:
						fh.write(f"{params_name}\n")
				any_success = True
		
		return any_success
	
	except Exception as e:
		print(f"ERROR in conformer generation: {e}")
		return False
	
	finally:
		shutil.rmtree(shared_work_dir, ignore_errors=True)


def run_rosetta_discovery(test_params_dir, target_pdb, anchor_residues, 
						  motifs_file, batch_dir, discovery_root,
						  atr, rep, ddg, extra_params=None):
	"""Submit Rosetta discovery jobs via bsub, wait for all to finish, and verify success.
	Returns True only if all jobs were submitted AND completed with exit code 0."""
	import re
	
	job_ids = []
	
	for residue in anchor_residues.split(","):
		residue = residue.strip()
		residue_dir = os.path.join(batch_dir, residue)
		os.makedirs(residue_dir, exist_ok=True)
		
		rosetta_script = os.path.join(discovery_root, "func", "discovery", "run_ligand_discovery_search.py")
		extra_flag = f"--extra-args-file {extra_params}" if extra_params else ""
		# Build Rosetta command
		job_cmd = (
			f"bsub -q long -W 96:00 -n 8 -R \"rusage[mem=10000]\" "
			f"\"python {rosetta_script} {target_pdb} {residue} {motifs_file} "
			f"{test_params_dir} {discovery_root} {atr} {rep} {ddg} {extra_flag}\""
		)
		
		print("Submitting Rosetta discovery job:", job_cmd)
		result = subprocess.run(job_cmd, shell=True, capture_output=True, text=True, cwd=residue_dir)
		if result.stdout:
			print(result.stdout, end="")
		if result.stderr:
			print(result.stderr, file=sys.stderr, end="")
		
		if result.returncode != 0:
			print(f"ERROR: bsub submission failed for residue {residue} (returncode={result.returncode})")
			return False
		
		# Parse job ID from bsub output: "Job <12345> is submitted to queue <long>."
		match = re.search(r"Job <(\d+)>", result.stdout)
		if match:
			job_id = match.group(1)
			job_ids.append(job_id)
			print(f"  Submitted job {job_id} for residue {residue}")
		else:
			print(f"ERROR: Could not parse job ID from bsub output for residue {residue}")
			return False
	
	if not job_ids:
		print("ERROR: No jobs were submitted")
		return False
	
	# Wait for all jobs to finish (ended = done or exited, regardless of status)
	print(f"Waiting for {len(job_ids)} job(s) to complete: {', '.join(job_ids)}")
	ended_conditions = " && ".join(f"ended({jid})" for jid in job_ids)
	wait_cmd = f"bwait -w '{ended_conditions}'"
	wait_result = subprocess.run(wait_cmd, shell=True, capture_output=True, text=True)
	if wait_result.returncode != 0:
		print(f"WARNING: bwait returned non-zero (some jobs may have failed)")
	
	# Check exit status of each job
	all_success = True
	for job_id in job_ids:
		# bjobs -o 'stat exit_code' -noheader <job_id>
		check_cmd = f"bjobs -o 'stat exit_code' -noheader {job_id} 2>/dev/null"
		check_result = subprocess.run(check_cmd, shell=True, capture_output=True, text=True)
		status_line = check_result.stdout.strip()
		print(f"  Job {job_id} status: {status_line}")
		
		if "DONE" not in status_line:
			all_success = False
			print(f"ERROR: Job {job_id} did not complete successfully: {status_line}")
	
	return all_success


def score_placements(work_dir, weights_file="", min_motif_ratio=0.0):
	"""Score placed ligands by parsing placement PDBs and extracting Rosetta metrics.
	Returns dict of {ligand_name: total_weighted_score}."""
	scores = {}
	weights = load_score_weights(weights_file) if weights_file else {}

	for root, dirs, files in os.walk(work_dir):
		for file in files:
			if not (file.endswith(".pdb") and "placed" in file and "minipose" not in file):
				continue

			pdb_path = os.path.join(root, file)
			parsed = parse_placement_scores(pdb_path)
			if not parsed:
				continue

			# Filter by minimum real motif ratio
			if parsed["real_motif_ratio"] < min_motif_ratio:
				continue

			total, _ = compute_weighted_total(parsed, weights)
			lig_name = file.replace(".pdb", "")
			scores[lig_name] = total

	return scores


def _deduplicate_heap():
	"""Remove duplicate base_ligand entries from the heap, keeping only the best score.
	Must be called while holding heap_lock."""
	seen = {}
	deduped = []
	count_removed = 0
	for neg_score, ligand_name, data in ligand_scores_heap:
		base = "_".join(ligand_name.split("_")[:-1])
		if base in seen:
			# Duplicate — keep the entry with the better (more negative) score
			existing_neg, existing_name, existing_data = seen[base]
			if neg_score < existing_neg:  # more negative = better score
				# New entry is better, replace
				_write_error(f"DEDUP_HEAP: Replacing {existing_name} (score={-existing_neg:.4f}) "
				             f"with {ligand_name} (score={-neg_score:.4f})")
				seen[base] = (neg_score, ligand_name, data)
				count_removed += 1
			else:
				_write_error(f"DEDUP_HEAP: Discarding duplicate {ligand_name} (score={-neg_score:.4f}), "
				             f"keeping {existing_name} (score={-existing_neg:.4f})")
				count_removed += 1
		else:
			seen[base] = (neg_score, ligand_name, data)
	for entry in seen.values():
		deduped.append(entry)
	ligand_scores_heap.clear()
	for entry in deduped:
		heapq.heappush(ligand_scores_heap, entry)
	if count_removed > 0:
		_write_error(f"DEDUP_HEAP: Removed {count_removed} duplicate entries, heap size now {len(ligand_scores_heap)}")


def update_top_ligands_heap(scores_dict, max_ligands=1000):
	"""Thread-safe update of top ligands heap.
	Duplicates (same base ligand) are logged and discarded; the best-scored entry is kept."""
	with heap_lock:
		# First, deduplicate the existing heap in case duplicates snuck in
		_deduplicate_heap()
		
		for ligand_name, score in scores_dict.items():
			# Deduplicate: only keep first occurrence of each ligand (no conformer duplicates)
			base_ligand = "_".join(ligand_name.split("_")[:-1])  # Remove conf number
			
			with ligand_seen_lock:
				if base_ligand in ligand_seen:
					_write_error(f"DUP_INSERT: Skipping duplicate ligand {ligand_name} "
					             f"(score={score:.4f}, base={base_ligand} already in heap)")
					continue
				ligand_seen.add(base_ligand)
			
			# Use negative score for max-heap behavior
			entry = (-score, ligand_name, {"score": score})
			
			if len(ligand_scores_heap) < max_ligands:
				heapq.heappush(ligand_scores_heap, entry)
			elif score > -ligand_scores_heap[0][0]:  # Better than worst in heap
				heapq.heapreplace(ligand_scores_heap, entry)


def process_batch(batch_list, target_pdb, anchor_residues, motifs_file,
				  output_dir, realm_location, enamine_path, atr, rep, ddg,
				  batch_id, max_ligands, num_conformers=150, license_key="", tmp_root=None, extra_params=""):
	"""Process a batch of ligands through the entire pipeline.
	Each step is gated: the next step only runs if the previous step succeeded."""
	print(f"\n=== Processing batch {batch_id} ({len(batch_list)} ligands) ===")
	
	# Work directly in tmp_root — no batch subdirectories
	batch_dir = Path(tmp_root) / f"batch_{batch_id}" if tmp_root else  f"batch_{batch_id}"
	os.makedirs(batch_dir, exist_ok=True)
	
	try:
		# Step 1: Extract and prepare conformer params from Enamine library
		print("Step 1: Extracting library conformer params...")
		test_params_dir = create_test_params_dir(
			batch_list, batch_dir, realm_location, enamine_path, tmp_root=tmp_root
		)
		if not test_params_dir:
			msg = f"ERROR: Step 1 (extract conformer params) failed for batch {batch_id}"
			print(msg)
			_write_error(msg)
			return False
		
		# Step 2: Run Rosetta discovery
		print("Step 2: Running Rosetta discovery search...")
		if not run_rosetta_discovery(
			test_params_dir, target_pdb, anchor_residues, motifs_file,
			batch_dir, realm_location, atr, rep, ddg, extra_params=extra_params
		):
			msg = f"ERROR: Step 2 (Rosetta discovery) failed for batch {batch_id}"
			print(msg)
			_write_error(msg)
			return False
		
		# Step 3: Score placements
		print("Step 3: Scoring placements...")
		scores = score_placements(batch_dir, weights_file=os.path.join(realm_location, "func", "discovery", "score_weights.json"))
		if not scores:
			msg = f"ERROR: Step 3 (scoring placements) produced no scores for batch {batch_id}"
			print(msg)
			_write_error(msg)
			return False
		
		# Step 4: Generate 150 conformers and parameterize them
		print("Step 4: Generating 150 conformers from SMILES...")
		if not generate_and_add_conformers_to_test_params(
			test_params_dir, batch_list, realm_location, enamine_path,
			num_conformers=num_conformers, license_key=license_key, tmp_root=tmp_root
		):
			msg = f"ERROR: Step 4 (conformer generation) failed for batch {batch_id}"
			print(msg)
			_write_error(msg)
			return False
		
		# Step 5: Rerun Rosetta discovery with new conformers to get scores for all conformers of each ligand
		print("Step 5: Rerunning Rosetta discovery with new conformers...")	
		if not run_rosetta_discovery(
			test_params_dir, target_pdb, anchor_residues, motifs_file,
			batch_dir, realm_location, atr, rep, ddg, extra_params=extra_params
		):
			msg = f"ERROR: Step 5 (Rosetta discovery rerun) failed for batch {batch_id}"
			print(msg)
			_write_error(msg)
			return False
		
		# Step 6: Update top ligands heap
		print(f"Step 6: Updating heap with {len(scores)} scored ligands...")
		if scores:
			update_top_ligands_heap(scores, max_ligands)
		
		print(f"Batch {batch_id} complete. Current heap size: {len(ligand_scores_heap)}")
		return True
	
	except Exception as e:
		msg = f"Batch {batch_id} failed: {e}"
		print(f"ERROR: {msg}")
		_write_error(msg)
		return False


def write_final_results(output_file, max_ligands=1000):
	"""Write final top ligands to output file."""
	with heap_lock:
		sorted_results = sorted(ligand_scores_heap, reverse=True)
	
	os.makedirs(os.path.dirname(output_file) if os.path.dirname(output_file) else ".", exist_ok=True)
	
	with open(output_file, "w") as fh:
		fh.write("rank,ligand_name,score\n")
		for rank, (neg_score, ligand_name, data) in enumerate(sorted_results, 1):
			score = -neg_score
			fh.write(f"{rank},{ligand_name},{score:.6f}\n")
	
	print(f"\nFinal results written to {output_file}")
	print(f"Total unique ligands: {len(sorted_results)}")
	return output_file


def main():
	parser = argparse.ArgumentParser(
		description="Unified discovery controller: conformers -> discovery -> scoring -> top ligands"
	)
	parser.add_argument("--shapedb-list", required=True, 
					   help="Path to shapedb results list (score,ligand_conf,chunk,subchunk)")
	parser.add_argument("--target-pdb", required=True, help="Target PDB file")
	parser.add_argument("--anchor-residues", required=True, help="Anchor residues (e.g., 79 or 11,79,55)")
	parser.add_argument("--motifs-file", required=True, help="Motifs file path")
	parser.add_argument("--output-dir", required=True, help="Output directory")
	parser.add_argument("--realm-location", default="/pi/summer.thyme-umw/Ji_rosetta_discovery",
					   help="Realm root directory")
	parser.add_argument("--enamine-path", default="/pi/summer.thyme-umw/enamine-REAL-2.6billion",
					   help="Path to Enamine library")
	parser.add_argument("--workers", type=int, default=4, help="Number of parallel workers")
	parser.add_argument("--max-ligands", type=int, default=1000, help="Max ligands to keep")
	parser.add_argument("--batch-size", type=int, default=100, help="Ligands per batch")
	parser.add_argument("--top-hits", type=int, default=0,
					   help="Only process top N hits from the shapedb list (0 = process all)")
	parser.add_argument("--num-conformers", type=int, default=150, 
					   help="Number of conformers to generate per ligand (default: 150)")
	parser.add_argument("--license-key", default="", 
					   help="Conformator license key string (optional)")
	parser.add_argument("--conformator-license", default="",
					   help="Path to file containing the conformator license key (read first line)")
	parser.add_argument("--extra-params", default="",
					   help="Path to file with additional Rosetta discovery arguments")
	parser.add_argument("--atr", type=float, default=-2.0, help="fa_atr cutoff")
	parser.add_argument("--rep", type=float, default=150.0, help="fa_rep cutoff")
	parser.add_argument("--ddg", type=float, default=-9.0, help="ddg cutoff")
	
	args = parser.parse_args()
	
	# Resolve license key: --conformator-license file takes precedence over --license-key
	license_key = args.license_key
	if args.conformator_license:
		try:
			with open(args.conformator_license, "r") as fh:
				license_key = fh.readline().strip()
			print(f"Read license key from {args.conformator_license}")
		except IOError as e:
			print(f"WARNING: Cannot read license file {args.conformator_license}: {e}")
	
	# Setup output directory and logging
	os.makedirs(args.output_dir, exist_ok=True)
	global ERROR_LOG
	ERROR_LOG = open(os.path.join(args.output_dir, "error_log.txt"), "w")
	
	# Create a temp root under realm_location instead of system /tmp
	tmp_root = os.path.join(args.realm_location, "tmp")
	os.makedirs(tmp_root, exist_ok=True)

	print(f"Using temp directory: {tmp_root}")
	
	# Read shapedb list and group into batches
	batches = []
	current_batch = []
	total_entries = 0
	top_hits = args.top_hits
	
	print("Reading shapedb list...")
	with open(args.shapedb_list, "r") as fh:
		for line in fh:
			line = line.strip()
			if not line:
				continue
			
			# Stop reading if --top-hits limit reached
			if top_hits > 0 and total_entries >= top_hits:
				break
			
			fields = line.split(",")
			if len(fields) < 4:
				print(f"Skipping malformed line: {line}")
				continue
			
			score, ligand_conf, chunk, subchunk = fields[:4]
			ligand_name = "_".join(ligand_conf.split("_")[:-1])
			conf_num = ligand_conf.split("_")[-1]
			
			current_batch.append((ligand_name, int(conf_num), chunk, subchunk))
			total_entries += 1
			
			if len(current_batch) >= args.batch_size:
				batches.append(current_batch)
				current_batch = []
	
	if current_batch:
		batches.append(current_batch)
	
	if top_hits > 0:
		print(f"Loaded {len(batches)} batches from shapedb list (limited to top {total_entries} hits)")
	else:
		print(f"Loaded {len(batches)} batches from shapedb list ({total_entries} total entries)")
	
	# Process batches in parallel
	futures = []
	with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as executor:
		for batch_id, batch in enumerate(batches):
			future = executor.submit(
				process_batch,
				batch,
				args.target_pdb,
				args.anchor_residues,
				args.motifs_file,
				args.output_dir,
				args.realm_location,
				args.enamine_path,
				args.atr,
				args.rep,
				args.ddg,
				batch_id,
				args.max_ligands,
				args.num_conformers,
				license_key,
				tmp_root,
				args.extra_params
			)
			futures.append(future)
	
	# Wait for completion
	completed = 0
	for future in concurrent.futures.as_completed(futures):
		try:
			result = future.result()
			if result:
				completed += 1
		except Exception as e:
			print(f"Batch failed: {e}")
	
	print(f"\nCompleted {completed}/{len(futures)} batches")
	
	# Write final results
	output_file = os.path.join(args.output_dir, "top_ligands.csv")
	write_final_results(output_file, args.max_ligands)
	
	if ERROR_LOG:
		ERROR_LOG.close()
	
	print("\nDiscovery pipeline complete!")


if __name__ == "__main__":
	main()
