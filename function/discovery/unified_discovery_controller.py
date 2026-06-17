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
		# Build Rosetta command with stdout/stderr redirected to files
		stdout_file = os.path.join(residue_dir, f"{residue}_rosetta.out")
		stderr_file = os.path.join(residue_dir, f"{residue}_rosetta.err")
		job_cmd = (
			f"bsub -q long -W 12:00 -n 8 -u \"\" -R \"rusage[mem=1000]\" "
			f"-o {stdout_file} -e {stderr_file} "
			f"\"python {rosetta_script} {target_pdb} {residue} {motifs_file} "
			f"{test_params_dir}  {discovery_root} {atr} {rep} {ddg} {extra_flag}\""
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
	
		# Poll for all jobs to finish with a 60-second interval
	import time
	pending = set(job_ids)
	poll_interval = 60
	print(f"Waiting for {len(pending)} job(s) to complete (polling every {poll_interval}s)...")
	while pending:
		time.sleep(poll_interval)
		still_pending = set()
		for jid in sorted(pending):
			check_cmd = f"bjobs -o 'stat' -noheader {jid} 2>/dev/null"
			check_result = subprocess.run(check_cmd, shell=True, capture_output=True, text=True)
			status = check_result.stdout.strip()
			if not status:
				# Job no longer tracked by LSF — treat as done
				print(f"  Job {jid}: no longer tracked (presumed done)")
				continue
			if "DONE" in status or "EXIT" in status:
				print(f"  Job {jid}: {status}")
				continue
			# Still running or pending
			still_pending.add(jid)
		pending = still_pending
		if pending:
			print(f"  Still waiting for {len(pending)} job(s): {', '.join(sorted(pending))}")
	
	# Check exit status of each job
	all_success = True
	for job_id in job_ids:
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
			if not (file.endswith(".pdb") and "minipose" not in file):
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
	batch_dir = Path(tmp_root) / f"batch_{batch_id:05d}" if tmp_root else  f"batch_{batch_id:05d}"
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
	parser.add_argument("--arg", help="arg file path")
	parser.add_argument("--motifs-file", required=True, help="Motifs file path")
	parser.add_argument("--output-dir", required=True, help="Output directory")
	parser.add_argument("--realm-location", default="/pi/summer.thyme-umw/Ji_rosetta_discovery",
					   help="Realm root directory")
	parser.add_argument("--enamine-path", default="/pi/summer.thyme-umw/enamine-REAL-2.6billion",
					   help="Path to Enamine library")
	parser.add_argument("--workers", type=int, default=4, help="Number of parallel workers")
	parser.add_argument("--max-ligands", type=int, default=1000, help="Max ligands to keep")
	parser.add_argument("--batch-size", type=int, default=10, help="Ligands per batch")
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
