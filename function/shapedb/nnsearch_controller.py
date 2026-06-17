#the purpose of this script is to be an automated controller script that runs run_nnsearch_hpc.py on the entire enamine library for a given ligand shape
#the user needs to at least give an aligned ligand molecule shape (mol2 or sdf), and can optionally give a path to a location to put the results in

#imports 
import argparse
import heapq
import os,sys
import re
import shutil
import subprocess
import time
from pathlib import Path

ERROR_LOG = None


def _write_error(msg):
	"""Write an error message to both stderr and the error log file."""
	print(msg, file=sys.stderr)
	if ERROR_LOG:
		ERROR_LOG.write(msg + "\n")
		ERROR_LOG.flush()


def _deduplicate_heap(heap, ligand_seen):
	"""Remove duplicate base-ligand entries from the heap, keeping only the best score.
	Returns (deduped_heap, updated_ligand_seen, count_removed)."""
	if not heap:
		return heap, ligand_seen, 0
	
	seen = {}  # base_ligand -> (score, conf_name, chunk_str, subchunk_str)
	count_removed = 0
	for entry in heap:
		conf_score, conf_name, chunk_str, subchunk_str = entry
		base = "_".join(conf_name.split("_")[:-1])
		if base in seen:
			existing_score, existing_name, existing_chunk, existing_sub = seen[base]
			if conf_score > existing_score:  # higher (more negative) = better shape similarity
				_write_error(f"DEDUP_HEAP: Replacing {existing_name} (score={existing_score:.4f}) "
				             f"with {conf_name} (score={conf_score:.4f})")
				seen[base] = (conf_score, conf_name, chunk_str, subchunk_str)
			else:
				_write_error(f"DEDUP_HEAP: Discarding duplicate {conf_name} (score={conf_score:.4f}), "
				             f"keeping {existing_name} (score={existing_score:.4f})")
			count_removed += 1
		else:
			seen[base] = (conf_score, conf_name, chunk_str, subchunk_str)
	
	deduped = list(seen.values())
	heapq.heapify(deduped)
	if count_removed > 0:
		_write_error(f"DEDUP_HEAP: Removed {count_removed} duplicate entries, heap size now {len(deduped)}")
	
	# Update ligand_seen to reflect the deduplicated set
	ligand_seen.clear()
	ligand_seen.update(seen.keys())
	
	return deduped, ligand_seen, count_removed

def run_cmd(cmd):
	result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
	if result.stdout:
		print(result.stdout, end="")
	if result.stderr:
		print(result.stderr, file=sys.stderr, end="")
	return result.returncode

def chunk_to_path(chunk):
	chunk_str = str(chunk).zfill(5)
	superchunk_str = str(int(chunk_str[:3]))
	return superchunk_str, chunk_str

def load_chunk_into_heap(max_ligands, heap, working_location, chunk, ligand_seen):
	"""Read the nn_filtered.txt files for a single chunk and merge into the running heap.
	Duplicates (same base ligand across conformers/chunks) are logged and discarded;
	only the best-scored entry per base ligand is kept."""
	superchunk_str, chunk_str = chunk_to_path(chunk)
	for subchunk in range(10):
		path = os.path.join(working_location, superchunk_str, chunk_str,
							f'{chunk_str}_{subchunk}_nn_filtered.txt')
		if not os.path.isfile(path):
			print(f'  missing: {path}')
			continue

		try:
			with open(path, 'r') as fh:
				for line in fh:
					line = line.strip()
					if not line:
						continue
					fields = line.split()
					if len(fields) < 2:
						continue
					conf_name = fields[0]
					try:
						conf_score = float(fields[1]) * -1
					except ValueError:
						continue
					
					# Deduplicate by base ligand name (strip conformer number)
					base_ligand = "_".join(conf_name.split("_")[:-1])
					if base_ligand in ligand_seen:
						_write_error(f"DUP_INSERT: Skipping duplicate ligand {conf_name} "
						             f"(score={conf_score:.4f}, base={base_ligand} already in heap)")
						continue

					entry = (conf_score, conf_name, chunk_str, str(subchunk))
					if len(heap) < max_ligands:
						heap.append(entry)
						ligand_seen.add(base_ligand)
						if len(heap) == max_ligands:
							heapq.heapify(heap)
					else:
						if conf_score > heap[0][0]:
							# Remove the worst entry's base from seen set
							worst_base = "_".join(heap[0][1].split("_")[:-1])
							ligand_seen.discard(worst_base)
							heapq.heapreplace(heap, entry)
							ligand_seen.add(base_ligand)
		except Exception as exc:
			print(f'  error reading {path}: {exc}')


def write_results(heap, working_location, max_ligands, min_chunk, max_chunk):
	output_dir = os.path.join(working_location, 'combined_results')
	os.makedirs(output_dir, exist_ok=True)
	output_file = os.path.join(
		output_dir,
		f'combined_best_{max_ligands}_chunks_{str(min_chunk).zfill(5)}_{str(max_chunk).zfill(5)}.txt'
	)
	sorted_results = sorted(heap, reverse=True)
	with open(output_file, 'w') as fh:
		for conf_score, conf_name, chunk_str, subchunk_str in sorted_results:
			fh.write(f'{conf_score},{conf_name},{chunk_str},{subchunk_str}\n')
	print(f'Wrote {len(sorted_results)} combined entries to {output_file}')
	return output_file


def save_heap_checkpoint(heap, working_location, last_chunk):
	"""Save the current heap to a checkpoint file after each chunk merge.
	This protects against data loss if the controller is interrupted."""
	checkpoint_dir = os.path.join(working_location, 'combined_results')
	os.makedirs(checkpoint_dir, exist_ok=True)
	checkpoint_file = os.path.join(checkpoint_dir, 'heap_checkpoint.txt')
	tmp_file = checkpoint_file + '.tmp'

	# sort descending so the best ligands are at the top when inspecting manually
	sorted_results = sorted(heap, reverse=True)
	with open(tmp_file, 'w') as fh:
		fh.write(f'# last_merged_chunk={str(last_chunk).zfill(5)} heap_size={len(sorted_results)}\n')
		for conf_score, conf_name, chunk_str, subchunk_str in sorted_results:
			fh.write(f'{conf_score},{conf_name},{chunk_str},{subchunk_str}\n')

	os.replace(tmp_file, checkpoint_file)  # atomic rename
	print(f'Checkpoint saved: {len(sorted_results)} entries to {checkpoint_file}')


def submit_chunk_job(chunk, target_molecule_file, working_location, realm_location):
	"""Submit a bsub job for a single chunk. Returns the LSF job ID."""
	chunk_str = str(chunk).zfill(5)
	superchunk_str, _ = chunk_to_path(chunk)

	# ensure the chunk directory exists before submitting the job
	chunk_dir = os.path.join(working_location, superchunk_str, chunk_str)
	os.makedirs(chunk_dir, exist_ok=True)

	queue_cmd = [
		f'bsub -q "short" -n 1 -W 1:00 -u "" -R "rusage[mem=6000]"',
		"python",
		f'{realm_location}/func/shapedb/run_nnsearch_hpc.py',
		chunk_str,
		target_molecule_file,
		working_location,
		realm_location,
	]
	cmd_str = " ".join(queue_cmd)
	print(cmd_str)
	result = subprocess.run(cmd_str, shell=True, capture_output=True, text=True)
	if result.stdout:
		print(result.stdout, end="")
	if result.stderr:
		print(result.stderr, file=sys.stderr, end="")

	# parse job ID from bsub output: "Job <12345> is submitted..."
	match = re.search(r'Job <(\d+)>', result.stdout)
	if match:
		return match.group(1)
	return None


def get_done_job_ids(pending_job_ids, batch_size=500):
	"""Return the subset of job_ids that are DONE/EXIT or no longer tracked by LSF.
	Queries bjobs in batches to avoid command-line length limits with many thousands of IDs."""
	if not pending_job_ids:
		return []

	still_tracked = set()
	for i in range(0, len(pending_job_ids), batch_size):
		batch = pending_job_ids[i:i + batch_size]
		result = subprocess.run(
			f'bjobs -o "jobid stat" -noheader {" ".join(batch)}',
			shell=True, capture_output=True, text=True
		)
		# bjobs outputs nothing for jobs that no longer exist (cleaned up = done)
		for line in result.stdout.strip().split('\n'):
			if not line.strip():
				continue
			parts = line.split()
			if len(parts) >= 2:
				jid, status = parts[0], parts[1]
				if 'DONE' not in status and 'EXIT' not in status:
					still_tracked.add(jid)

	return [jid for jid in pending_job_ids if jid not in still_tracked]


def process_finished_chunk(chunk, max_ligands, combined_heap, working_location,
						   superchunk_remaining, ligand_seen, checkpoint=True):
	"""Load chunk results into heap, delete the chunk directory, and
	delete the superchunk directory once all its chunks in range are done."""
	chunk_str = str(chunk).zfill(5)
	print(f'Loading chunk {chunk_str} results into heap...')
	load_chunk_into_heap(max_ligands, combined_heap, working_location, chunk, ligand_seen)
	print(f'Heap now has {len(combined_heap)} entries (unique ligands: {len(ligand_seen)})')
	if checkpoint:
		save_heap_checkpoint(combined_heap, working_location, chunk)

	# delete the chunk directory to free disk space
	superchunk_str, _ = chunk_to_path(chunk)
	chunk_dir = os.path.join(working_location, superchunk_str, chunk_str)
	if os.path.isdir(chunk_dir):
		print(f'Deleting chunk directory: {chunk_dir}')
		shutil.rmtree(chunk_dir, ignore_errors=True)

	# delete the superchunk directory once all its jobs are done
	superchunk_remaining[superchunk_str] -= 1
	if superchunk_remaining[superchunk_str] == 0:
		superchunk_dir = os.path.join(working_location, superchunk_str)
		if os.path.isdir(superchunk_dir):
			try:
				os.rmdir(superchunk_dir)
				print(f'Deleted superchunk directory: {superchunk_dir}')
			except OSError:
				pass
		del superchunk_remaining[superchunk_str]


def parse_args():
	parser = argparse.ArgumentParser(
		description='Automated controller that runs ShapeDB NNSearch across enamine library chunks '
					'with a dynamic worker pool, keeping the top ligands via a min-heap.'
	)
	parser.add_argument(
		'target_molecule',
		help='Path to the aligned ligand molecule shape file (mol2 or sdf)'
	)
	parser.add_argument(
		'-w', '--working-dir',
		default='.',
		help='Directory for chunk output and combined results (default: .)'
	)
	parser.add_argument(
		'-r', '--realm-dir',
		default='/pi/summer.thyme-umw/Ji_rosetta_discovery',
		help='Path to the Ji_rosetta_discovery realm root'
	)
	parser.add_argument(
		'-n', '--max-ligands',
		type=int, default=100000,
		help='Number of top ligands to keep (default: 100000)'
	)
	parser.add_argument(
		'-j', '--workers',
		type=int, default=100,
		help='Number of concurrent LSF jobs (default: 100)'
	)
	parser.add_argument(
		'--min-chunk',
		type=int, default=0,
		help='Minimum chunk index to process, inclusive (default: 0)'
	)
	parser.add_argument(
		'--max-chunk',
		type=int, default=53084,
		help='Maximum chunk index to process, inclusive (default: 53084)'
	)
	parser.add_argument(
		'--poll-interval',
		type=int, default=5,
		help='Seconds between polling bjobs for completion (default: 5)'
	)
	return parser.parse_args()


def main():
	args = parse_args()

	target_molecule_file = args.target_molecule
	working_location = args.working_dir
	realm_location = args.realm_dir
	max_ligands = args.max_ligands
	num_workers = args.workers
	min_chunk = args.min_chunk
	max_chunk = args.max_chunk

	#track the top N ligands across all chunks with a min-heap
	combined_heap = []
	ligand_seen = set()  # track unique base ligands to prevent duplicates

	# Set up error logging
	global ERROR_LOG
	error_log_path = os.path.join(working_location, 'combined_results', 'dedup_error_log.txt')
	os.makedirs(os.path.join(working_location, 'combined_results'), exist_ok=True)
	ERROR_LOG = open(error_log_path, 'w')
	_write_error(f"=== Dedup error log started at {time.ctime()} ===")

	# count how many chunks per superchunk are in our range, for cleanup tracking
	from collections import defaultdict
	superchunk_remaining = defaultdict(int)
	for c in range(min_chunk, max_chunk + 1):
		sc, _ = chunk_to_path(c)
		superchunk_remaining[sc] += 1

	# --- dynamic worker pool ---
	# pending: dict mapping job_id -> chunk number
	pending = {}
	chunk_iter = iter(range(min_chunk, max_chunk + 1))

	# submit initial batch
	for _ in range(num_workers):
		try:
			chunk = next(chunk_iter)
		except StopIteration:
			break
		job_id = submit_chunk_job(chunk, target_molecule_file, working_location, realm_location)
		if job_id:
			pending[job_id] = chunk
			print(f'  submitted chunk {str(chunk).zfill(5)} as job {job_id} ({len(pending)} active)')

	print(f'Initial batch submitted: {len(pending)} jobs running')

	# poll until all jobs done and all chunks processed
	has_more_chunks = True
	while pending or has_more_chunks:
		time.sleep(args.poll_interval)

		# --- check which jobs finished ---
		done_job_ids = get_done_job_ids(list(pending.keys()))

		# --- process ALL finished jobs first ---
		for job_id in done_job_ids:
			chunk = pending.pop(job_id)
			print(f'Job {job_id} (chunk {str(chunk).zfill(5)}) finished')
			process_finished_chunk(chunk, max_ligands, combined_heap, working_location,
								   superchunk_remaining, ligand_seen,
								   checkpoint=(chunk % 100 == 0))

		# --- bulk-fill the pool back to num_workers ---
		#  instead of 1-in-1-out, submit as many as needed to max out the pool
		slots = num_workers - len(pending)
		new_submissions = 0
		for _ in range(slots):
			try:
				next_chunk = next(chunk_iter)
			except StopIteration:
				has_more_chunks = False
				break
			job_id = submit_chunk_job(next_chunk, target_molecule_file, working_location, realm_location)
			if job_id:
				pending[job_id] = next_chunk
				new_submissions += 1

		if done_job_ids or new_submissions:
			print(f'  pool: {len(pending)}/{num_workers} active, '
				  f'{new_submissions} new, heap: {len(combined_heap)}')
		elif not has_more_chunks and not pending:
			pass  # all done, will exit loop
		else:
			# heartbeat every ~2.5 min (poll_interval=5 * 30 cycles)
			pass

	#write the final combined results
	write_results(combined_heap, working_location, max_ligands, min_chunk, max_chunk)
	
	if ERROR_LOG:
		_write_error(f"=== Dedup error log ended at {time.ctime()} ===")
		ERROR_LOG.close()


if __name__ == '__main__':
	main()

