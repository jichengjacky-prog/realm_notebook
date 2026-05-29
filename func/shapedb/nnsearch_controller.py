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

def load_chunk_into_heap(max_ligands, heap, working_location, chunk):
	"""Read the nn_filtered.txt files for a single chunk and merge into the running heap."""
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

					entry = (conf_score, conf_name, chunk_str, str(subchunk))
					if len(heap) < max_ligands:
						heap.append(entry)
						if len(heap) == max_ligands:
							heapq.heapify(heap)
					else:
						if conf_score > heap[0][0]:
							heapq.heapreplace(heap, entry)
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


def submit_chunk_job(chunk, target_molecule_file, working_location, realm_location):
	"""Submit a bsub job for a single chunk. Returns the LSF job ID."""
	chunk_str = str(chunk).zfill(5)
	superchunk_str, _ = chunk_to_path(chunk)

	# ensure the chunk directory exists before submitting the job
	chunk_dir = os.path.join(working_location, superchunk_str, chunk_str)
	os.makedirs(chunk_dir, exist_ok=True)

	queue_cmd = [
		f'bsub -q "long short" -n 1 -W 1:00 -u "" -R "rusage[mem=6000]"',
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


def get_done_job_ids(pending_job_ids):
	"""Return the subset of job_ids that are DONE/EXIT or no longer tracked by LSF.
	Queries all IDs in a single bjobs call for efficiency."""
	if not pending_job_ids:
		return []

	result = subprocess.run(
		f'bjobs -o "jobid stat" -noheader {" ".join(pending_job_ids)}',
		shell=True, capture_output=True, text=True
	)
	# bjobs outputs nothing for jobs that no longer exist (cleaned up = done)
	still_tracked = set()
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
						   superchunk_remaining):
	"""Load chunk results into heap, delete the chunk directory, and
	delete the superchunk directory once all its chunks in range are done."""
	chunk_str = str(chunk).zfill(5)
	print(f'Loading chunk {chunk_str} results into heap...')
	load_chunk_into_heap(max_ligands, combined_heap, working_location, chunk)
	print(f'Heap now has {len(combined_heap)} entries')

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
		type=int, default=1000,
		help='Number of top ligands to keep (default: 1000)'
	)
	parser.add_argument(
		'-j', '--workers',
		type=int, default=30,
		help='Number of concurrent LSF jobs (default: 30)'
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
		type=int, default=15,
		help='Seconds between polling bjobs for completion (default: 15)'
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
	while pending:
		time.sleep(args.poll_interval)

		done_job_ids = get_done_job_ids(list(pending.keys()))

		for job_id in done_job_ids:
			chunk = pending.pop(job_id)
			print(f'Job {job_id} (chunk {str(chunk).zfill(5)}) finished')
			process_finished_chunk(chunk, max_ligands, combined_heap, working_location,
								   superchunk_remaining)

			# submit next chunk if any remain
			try:
				next_chunk = next(chunk_iter)
			except StopIteration:
				continue
			new_job_id = submit_chunk_job(next_chunk, target_molecule_file, working_location, realm_location)
			if new_job_id:
				pending[new_job_id] = next_chunk
				print(f'  submitted chunk {str(next_chunk).zfill(5)} as job {new_job_id} ({len(pending)} active)')

		if done_job_ids:
			print(f'Active jobs: {len(pending)}, heap size: {len(combined_heap)}')

	#write the final combined results
	write_results(combined_heap, working_location, max_ligands, min_chunk, max_chunk)


if __name__ == '__main__':
	main()

