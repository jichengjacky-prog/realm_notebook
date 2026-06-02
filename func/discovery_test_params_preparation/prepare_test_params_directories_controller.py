#the purpose of this script is to run instances of the prepare_test_params_directories_sub.py script on up to 100 conformers per job, to control and parallelize the process
#this script needs inputs of a location to make the test_params directories and a shapedb list formatted like the following:
"""
-0.489485800266,PV-005633531035_4,29811,1
-0.492682218552,Z4448503755_1,01035,7
-0.492682218552,Z3789042056_8,01035,7
-0.498856127262,PV-004964373284_3,38709,7
-0.498856127262,PV-004964373283_11,38709,7
-0.498872220516,PV-006084859944_3,31726,1
-0.499529778957,PV-005674048644_3,25233,5
-0.50066614151,Z4211293263_2,01366,2
"""

#to avoid overflowing a single directory, directorie will be nested where 100 test_params directories will be made in a single sub-location beneath the working location for all directories
#i.e. if there are 3 million ligands in the input file, there will be 300 top-level directories, each having 100 sub-directories of test_params, with 100 ligand conformers per test param directory
#(300 top-level directories * 100 lower directories per top-level * 100 conformers per lower level directory = 3,000,000)


#imports
import argparse
import concurrent.futures
import os
import subprocess
import sys


def run_cmd(cmd, cwd=None):
	result = subprocess.run(cmd, shell=True, capture_output=True, text=True, cwd=cwd)
	if result.stdout:
		print(result.stdout, end="")
	if result.stderr:
		print(result.stderr, file=sys.stderr, end="")
	return result.returncode


def submit_job(batch_dir, script_path):
	cmd = f'bsub -q long -R "rusage[mem=1024]" -n 1 -W 8:00 -u "" "python {script_path}"'
	return run_cmd(cmd, cwd=batch_dir)


def increment_directory_indices(top_level, sub):
	sub += 1
	if sub == 100:
		sub = 0
		top_level += 1
	return top_level, sub


def create_batch_dir(base_dir, top_level, sub):
	path = os.path.join(base_dir, str(top_level), str(sub))
	os.makedirs(path, exist_ok=True)
	return path


# parse CLI arguments
parser = argparse.ArgumentParser(
	description="Create nested test_params directories in batches from a shapedb list."
)
parser.add_argument("working_location", help="Base directory for generated test_params directories")
parser.add_argument("master_list", help="Path to the shapedb-style ligand conformer list file")
parser.add_argument(
	"--workers",
	type=int,
	default=4,
	help="Number of threads used to submit jobs in parallel"
)
args = parser.parse_args()

#get the working location and initial list file
working_location = os.path.abspath(args.working_location)
master_list = os.path.abspath(args.master_list)
script_path = os.path.abspath(os.path.join(os.path.dirname(__file__), "prepare_test_params_directories_sub.py"))

#declare counters for the top level directories and sub directories
top_level_dirs = 0
sub_dirs = 0

#declare a working list to hold batches of up to 100 ligands to make the lists and jobs
small_confs_list = []

# job futures tracked for parallel submission
futures = []

with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as executor:
	batch_dir = create_batch_dir(working_location, top_level_dirs, sub_dirs)
	conf_list_path = os.path.join(batch_dir, "conf_list.csv")
	write_file = open(conf_list_path, "w")

	with open(master_list, "r") as master_list_file:
		for line in master_list_file:
			line = line.rstrip("\n")
			if not line:
				continue

			#extract the ligand with conformer to ensure it is not in the working small_confs list
			fields = line.split(",")
			if len(fields) < 2:
				print(f"WARNING: malformed line skipped: {line}")
				continue
			ligconf = fields[1]

			#continue to avoid repeating, although this should generally be impossible. having the same ligand in the same call to rosetta twice could break things
			if ligconf in small_confs_list:
				print("WARNING: Encountered repeat of ligand conformer of " + ligconf)
				continue

			#assuming we can move forward with the ligand, add it to the small confs list and write it to the write file
			write_file.write(line + "\n")
			small_confs_list.append(ligconf)

			#if the size of the small confs list reaches 100, prepare the run the sub script on this directory and then move on to the next
			if len(small_confs_list) == 100:
				write_file.close()
				futures.append(executor.submit(submit_job, batch_dir, script_path))
				small_confs_list = []

				#increment directories for the next batch
				top_level_dirs, sub_dirs = increment_directory_indices(top_level_dirs, sub_dirs)
				batch_dir = create_batch_dir(working_location, top_level_dirs, sub_dirs)
				conf_list_path = os.path.join(batch_dir, "conf_list.csv")
				write_file = open(conf_list_path, "w")

	#end behavior for final list
	write_file.close()
	if small_confs_list:
		futures.append(executor.submit(submit_job, batch_dir, script_path))

	# wait for all submissions to complete
	for future in concurrent.futures.as_completed(futures):
		exit_code = future.result()
		if exit_code != 0:
			print(f"Job submission failed with exit code {exit_code}", file=sys.stderr)

print(f"Submitted {len(futures)} batch jobs with up to {args.max_workers} concurrent workers.")