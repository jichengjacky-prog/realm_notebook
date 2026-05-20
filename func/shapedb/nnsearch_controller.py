#the purpose of this script is to be an automated controller script that runs runn_nnsearch_hpc.py on the entire enamine library for a given ligand shape
#the user needs to at least give an aligned ligand molecule shape (mol2 or sdf), and can optionally give a path to a location to put the results in

#imports 
import os,sys
import subprocess
from pathlib import Path

def run_cmd(cmd):
	result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
	if result.stdout:
		print(result.stdout, end="")
	if result.stderr:
		print(result.stderr, file=sys.stderr, end="")
	return result.returncode

#get the ligand file
target_molecule_file = sys.argv[1]

#create the result output in a specified location, or use location where this was called from otherwise
working_location = "."
if len(sys.argv) > 2:
	working_location = sys.argv[2]

realm_location = "/pi/summer.thyme-umw/Ji_rosetta_discovery"
if len(sys.argv) > 3:
	realm_location = sys.argv[3]

#iterate over every chunk, and call the search script on each subchunk within the chunk
for i in range(0,53085):
#for i in range(5001,53085):
	#build the chunk string
	chunk_str = str(i)

	#append leading zeroes until the string is 5 digits long
	while len(chunk_str) < 5:
		chunk_str = "0" + chunk_str

	
	#submit a bsub job that runs the nnsearch python script
	queue_cmd =["bsub -q \"long short\" -n 1 -W 1:00  -u \"\" -R \"rusage[mem=5000]\" " ,
			 "python",
			 realm_location + "/func/shapedb/run_nnsearch_hpc.py ",
             chunk_str,
			 target_molecule_file,
			 working_location,
			 realm_location]
	
	print(" ".join(queue_cmd))
	run_cmd(" ".join(queue_cmd))
	##throttle the job submission to avoid overloading the system, we will check the length of the bjobs queue and if it is above 100, we will wait until it is below 100 before submitting the next job
	os.system("sleep 0.1s")



	# #adding 500 job throttle
	# #bsub job throttle to make sure we do not exceed our local limit
	# #write the length of the bjobs queue to this current location
	# os.system("bjobs | wc -l > " + working_location + "/bjobs_length.txt")
	# job_count = 0
	# with open(working_location + "/bjobs_length.txt") as f:
	# 	job_count = int(f.read().strip())
	# while job_count % 100 == 0 and job_count != 0:
	# 	#sleep for 1 second to not overburden the system
	# 	os.system("sleep 10s")
	# 	os.system("bjobs | wc -l > " + working_location + "/bjobs_length.txt")
	# 	with open(working_location+"bjobs_length.txt") as f:
	# 		job_count = int(f.read().strip())
	# #remove the length file to avoid clutter
	# os.system("rm " + working_location + "/bjobs_length.txt")