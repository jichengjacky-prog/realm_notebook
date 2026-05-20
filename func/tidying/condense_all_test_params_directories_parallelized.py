#quick script to parallelize condensing all test_params directories from a working location for archiving

#imports 
import os,sys

working_location = os.getcwd()

dire_counter = 0

for r,d,f in os.walk(working_location):
	for dire in d:
		if dire == "test_params":

			dire_counter = dire_counter + 1

			print(dire_counter)

			#move to the location
			os.chdir(r)

			#pause 0.1 seconds to help avoid overload
			#os.system("sleep 0.1")

			command = "bsub -q long -W 1:00 -u \"\" -R \"rusage[mem=5000]\" \"tar -czf test_params.tar.gz test_params && rm -drf test_params \""

			#submit a bsub job that runs the nnsearch python script
			print(command)
			os.system(command)

			#adding 500 job throttle
			#bsub job throttle to make sure we do not exceed our local limit
			#write the length of the bjobs queue to this current location
			os.system("bjobs | wc -l > bjobs_length.txt")
			job_count = 0
			with open("bjobs_length.txt") as f:
				job_count = int(f.read().strip())
			while job_count > 300:
				#sleep for 1 second to not overburden the system
				os.system("sleep 1")
				os.system("bjobs | wc -l > bjobs_length.txt")
				with open("bjobs_length.txt") as f:
					job_count = int(f.read().strip())
			#remove the length file to avoid clutter
			os.system("rm bjobs_length.txt")