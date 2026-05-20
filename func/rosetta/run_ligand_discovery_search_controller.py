#the purpose of this script is to be a controller that runs the run_ligand_discovery_search.py script on a database of test_params directories to perform rosetta ligand discovery search
#this script requires the following inputs, some of which require knowledge of how Rosetta reads in the target pdb to know the indices of residues as Rosetta perceives it:
#the script needs the following
#1. target pdb
#2. an anchor residue  string of all positions to be used as anchors, which is comma separated if there is more than 1 (i.e. 79 for 1 or 11,79,55,403 for multiple); this residue needs to use rosetta indexing, which needs to be determined before running the pipeline; for each anchor residue, a unique job will be made to improve runtime and avoid the risk of clobbering placement files (since files from the same ligand have a risk of being clobbered)
#3. a motifs file (by default, will use the main 1.6M motif library, but if another is specified, will be used instead)
#4. a location to look down where there are test_params folders. the script will look down at all test_params folders from the given location

#additional arguments that we will treat as mandatory (5-7)
#score cutoff overrides for fa_atr, fa_rep, and ddg

#an option to bring in a text file containing arguments to be appended to the discovery arguments for things like a space fill region or to print the space fill matrix (not recommended forr memory and speed)
#dimensions and center for a space fill matrix (which can speed up discovery by filtering away placements that are not close enough to the desired binding region)

#by default, this script will try to run individual jobs for up to 7 days. if this really needs to be changed, I can make an option or something, but jobs shouldn't be running for that long anyway

#imports
import os,sys

#take in arguments
#target pdb
target_pdb = sys.argv[1]
#anchor residue(s)
anchor_residue_string = sys.argv[2]

#break apart the string by commas into a list
anchor_residue_string_list = anchor_residue_string.split(",")

#motifs file (likely want /pi/summer.thyme-umw/enamine-REAL-2.6billion/FINAL_motifs_list_filtered_2_3_2023.motifs unless you know what you are doing)
motifs_file = sys.argv[3]

#discovery directory root, will look down from here for all test params directories and attempt discovery for each anchor residue at the same level as each test_params directory
discovery_directory_root = sys.argv[4]

#atr, rep, ddg cutoffs
atr = sys.argv[5]
rep = sys.argv[6]
ddg = sys.argv[7]

#determine whether to clobber existing COMPLETE discovery directories (set as binary 1-true and do clobber or 0-false and do not clobber)
#complete is defined by a directory having a compressed placements directory (placements.tar.gz) and raw scores directory (raw_scores.csv)
clobber = int(sys.argv[8])

#if there is an extra args file, take it in, so we can pass it down to all discovery jobs
extra_args_file = ""
if len(sys.argv) >= 10:
	extra_args_file = sys.argv[9]

#look over the discovery directory root to identify test_params directories
for r,d,f in os.walk(discovery_directory_root):
	for dire in d:
		#only look at the test_params directories
		if dire == "test_params":

			#store the root
			tp_root = r

			#for each test_params directory, go over the anchor residue string and prepare to run discovery for each unique residue in the list
			for residue in anchor_residue_string_list:
				#go to the root
				os.chdir(tp_root)

				#determine whether there was a completed run that we do not want to clobber and avoid it
				if clobber == 0:
					#check the directory for if there is a placements.tar.gz file and a raw_scores.csv file. if there are both, do not clobber and simply continue
					has_placements_tar = False
					has_raw_scores = False
					for r2,d2,f2 in os.walk(tp_root + "/" + str(residue)):
						for file2 in f2:
							if file2 == "placements.tar.gz":
								has_placements_tar = True
							if file2 == "raw_scores.csv":
								has_raw_scores = True

					#if the directory has both, do not clobber and simply continue
					if has_placements_tar and has_raw_scores:
						continue


				#clobber the existing directory for a fresh run
				os.system("rm -drf " + str(residue))

				#make a directory for the anchor residue
				os.system("mkdir -p " + str(residue))

				#enter the directory
				os.chdir(str(residue))

				#determine whether to send the next job to the long or large queue, based on the number of jobs running in each queue
				#if there are under 1600 long jobs, send to long
				#else if there are under 900 large jobs, send to large
				#otherwise if both buffered queues are running and full, throttle until space opens up

				#initialize
				queue = "None"

				#add a throttle at 1600 parallel jobs to avoid overflowing the system (1600 so extra jobs can be queued outside of the 1500 allowed to run)
				#bsub job throttle to make sure we do not exceed our local limit
				#write the length of the bjobs queue to this current location
				os.system("bjobs | awk '{print $4}' | grep long | wc -l > bjobs_long_length.txt")
				os.system("bjobs | awk '{print $4}' | grep large | wc -l > bjobs_large_length.txt")
				long_job_count = 0
				large_job_count = 0
				with open("bjobs_long_length.txt") as f:
					long_job_count = int(f.read().strip())
				with open("bjobs_large_length.txt") as f:
					large_job_count = int(f.read().strip())

				#queue assignment
				if long_job_count < 1600:
					queue = "long"

				#changing from 900 to -1 so we don't try to use the large queue anymore
				if long_job_count > 1600 and large_job_count < -1:
					queue = "large"


				#hit the throttle if both queues are stuffed
				#while long_job_count > 1600 and large_job_count > 900 and queue == "None":
				while True:
					
					os.system("bjobs | awk '{print $4}' | grep long | wc -l > bjobs_long_length.txt")
					os.system("bjobs | awk '{print $4}' | grep large | wc -l > bjobs_large_length.txt")
					long_job_count = 0
					large_job_count = 0
					with open("bjobs_long_length.txt") as f:
						long_job_count = int(f.read().strip())
					with open("bjobs_large_length.txt") as f:
						large_job_count = int(f.read().strip())

					#queue assignment
					if long_job_count < 1600:
						queue = "long"
						break

					#changing from 900 to -1 so we don't try to use the large queue anymore
					elif large_job_count < -1:
						queue = "large"
						break

					#sleep for 1 second to not overburden the system
					os.system("sleep 1")

				#remove the length file to avoid clutter
				os.system("rm bjobs_long_length.txt bjobs_large_length.txt")				

				#prepare to run discovery on test params for this residue
				#start the command
				#cmd = "bsub -q long -W 168:00 -u \"\" -o output -e error -R \"rusage[mem=8000]\" \"python /pi/summer.thyme-umw/enamine-REAL-2.6billion/umass_chan_REAL-M_platform/rosetta/run_ligand_discovery_search.py "
				#removing the output and error std out
				cmd = "bsub -q " + queue + " -W 96:00 -u \"\" -R \"rusage[mem=10000]\" \"python /pi/summer.thyme-umw/enamine-REAL-2.6billion/umass_chan_REAL-M_platform/rosetta/run_ligand_discovery_search.py "
				#target pdb
				cmd = cmd + target_pdb + " "
				#anchor residue
				cmd = cmd + str(residue) + " "
				#motifs file
				cmd = cmd + motifs_file + " "
				#test_params directory
				cmd = cmd + tp_root + "/test_params/ "
				#atr cutoff, rep cutoff, ddg cutoff
				cmd = cmd + str(atr) + " " + str(rep) + " " + ddg + " "
				#extra args file (if any, and if not will append a blank string)
				cmd = cmd + extra_args_file
				#cap off
				cmd = cmd + "\""

				print(cmd)
				os.system(cmd)

				

