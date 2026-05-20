#the purpose of this script is to runthe dehydration script on all cases of there being the placements.tar.gz and raw_scores.csv in the same folder

import os,sys
from pathlib import Path

#arguments
#take in the root of the folder that you want to look down
cleaner_root = sys.argv[1]

#take in the reference pdb to make a skeleton and remove residues per pdb that did not move significantly
reference_pdb = sys.argv[2]

throttle_limit = 100

#optional argument for an override on the throttle limit
if len(sys.argv) == 4:
	throttle_limit = int(sys.argv[3])

#move to and iterate down the cleaner root
os.chdir(cleaner_root)

for r,d,f in os.walk(cleaner_root):
	#if there are cases where there is a placements.tar.gz and a raw_Scores.csv, then we know that we completed discovery and can try to dehydrate
	for file in f:
		if file == "placements.tar.gz":
			#heck if there is a raw_scores.csv
				if Path(r + "/raw_scores.csv").exists():
					#we have a placements and raw scores directory where we can try to ruthe dehydrator
					#run it in a lsf job and throttle
					print(r + "/raw_scores.csv " + r + "/placements.tar.gz")

					#try to fire off a lsf job to run the dehydrator
					cmd = "bsub -q short -W 1:00 -u \"\" -R \"rusage[mem=5000]\" \"python /pi/summer.thyme-umw/enamine-REAL-2.6billion/umass_chan_REAL-M_platform/tidying/shrink_placement_pdbs_to_placement_and_surrounding_residues.py "

					#reference pdb
					cmd = cmd + reference_pdb + " "

					#root location where the complete placements file is
					cmd = cmd + r + " "



					#cap off
					cmd = cmd + "\""

					print(cmd)
					os.system(cmd)

					#adding 10 job throttle
					#bsub job throttle to make sure we do not exceed our local limit
					#write the length of the bjobs queue to this current location
					os.system("bjobs | grep short |  wc -l > bjobs_length.txt")
					job_count = 0
					with open("bjobs_length.txt") as f:
						job_count = int(f.read().strip())
					while job_count > throttle_limit:
						#sleep for 1 second to not overburden the system
						os.system("sleep 1")
						os.system("bjobs | grep short| wc -l > bjobs_length.txt")
						with open("bjobs_length.txt") as f:
							job_count = int(f.read().strip())
					#remove the length file to avoid clutter
					os.system("rm bjobs_length.txt")			

