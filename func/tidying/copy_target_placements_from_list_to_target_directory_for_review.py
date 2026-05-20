#the purpose of this script is to look through a placements list and take all placements in the list and put them in a new location for review
#this script assumes that the files are in compressed and dehydrated states, and will copy the full placements directory, attempt to rehydrate the directory, extract the file, and then delete the copy
#this will be parallelized across jobs


#imports
import os,sys
from subprocess import run
import shlex

#arguments
#the ligand list file
list_file = sys.argv[1]

#target_destination
target_destination = sys.argv[2]

#add end / if it does not exist for handling in this script
if target_destination.endswith("/") == False:
	target_destination = target_destination + "/"

#make the target destination
os.makedirs(target_destination, exist_ok=True)

#move into the location
os.chdir(target_destination)

#have a counter for the number of files and directories, we will try to restrict to 100 files per sub-directory to mitigate potential flooding
file_counter = 0
dir_counter = 0

os.makedirs(str(dir_counter), exist_ok=True)
#enter the directory
os.chdir(str(dir_counter))

#read through each line in the file
read_list_file = open(list_file,"r")

for line in read_list_file.readlines():
	#skip the header if there is one
	if line.startswith("file,ddg,"):
		continue
	else:
		#identify the file of interest and work on extracting it
		full_file = line.split(",")[0]
		#path to the tar placements file (placements should exist as compressed .tar.gz to save memory overhead)
		placements_file = line.split("/placement/")[0] + "/placements.tar.gz"
		#file root for folder naming for operations
		file_root = line.split("/placement/")[1].split(".pdb")[0]

		#increment the file counter and make a new subdir for every 100
		file_counter = file_counter + 1
		if file_counter % 100 == 0:
			dir_counter = dir_counter + 1
			os.chdir(target_destination)
			os.makedirs(str(dir_counter), exist_ok=True)
			os.chdir(str(dir_counter))

		#make a directory to work with for extracting the file
		os.makedirs(file_root, exist_ok=True)

		#enter the directory
		os.chdir(file_root)
		"""
		#make a bsub job to handle this
		#bsub parameters
		cmd = ["bsub","-q","short","-W","2:00","-u","\"\"","-R","\"rusage[mem=5000]\""]

		#start of command
		cmd.extend(["\""])

		#copy the placements tar
		cmd.extend(["cp",placements_file,"."])

		#couple
		cmd.extend(["&&"])

		#rehydrate the placements (this will result in still having a placements.tar.gz)
		cmd.extend(["python","/pi/summer.thyme-umw/enamine-REAL-2.6billion/umass_chan_REAL-M_platform/tidying/rehydrate_reduced_pdbs_with_skeleton_gpt.py"])

		#couple
		cmd.extend(["&&"])

		#extract the file of interest		
		cmd.extend(["tar","-xzf","--strip-components=1","placements.tar.gz","placements/" + file_root + ".pdb"])

		#couple
		cmd.extend(["&&"])

		#move the file up
		cmd.extend(["mv",file_root + ".pdb",".."])

		#couple
		cmd.extend(["&&"])

		#delete the placements
		cmd.extend(["rm","-drf","placements*"])

		#couple
		cmd.extend(["&&"])

		#move up a directory
		cmd.extend(["cd",".."])

		#couple
		cmd.extend(["&&"])

		#delete the deporary directory
		cmd.extend(["rm","-drf",file_root])

		#end of command
		cmd.extend(["\""])
		"""

		job_cmd = f"""
		cp {shlex.quote(placements_file)} . &&
		python /pi/summer.thyme-umw/enamine-REAL-2.6billion/umass_chan_REAL-M_platform/tidying/rehydrate_reduced_pdbs_with_skeleton_gpt.py &&
		tar -xzf placements.tar.gz --strip-components=1 {shlex.quote(f"placements/{file_root}.pdb")} &&
		mv {shlex.quote(file_root + ".pdb")} .. &&
		rm -rf placements* &&
		rm -drf ../{file_root}
		"""

		cmd = [
			"bsub",
			"-q", "short",
			"-W", "2:00",
			"-u", "",
			"-R", "rusage[mem=5000]",
			"bash", "-lc", job_cmd
		]

		#adding job throttle
		#bsub job throttle to make sure we do not exceed our local limit
		#write the length of the bjobs queue to this current location
		os.system("bjobs | grep short |  wc -l > bjobs_length.txt")
		job_count = 0
		with open("bjobs_length.txt") as f:
			job_count = int(f.read().strip())
		while job_count > 500:
			#sleep for 1 second to not overburden the system
			os.system("sleep 1")
			os.system("bjobs | grep short| wc -l > bjobs_length.txt")
			with open("bjobs_length.txt") as f:
				job_count = int(f.read().strip())
		#remove the length file to avoid clutter
		os.system("rm bjobs_length.txt")	

		print(cmd)
		#try:
		run(cmd, check=True)

		#move back up
		os.chdir("..")
