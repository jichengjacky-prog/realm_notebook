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
import os,sys

#get the working location and initial list file
working_location = sys.argv[1]
master_list = sys.argv[2]

#move to the working location
os.chdir(working_location)

#declare counters for the top level directories and sub directories
top_level_dirs = 0
sub_dirs = 0

#declare a working list to hold batches of up to 100 ligands to make the lists and jobs
small_confs_list = []

#create a starting directory for the first directory to be made
os.system("mkdir -p " + str(top_level_dirs))
os.system("mkdir -p " + str(top_level_dirs) + "/" + str(sub_dirs))

#move to the created location
os.chdir(str(top_level_dirs) + "/" + str(sub_dirs))

#create a file to write the small list to
write_file = open("conf_list.csv", "w")

#read through the master list, collect conformers, and make test_params directories every 100 encountered (or when at the end)
master_list_file = open(master_list,"r")
for line in master_list_file.readlines():
	#extract the ligand with conformer to ensure it is not in the working small_confs list
	ligconf = line.split(",")[1]
	#continue to avoid repeating, although this should generally be impossible. having the same ligand in the same call to rosetta twice could break things
	if ligconf in small_confs_list:
		print("WARNING: Encountered repeat of ligand conformer of " + ligconf)
		continue

	#assuming we can move forward with the ligand, add it to the small confs list and write it to the write file
	write_file.write(line)

	small_confs_list.append(ligconf)

	#if the size of teh small confs list reaches 100, prepare the run the sub script on this directory and then move on to the next
	if len(small_confs_list) == 100:
		#close the write stream
		write_file.close()

		#reset the small confs list
		small_confs_list = []

		#fire off a job to process the concluding batch
		#os.system("bsub -q long -W 8:00 -u \"\" \"python /pi/summer.thyme-umw/enamine-REAL-2.6billion/umass_chan_REAL-M_platform/shapedb/run_nnsearch_hpc.py " + chunk_str + " " +  target_molecule_file + " " + working_location + " \"")
		os.system("bsub -q long -W 8:00 -u \"\" \"python /pi/summer.thyme-umw/enamine-REAL-2.6billion/umass_chan_REAL-M_platform/discovery_test_params_preparation/prepare_test_params_directories_sub.py\"")

		#throttle progress if there are too many job in queue
		#bsub job throttle to make sure we do not exceed our local limit
		#write the length of the bjobs queue to this current location
		os.system("bjobs | wc -l > bjobs_length.txt")
		job_count = 0
		with open("bjobs_length.txt") as f:
			job_count = int(f.read().strip())
		while job_count > 1500:
			#sleep for 1 second to not overburden the system
			os.system("sleep 1")
			os.system("bjobs | wc -l > bjobs_length.txt")
			with open("bjobs_length.txt") as f:
				job_count = int(f.read().strip())
		#remove the length file to avoid clutter
		os.system("rm bjobs_length.txt")		

		#move back to the top and make a new directory set and write file
		#move to the working location
		os.chdir(working_location)

		#increment the lower directory and work from there
		sub_dirs = sub_dirs + 1

		#behavior for if we make 100 subdirs
		if sub_dirs == 100:
			sub_dirs = 0
			top_level_dirs = top_level_dirs + 1

		#reset the small list, make directories, and make the next small list file
		#create a starting directory for the first directory to be made
		os.system("mkdir -p " + str(top_level_dirs))
		os.system("mkdir -p " + str(top_level_dirs) + "/" + str(sub_dirs))

		#move to the created location
		os.chdir(str(top_level_dirs) + "/" + str(sub_dirs))

		#create a file to write the small list to
		write_file = open("conf_list.csv", "w")		

#end behavior for final list
#close the final file
write_file.close()

#submit the final job
os.system("bsub -q long -W 8:00 -u \"\" \"python /pi/summer.thyme-umw/enamine-REAL-2.6billion/umass_chan_REAL-M_platform/discovery_test_params_preparation/prepare_test_params_directories_sub.py\"")