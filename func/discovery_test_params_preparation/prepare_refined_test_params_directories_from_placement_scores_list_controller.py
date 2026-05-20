#the purpose of this script is to act as a wrapper controller to run the prepare_refined_test_params_directories_from_placement_scores_list.py script
#this script will take in an individual file of placement results
#for each placement/file (line) in the file, the smiles of the ligand will be extracted and used to make up to 250 conformers of the ligand with conformator to make params for another round of Rosetta discovery
#this will fire off a job per placement using the helper script to run the preparations. After this, the run rosetta script needs to be run to actually run discovery
#this script will operate in a given location
#it will make numbered folders that hold up to 100 directories for respective ligands

import os,sys
#from rdkit import Chem
#arguments

#placements summary file
summary_file = open(sys.argv[1],"r")

#working location
working_location = sys.argv[2]
os.chdir(working_location)

#optionally, take in the conformator license string (as an active license string will be mandatory after the build in the container expires in order for conformator to work)
license_key = ""
if len(sys.argv) == 4:
	license_key = sys.argv[3]

#declare counters and make the initial directory for the start of the setup
dir_counter = 0

#make a list that holds a list of ligand names that have been processed to avoid unecessary duplicates
ligand_name_list = []

os.system("mkdir " + str(dir_counter))

#/pi/summer.thyme-umw/rosetta_discovery_space/pth2/shapedb_results/top_hits/lower/lower_res1_11_shifted/0/0/122/placements/res122_7F16_receptor_only_PV-001704743841_6_0.pdb

#read over the summary file
for line in summary_file.readlines():
	#skip the header if the file has one
	if line.startswith("file,ddg,total_motifs,significant_motifs"):
		continue

	#otherwise, break down the pathing to the file to get the ligand name, fill file name to be extracted, and location of where the compressed placements file is
	full_file = line.split(",")[0]

	#path to compressed placements; will exist in the same place where the listed uncompressed placements location is
	compressed_placements_path = full_file.split("/placements/")[0]

	#placement file, which is after the last slash in the full name
	placement_file = full_file.split("/")[-1]

	#ligand name, which comes after the last / and is between underscores in the 3rd to last spot
	ligand_name = placement_file.split("_")[-3]

	#check if the ligand name is in the list
	#continue if it is, otherwise add the name
	if ligand_name in ligand_name_list:
		continue
	else:
		ligand_name_list.append(ligand_name)

	#make a folder in the numbered directory named after the ligand
	os.system("mkdir " + str(dir_counter) + "/" + ligand_name)

	#move into the location
	os.chdir(str(dir_counter) + "/" + ligand_name)

	#now that we know the file info, extract the placement
	os.system("tar -xzf " + compressed_placements_path + "/placements.tar.gz placements/" + placement_file + " --strip-components=1 -C .")

	#we have the file, extract the ligand
	os.system("grep HETATM " + placement_file + " > ligand.pdb")
	#os.system("echo TER >> ligand.pdb")
	#os.system("grep CONECT " + placement_file + " >> ligand.pdb")

	#extract the ligand smiles from the pdb file
	os.system("obabel -ipdb ligand.pdb -osmi -O ligand.smi -xn")

	#read the smiles file to get the smiles string
	smiles_file = open("ligand.smi", "r")

	smiles_string = ""
	for line in smiles_file.readlines():
		smiles_string = line.strip()

	#we now have the smiles string, clear the pdbs and smi files
	os.system("rm *pdb *smi")

	#now, run the prepare script in bsub
	#/pi/summer.thyme-umw/enamine-REAL-2.6billion/umass_chan_REAL-M_platform
	cmd = "bsub -q short -W 1:00 -u \"\" -R \"rusage[mem=4000]\" -o log.out -e log.err  \'python /pi/summer.thyme-umw/enamine-REAL-2.6billion/umass_chan_REAL-M_platform/discovery_test_params_preparation/prepare_refined_test_params_directories_from_placement_scores_list.py "
	#target pdb
	cmd = cmd + "\"" + smiles_string + "\" "
	#anchor residue
	cmd = cmd + str(ligand_name) + " "
	#motifs file
	cmd = cmd + "\"" + license_key + "\" "
	#cap off
	cmd = cmd + "\'"

	print(cmd)
	os.system(cmd)

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

	#at end, return to working location
	os.chdir(working_location)

	#tick up the directory counter every time the length of ligand name list mod 100 is 0
	if len(ligand_name_list) % 100 == 0:
		dir_counter = dir_counter + 1
		os.system("mkdir " + str(dir_counter))