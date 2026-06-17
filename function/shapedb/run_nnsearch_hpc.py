#the purpose of this script is to run the shapedb NNSearch out of the shapedb container on a given chunk of the enamine library
#this will write shape similarity scores to output text files for processing to select te most similar ligands
#note, the target shape must be pre-aligned against the rest of the enamine library (ideally, aligned against suvo_shifted.sdf) or else you are almost guaranteed to get all 1 values, which is useless

#ideally, this script can be called by nnsearch_controller.py to automate running the script on the entire enamine library

#this will operate on a per-chunk basis, hitting all (up to) 10 subchunks within the chunk, so this can be compartmentalized into one chunk per job

#importantly, this script will also be a preliminary filter to remove any blacklisted ligands that are removed for having unwanted chemical fingerprints (or conformer shape, if we ever do that)

#imports
import os, sys, subprocess
from pathlib import Path


def run_cmd(cmd):
	result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
	if result.stdout:
		print(result.stdout, end="")
	if result.stderr:
		print(result.stderr, file=sys.stderr, end="")
	return result.returncode

def get_folder_size_pathlib(folder_path):
    root_directory = Path(folder_path)
    # rglob('*') finds all files recursively
    return sum(f.stat().st_size for f in root_directory.rglob('*') if f.is_file())


#get the chunk as a 5 digit code (i.e. the first chunk is 00000, not 0)
working_chunk = str(sys.argv[1])

#derive the superchunk that this chunk belongs in for safer result storage (so we don't explode a directory with 53k directories)
superchunk_str = working_chunk[0:3]

#superchunk does not have preceeding zeroes, so cut any off, doing via casting
superchunk_str = int(superchunk_str)
superchunk_str = str(superchunk_str)

#get the target shape molecule (Can have a path to it, must be mol2 or sdf, or shapedb likely won't be able to read it)
target_molecule_file = sys.argv[2]

target_molecule_name = target_molecule_file.split("/")[len(target_molecule_file.split("/")) - 1]

#create the result output in a specified location, or use location where this was called from otherwise
working_location = "."
if len(sys.argv) > 3:
	working_location = sys.argv[3]

realm_location = "/pi/summer.thyme-umw/Ji_rosetta_discovery"
if len(sys.argv) > 4:
	realm_location = sys.argv[4]

#move to the working location
os.chdir(working_location)

#make a directory to perform the analysis in, based on the superchunk name
run_cmd("mkdir -p " + superchunk_str)
os.chdir(superchunk_str)
#clobber any existing workign chunk
run_cmd("rm -drf " + working_chunk)
run_cmd("mkdir -p " + working_chunk)
os.chdir(working_chunk)

i=0
try:
	#iteratively decompress the params and db data to the working chunk location
	for i in range (0,10):


		run_cmd("tar -xzf /pi/summer.thyme-umw/enamine-REAL-2.6billion/" + superchunk_str + "/" + working_chunk + "/condensed_params_and_db_" + str(i) + ".tar.gz condensed_params_and_db_" + str(i) + "/db.db -C .")

		#run shapedb out of the container on the subchunk database, executed via singularity
		cmd = "singularity exec --bind condensed_params_and_db_" + str(i) + "/db.db:/input/db.db --bind " + target_molecule_file + ":/input/" + target_molecule_name + " " + realm_location + "/sif/shapedb_container.sif /pharmit/src/build/shapedb -NNSearch -k 100000 -ligand /input/" + target_molecule_name + " -db /input/db.db -print > " + working_chunk + "_" + str(i) + "_nn.txt"
		print(cmd)
		run_cmd(cmd)

		#run through the shapedb list and remove any ligands on the blacklist
		#TODO
		#create a write stream for the filered file
		filtered_file = open(working_chunk + "_" + str(i) + "_nn_filtered.txt", "w")

		#create read streams for the initial nn file and the blacklist
		initial_file = open(working_chunk + "_" + str(i) + "_nn.txt", "r")

		try:
			blacklist_file = open("/pi/summer.thyme-umw/enamine-REAL-2.6billion/" + superchunk_str + "/" + working_chunk + "/blacklist_file.csv" , "r")
		except:
			print("No blacklist file found for subchunk " + str(i) + " of chunk " + working_chunk + ", skipping blacklist filtering for this subchunk")
			blacklist_file = None

		blacklist_ligands = []
		if blacklist_file is not None:	
			#port the blacklist file ligands to a list
			for line in blacklist_file.readlines():
				if line.startswith("ligand,") == False:
					#to greatly improve speed and not waste time on competing subchunks, only take ligands from the subchunk
					#derive the subchunk, which is right before the ".sdf" anbd after the last underscore after that split
					subchunk = line.split(".sdf")[0].split("_")[len(line.split(".sdf")[0].split("_")) - 1]

					if str(i) == str(subchunk):
						blacklist_ligands.append(line.split(",")[0])

		#iterate over the initial file, and write lines to the filtered file if the ligand is not on the blacklist
		for line in initial_file.readlines():
			#derive the ligand in the initial file
			cur_lig = line.split("_")[0]

			#add if the ligand is not in the blacklist
			if cur_lig not in blacklist_ligands:
				filtered_file.write(line)

except:
	print("Something is wrong with subchunk " + str(i) + " of chunk " + working_chunk)
	#delete the unfiltered file

	
for i in range (0,10):
	run_cmd("rm " + working_chunk + "_" + str(i) + "_nn.txt")
		#delete the decompressed folder
	run_cmd("rm -drf condensed_params_and_db_" + str(i))




