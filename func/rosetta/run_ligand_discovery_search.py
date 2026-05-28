#The purpose of this script is to take in specified inputs to run Rosetta ligand discovery search
#This is meant to run a single call of liganddiscoverysearch and be called by a controller script that controls calling larger batches of the discovery search
#mandatory inputs include:
#1. target pdb
#2. an anchor residue  string of all positions to be used as anchors, which is comma separated if there is more than 1 (i.e. 79 for 1 or 11,79,55,403 for multiple); this residue needs to use rosetta indexing, which needs to be determined before running the pipeline; this program itself if called by the controller will only use a single anchor residue in the call
#3. a motifs file (by default, will use the main 1.6M motif library, but if another is specified, will be used instead)
#4. location to a test_params folder containing ligands of interest to be docked

#additional arguments that we will treat as mandatory (5-7)
#score cutoff overrides for fa_atr, fa_rep, and ddg

#an option to bring in a text file containing arguments to be appended to the discovery arguments for things like a space fill region or to print the space fill matrix (not recommended forr memory and speed)
#dimensions and center for a space fill matrix (which can speed up discovery by filtering away placements that are not close enough to the desired binding region)

#Rosetta will be called to run in the location where this script is called

#example command
#python /pi/summer.thyme-umw/enamine-REAL-2.6billion/umass_chan_REAL-M_platform/rosetta/run_ligand_discovery_search.py /pi/summer.thyme-umw/rosetta_discovery_space/pth2/thymelab_pth2_discovery/pth2_structures/7F16_receptor_only.pdb 63,87,96,179 /pi/summer.thyme-umw/enamine-REAL-2.6billion/FINAL_motifs_list_filtered_2_3_2023.motifs /pi/summer.thyme-umw/rosetta_discovery_space/pth2/shapedb_results/top_hits/upper/upper_res29_34_shifted/test/test_params/ -2 150 -9 /pi/summer.thyme-umw/rosetta_discovery_space/pth2/shapedb_results/top_hits/upper/upper_res29_34_shifted/test/extra_args
#which calls the container like
#singularity exec --bind test_params:/input/test_params --bind test_args:/input/test_args --bind /pi/summer.thyme-umw/rosetta_discovery_space/pth2/thymelab_pth2_discovery/pth2_structures/7F16_receptor_only.pdb:/input/7F16_receptor_only.pdb --bind /pi/summer.thyme-umw/2024_intern_lab_space/FINAL_motifs_list_filtered_2_3_2023.motifs:/input/FINAL_motifs_list_filtered_2_3_2023.motifs /pi/summer.thyme-umw/2024_intern_lab_space/ari_work/containers/rosetta_condensed_6_25_2024.sif /rosetta/source/bin/ligand_discovery_search_protocol.linuxgccrelease @/input/test_args

#imports
import argparse
import os
import subprocess
import sys
from pathlib import Path

def run_cmd(cmd):
	result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
	if result.stdout:
		print(result.stdout, end="")
	if result.stderr:
		print(result.stderr, file=sys.stderr, end="")
	return result.returncode

# parse CLI arguments
parser = argparse.ArgumentParser(
	description="Run Rosetta ligand discovery search for a single test_params directory."
)
parser.add_argument("target_pdb", help="Target PDB file path")
parser.add_argument(
	"anchor_residue_string",
	help="Comma-separated Rosetta-indexed anchor residue(s), e.g. 79 or 11,79,55,403"
)
parser.add_argument("motifs_file", help="Motifs file path")
parser.add_argument(
	"test_params_dir",
	help="Path to the test_params directory",
)
parser.add_argument(
	"discovery_directory_root",
	help="Root directory to search for test_params directories"
)

parser.add_argument("atr", help="fa_atr cutoff")
parser.add_argument("rep", help="fa_rep cutoff")
parser.add_argument("ddg", help="ddg cutoff")
parser.add_argument(
	"--extra-args-file",
	dest="extra_args_file",
	default="",
	help="Optional file with additional Rosetta discovery arguments"
)
args = parser.parse_args()

# target pdb
target_pdb = args.target_pdb
# anchor residue(s)
anchor_residue_string = args.anchor_residue_string
# motifs file (likely want /pi/summer.thyme-umw/enamine-REAL-2.6billion/FINAL_motifs_list_filtered_2_3_2023.motifs unless you know what you are doing)
motifs_file = args.motifs_file
# test_params directory (needs to end with a /, also needs to be explicitly named "test_params" due to how functions in Rosetta work)
test_params_dir = args.test_params_dir
discovery_directory_root = args.discovery_directory_root
if test_params_dir.endswith("/") == False:
	test_params_dir = test_params_dir + "/"
# atr, rep, ddg cutoffs
atr = args.atr
rep = args.rep
ddg = args.ddg

# if there is an extra args file, take it in
extra_args_file = args.extra_args_file

#for the input files/directories, they need to be mapped to a container location for the Rosetta container (by default, will map to the /input location for the singularity call)
#here we will make strings for the mapping, which will go in the rosetta args file (the actual paths will be used in executing in the singularity image)
#input_target_pdb = "/input/" + target_pdb.split("/")[len(target_pdb.split("/")) - 1]
#input_motifs_file = "/input/" + motifs_file.split("/")[len(motifs_file.split("/")) - 1]
#input_test_params_dir = "/input/test_params/"
#at least for now, if there are other input files that are added in the extra args, they can't be supported unless I improve the mapping logic in the call to rosetta executaion in the container (since files would have to be recognized and mapped)

#now that we have all the args, compose an args file for discovery in the location where this is called
args_file = open("args","w")

#these first few args are hard coded for housekeeping, but can be removed later if we really want these to be mutable/not included
#if someone really wants to change these, they can include the arg again in the extra args file, and the arg will get overwritten with the desired value
args_file.write("#keep seed constant\n")
args_file.write("-constant_seed 1\n")
args_file.write("#ignore unrecognized residues to help mitigate crashes\n")
args_file.write("-ignore_unrecognized_res\n")
args_file.write("#handle ligand repeats if using multiple anchor residues, will otherwise crash without this flag\n")
args_file.write("-in::file::override_database_params true\n")
args_file.write("#constrain coordinates\n")
args_file.write("-constrain_relax_to_start_coords\n")
args_file.write("#keep all placements; 0 means keep all, any other integer means keep up to that integer\n")
args_file.write("-best_pdbs_to_keep 0\n")

#user input dependent
args_file.write("#mapped protein system\n")
args_file.write("-s " + "/input/" + str(Path(target_pdb).name) + "\n")
args_file.write("#mapped motifs file\n")
args_file.write("-motif_filename " + "/input/" + str(Path(motifs_file).name) + "\n")
args_file.write("#mapped test_params directory\n")
args_file.write("-params_directory_path " + "/input/" + str(Path(test_params_dir).name) +"/" "\n")
args_file.write("#rosetta-indexed anchor residue index/indices\n")
args_file.write("-protein_discovery_locus " + anchor_residue_string + "\n")
args_file.write("#fa_atr cutoff\n")
args_file.write("-fa_atr_cutoff = " + atr + "\n")
args_file.write("#fa_rep cutoff\n")
args_file.write("-fa_rep_cutoff = " + rep + "\n")
args_file.write("#ddg cutoff\n")
args_file.write("-ddg_cutoff = " + ddg + "\n")

#if the user wanted to add extra args, add them from the extra args file
if extra_args_file != "":
	args_file.write("###################################################\n")
	args_file.write("#extra user args from: " + extra_args_file + "\n")

	#open the file, read the lines and write the the working args file
	read_file = open(extra_args_file,"r")
	for line in read_file.readlines():
		args_file.write(line)

args_file.close()

print(test_params_dir)

rosetta_cmd =[ "singularity exec " ,
		"--bind " + test_params_dir +":" + "/input/test_params/" ,
		" --bind " + os.getcwd() + "/args:/input/args" ,
		" --bind " + target_pdb + ":" + "/input/" + str(Path(target_pdb).name) ,
		" --bind " + motifs_file + ":" + "/input/" + str(Path(motifs_file).name),
		discovery_directory_root + "/sif/rosetta_condensed_6_25_2024.sif",
		"/rosetta/source/bin/ligand_discovery_search_protocol.linuxgccrelease @/input/args"]


#we now have the args file written, now call Rosetta discovery
run_cmd(' '.join(rosetta_cmd))
#move all pdb files to a placements directory
run_cmd("mkdir placements")

run_cmd("mv *pdb placements")

os.chdir("placements")

#rename each pdb file by prepending the anchor residue string used by this script
for r,d,f in os.walk(os.getcwd()):
	for file in f:
		if file.endswith(".pdb") and r == os.getcwd():
			os.system("mv " + file + " res" + anchor_residue_string + "_" + file)

#now, call the placement analysis script
run_cmd("python " + discovery_directory_root + "/func/rosetta/score_placed_ligands_with_filtering.py")

#copy the csv files up a level for easy accession outside of the to-be compressed placements directory
run_cmd("cp *csv ..")

#then, compress the placement files (this will move all pdb files to a directory called placements, so do not keep any important pdbs in here)
os.chdir("..")

run_cmd("tar -czf placements.tar.gz placements")

run_cmd("rm -drf placements")

#run the dehydrate script to minimize overhead (until it is time to process the discovery results)
run_cmd("python " + discovery_directory_root + "/func/tidying/shrink_placement_pdbs_to_placement_and_surrounding_residues.py " + target_pdb)