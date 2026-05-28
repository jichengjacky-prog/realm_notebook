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
import argparse
import os
import subprocess
import sys


def run_cmd(cmd):
	result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
	if result.stdout:
		print(result.stdout, end="")
	if result.stderr:
		print(result.stderr, file=sys.stderr, end="")
	return result.returncode

# parse CLI arguments
parser = argparse.ArgumentParser(
	description="Run Rosetta ligand discovery search across a tree of test_params directories."
)
parser.add_argument("target_pdb", help="Target PDB file path")
parser.add_argument(
	"anchor_residue_string",
	help="Comma-separated anchor residue indices (Rosetta indexing), e.g. 79 or 11,79,55,403"
)
parser.add_argument("motifs_file", help="Motifs file path")
parser.add_argument(
	"discovery_directory_root",
	help="Root directory to search for test_params directories"
)
parser.add_argument(
	"shapedb_output_dir",
	help="Directory to store shapedb output"
)
parser.add_argument("atr", help="fa_atr cutoff")
parser.add_argument("rep", help="fa_rep cutoff")
parser.add_argument("ddg", help="ddg cutoff")
def str2bool(v):
	if isinstance(v, bool):
		return v
	if v.lower() in ('yes', 'true', 't', '1'):
		return True
	elif v.lower() in ('no','false', 'f', '0'):
		return False
	else:
		raise argparse.ArgumentTypeError('Boolean value expected.')

parser.add_argument(
	"clobber",
	type=str2bool,
	help="True/False or 1/0 to clobber existing completed discovery directories"
)
parser.add_argument(
	"--extra-args-file",
	dest="extra_args_file",
	default="",
	help="Optional file containing additional Rosetta discovery arguments"
)
args = parser.parse_args()

# target pdb
target_pdb = args.target_pdb
# anchor residue(s)
anchor_residue_string = args.anchor_residue_string
# break apart the string by commas into a list
anchor_residue_string_list = anchor_residue_string.split(",")
# motifs file
motifs_file = args.motifs_file
# discovery directory root
discovery_directory_root = args.discovery_directory_root
# shapedb output directory
shapedb_output_dir = args.shapedb_output_dir
# atr, rep, ddg cutoffs
atr = args.atr
rep = args.rep
ddg = args.ddg
# clobber existing output directories
clobber = args.clobber
# optional extra args file
extra_args_file = args.extra_args_file

#look over the discovery directory root to identify test_params directories
for r,d,f in os.walk(shapedb_output_dir):
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
				if clobber:
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
				run_cmd("rm -drf " + str(residue))

				#make a directory for the anchor residue
				run_cmd("mkdir -p " + str(residue))

				#enter the directory
				os.chdir(str(residue))

				#determine whether to send the next job to the long or large queue, based on the number of jobs running in each queue
				#if there are under 1600 long jobs, send to long
				#else if there are under 900 large jobs, send to large
				#otherwise if both buffered queues are running and full, throttle until space opens up
				os.system("sleep 0.5")

				#prepare to run discovery on test params for this residue
				#start the command
				#removing the output and error std out
				cmd = ["bsub -q long -W 96:0 -n 1 -U \"\" -R \"rusage[mem=10000]\"", 
		   			"python", discovery_directory_root +"/func/rosetta/run_ligand_discovery_search.py",
				target_pdb,
				str(residue),
			   	motifs_file,
			   	tp_root + "/test_params/",
				discovery_directory_root,
			  	atr,
			   	rep,
			   	ddg,
			   	extra_args_file]
				

				print(' '.join(cmd))
				run_cmd(' '.join(cmd))

				

