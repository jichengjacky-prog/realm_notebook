#the purpose of this script is to prepare a directory of up to 100 ligand conformers selected from shabedb as params in a test_params directory to work in Rosetta
#the script will be given a sub-list of up to 100 ligand conformers with accession data, and create a test_params directory that Rosetta liganddiscoverysearch can use for discovrey against a target
#the script will use supporting scripts to extract and decompress condensed params data from the 2.6B enamine library

#the script is intended to operate where it is called, and will look for a file named "conf_list.csv"
#the file will contain data formatted like this, thich has the shapedb score, ligand name and conformer, accession chunk, and subchunk
# ###
# -0.489485800266,PV-005633531035_4,29811,1
# -0.492682218552,Z4448503755_1,01035,7
# -0.492682218552,Z3789042056_8,01035,7
# -0.498856127262,PV-004964373284_3,38709,7
# -0.498856127262,PV-004964373283_11,38709,7
# -0.498872220516,PV-006084859944_3,31726,1
# -0.499529778957,PV-005674048644_3,25233,5
# -0.50066614151,Z4211293263_2,01366,2
###


#imports
import os,sys
import argparse
import subprocess
error_log = None

def run_cmd(cmd):
	global error_log
	result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
	if result.stdout:
		print(result.stdout, end="")
	if result.stderr:
		print(result.stderr, file=sys.stderr, end="")
		if error_log is not None:
			error_log.write("COMMAND: " + cmd + "\n")
			error_log.write(result.stderr)
	return result.returncode


def remove_error_line_from_confs_file(line_index, path):
	global error_log
	try:
		with open(path, "r") as f:
			lines = f.readlines()
		if 0 <= line_index < len(lines):
			del lines[line_index]
		with open(path, "w") as f:
			f.writelines(lines)
	except Exception as exc:
		message = f"Warning: failed to remove line {line_index} from {path}: {exc}\n"
		print(message, file=sys.stderr)
		if error_log is not None:
			error_log.write(message)


#no arguments needed since this is intended to operate in a preset directory
#clobber an existing test_params directory if there is one
run_cmd("rm -drf test_params")

#make initial test_params directory preparations
run_cmd("mkdir -p test_params")

#open the list file stream
orig_dir = os.getcwd()
confs_path = os.path.join(orig_dir, "conf_list.csv")
with open(confs_path, "r") as confs_file:
	confs_lines = confs_file.readlines()

#enter the directory
os.chdir("test_params")

#make necessary empty files that Rosetta needs to operate
run_cmd("touch exclude_pdb_component_list.txt patches.txt")

#open the error log for runtime command failures
error_log = open("error_log.txt", "w")

##open write stream to write teh residue types file
res_types_file = open("residue_types.txt", "w")

#write header
res_types_file.write("## the atom_type_set and mm-atom_type_set to be used for the subsequent parameter\n")
res_types_file.write("TYPE_SET_MODE full_atom\n")
res_types_file.write("ATOM_TYPE_SET fa_standard\n")
res_types_file.write("ELEMENT_SET default\n")
res_types_file.write("MM_ATOM_TYPE_SET fa_standard\n")
res_types_file.write("ORBITAL_TYPE_SET fa_standard\n")
res_types_file.write("## Params files\n")

#iterate over each ligand in the list file
#os.system("tar -xzf /pi/summer.thyme-umw/enamine-REAL-2.6billion/" + superchunk_str + "/" + working_chunk + "/condensed_params_and_db_" + str(i) + ".tar.gz condensed_params_and_db_" + str(i) + "/db.db -C .")
#os.system("tar -xzf /pi/summer.thyme-umw/enamine-REAL-2.6billion/0/00000/condensed_params_and_db_0.tar.gz single_conf_params/Z1020538478_shorthand_params.txt --strip-components=2 -C .")
for line_index, line in enumerate(confs_lines):
	#break up the line into components to work on accessing the conformer
	#shapedb score (useless at this point, it has served its purpose)
	shapedb_score = line.split(",")[0]
	#ligand name, separate from the conf number
	ligname = line.split(",")[1].split("_")[0]
	#ligand conformer number, separate from the ligand
	conf_num = line.split(",")[1].split("_")[1]
	#chunk and subchunk for library accession
	chunk = line.split(",")[2]
	subchunk = line.strip().split(",")[3]

	#derive teh superchunk, which is needed for library accession
	#derive the superchunk that this chunk belongs in for safer result storage (so we don't explode a directory with 53k directories)
	superchunk_str = chunk[0:3]

	#superchunk does not have preceeding zeroes, so cut any off, doing via casting
	superchunk_str = int(superchunk_str)
	superchunk_str = str(superchunk_str)

	#test print
	print(ligname + " " + conf_num)

	#extract the working file to the current location
	tar_cmd = (
		"tar -xzf /pi/summer.thyme-umw/enamine-REAL-2.6billion/" + superchunk_str + "/" + chunk + "/condensed_params_and_db_" + subchunk + ".tar.gz "
		+ "condensed_params_and_db_" + subchunk + "/single_conf_params/" + ligname + "_shorthand_params.txt --strip-components=2 -C ."
	)
	rc = run_cmd(tar_cmd)
	if rc != 0:
		print("Error: failed to extract shorthand params for " + ligname + "_" + conf_num + " from archive " + subchunk + " (chunk " + chunk + ").", file=sys.stderr)
		remove_error_line_from_confs_file(line_index, confs_path)
		continue

	#extract and clean the conformer params from the shorthand file
	extract_cmd = (
		"python " + os.path.abspath(os.path.join(os.path.dirname(__file__), "extract_single_param_from_condensed_file.py"))
		+ " " + ligname + "_shorthand_params.txt " + conf_num + " " + ligname + "_" + conf_num
	)
	rc = run_cmd(extract_cmd)
	if rc != 0:
		print("Error: extract_single_param_from_condensed_file failed for " + ligname + "_" + conf_num + ".", file=sys.stderr)
		# cleanup shorthand if present and skip this entry
		run_cmd("rm -f " + ligname + "_shorthand_params.txt")
		remove_error_line_from_confs_file(line_index, confs_path)
		continue

	#clean the spacing
	fix_cmd = (
		"python " + os.path.abspath(os.path.join(os.path.dirname(__file__), "fix_condensed_param_file_spacing.py"))
		+ " " + ligname + "_" + conf_num + ".params"
	)
	rc = run_cmd(fix_cmd)
	if rc != 0:
		print("Error: fix_condensed_param_file_spacing failed for " + ligname + "_" + conf_num + ".", file=sys.stderr)
		remove_error_line_from_confs_file(line_index, confs_path)
		continue

	#overwrite the fixed file over the bad spacing file
	mv_cmd = "mv fixed_" + ligname + "_" + conf_num + ".params " + ligname + "_" + conf_num + ".params"
	rc = run_cmd(mv_cmd)
	if rc != 0:
		print("Error: failed to move fixed params for " + ligname + "_" + conf_num + ".", file=sys.stderr)
		remove_error_line_from_confs_file(line_index, confs_path)
		continue

	#delete the shorthand params files from the working location

	run_cmd("rm -f " + ligname + "_shorthand_params.txt")

	#write the ligand and conformer params file to the residue_types file
	res_types_file.write(ligname + "_" + conf_num + ".params\n")


#close output files
if error_log is not None:
	error_log.close()
res_types_file.close()
