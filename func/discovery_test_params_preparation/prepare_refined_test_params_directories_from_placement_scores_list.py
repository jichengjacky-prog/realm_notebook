#the purpose of this script is to prepare a directory of up to 250 conformers (using conformator; 250 is default) from a single ligand from its smiles string, which is obtained from a placement file

#this will create a test params directory (named test_params) for the ligand wherever this script is called

#imports
import os
import sys
import argparse
import subprocess


def parse_args():
	parser = argparse.ArgumentParser(
		description='Prepare conformer directory from SMILES string using Conformator'
	)
	parser.add_argument(
		'smiles',
		help='SMILES string of the ligand'
	)
	parser.add_argument(
		'ligand_name',
		help='Name of the ligand (used for file naming)'
	)
	parser.add_argument(
		'--license-key',
		default='',
		help='Optional Conformator license key (if needed for container)'
	)
	parser.add_argument(
	    'realm_location',
		default='None',
		help='Optional Conformator license key (if needed for container)'
	)
	return parser.parse_args()


def run_cmd(cmd):
	result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
	if result.stdout:
		print(result.stdout, end="")
	if result.stderr:
		print(result.stderr, file=sys.stderr, end="")
	return result.returncode


args = parse_args()
lig_smiles = args.smiles
lig_name = args.ligand_name
license_key = args.license_key
if args.realm_location != "None":
	realm_location = args.realm_location
else:
	try:
		realm_location = os.getcwd()
	except Exception as e:
		print(f"Error getting current working directory: {e}")
		sys.exit(1)	
	

#no arguments needed since this is intended to operate in a preset directory
#clobber an existing test_params directory if there is one
run_cmd("rm -drf test_params")

#make initial test_params directory preparations
run_cmd("mkdir -p test_params")


#enter the directory
os.chdir("test_params")

#make necessary empty files that Rosetta needs to operate
run_cmd("touch exclude_pdb_component_list.txt patches.txt")

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

#make a smiles file from the smiles string and ligand name
#COc1cc(C)cnc1C(=O)NC[C@H](C)N(C)C(=O)c1ccc2ccccc2n1
smiles_file = open(lig_name + ".smi", "w")
smiles_file.write(lig_smiles)
smiles_file.close()

#run conformator out of the conformator container
#determine whether to run command with or without using license activation
if license_key != "":
	cmd = ["singularity exec", "sif/conformator_container.sif", 
		"bash", "-lc", 
		"\"/conformator_for_container/conformator_1.2.1/conformator --license \'" + license_key + "\' && /conformator_for_container/conformator_1.2.1/conformator -i " + lig_name + ".smi  -o " + lig_name + "_confs.sdf --keep3d --hydrogens -v 0\""]
	run_cmd(" ".join(cmd))
else:
	cmd = ["singularity exec", "sif/conformator_container.sif", 
		"bash", "-lc", 
		"\"/conformator_for_container/conformator_1.2.1/conformator -i " + lig_name + ".smi  -o " + lig_name + "_confs.sdf --keep3d --hydrogens -v 0\""]
	run_cmd(" ".join(cmd))
	
#use obabel to split the conformers file into individual conformer files
run_cmd("obabel -isdf " + lig_name + "_confs.sdf -O " + lig_name + "_.sdf -m")

#make a params file of each generated single params file
for r,d,f in os.walk(os.getcwd()):
	for file in f:
		#if it is a ligand conformer sdf
		if lig_name in file and file.endswith(".sdf") and file.endswith("_confs.sdf") == False:
			#run molfile to params
			#make a params file of the unique file
			cmd = ["singularity exec", "sif/conformator_container.sif", "python", "/conformator_for_container/molfile_to_params.py", file, "-n", file.split(".sdf")[0], "--keep-names", "--long-names", "--clobber", "--no-pdb"]
			run_cmd(" ".join(cmd))

			#add the params to the residue_types list
			res_types_file.write(file.split(".sdf")[0] + ".params\n")

#cleanup by deleting sdf and smi files
run_cmd("rm -drf *smi *sdf")
