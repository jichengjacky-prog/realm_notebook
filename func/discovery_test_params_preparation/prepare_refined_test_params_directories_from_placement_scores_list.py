#the purpose of this script is to prepare a directory of up to 250 conformers (using conformator; 250 is default) from a single ligand from its smiles string, which is obtained from a placement file

#this will create a test params directory (named test_params) for the ligand wherever this script is called

#imports
import os,sys

#the smiles string of the placed ligand
lig_smiles = sys.argv[1]

#the name of the ligand
lig_name = sys.argv[2]

license_key = ""

#optionally, include a license key string for conformator (which may be needed after the existing license for conformator in the container expires)
if len(sys.argv) == 4:
	license_key = sys.argv[3]

#no arguments needed since this is intended to operate in a preset directory
#clobber an existing test_params directory if there is one
os.system("rm -drf test_params")

#make initial test_params directory preparations
os.system("mkdir -p test_params")


#enter the directory
os.chdir("test_params")

#make necessary empty files that Rosetta needs to operate
os.system("touch exclude_pdb_component_list.txt patches.txt")

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
	os.system("singularity exec /pi/summer.thyme-umw/enamine-REAL-2.6billion/conformator_container.sif bash -lc \"/conformator_for_container/conformator_1.2.1/conformator --license \'" + license_key + "\' && /conformator_for_container/conformator_1.2.1/conformator -i " + lig_name + ".smi  -o " + lig_name + "_confs.sdf --keep3d --hydrogens -v 0\"")
else:
	os.system("singularity exec /pi/summer.thyme-umw/enamine-REAL-2.6billion/conformator_container.sif /conformator_for_container/conformator_1.2.1/conformator -i " + lig_name + ".smi  -o " + lig_name + "_confs.sdf --keep3d --hydrogens -v 0")

#use obabel to split the conformers file into individual conformer files
os.system("obabel -isdf " + lig_name + "_confs.sdf -O " + lig_name + "_.sdf -m")

#make a params file of each generated single params file
for r,d,f in os.walk(os.getcwd()):
	for file in f:
		#if it is a ligand conformer sdf
		if lig_name in file and file.endswith(".sdf") and file.endswith("_confs.sdf") == False:
			#run molfile to params
			#make a params file of the unique file
			os.system("singularity exec /pi/summer.thyme-umw/enamine-REAL-2.6billion/conformator_container.sif python /conformator_for_container/molfile_to_params.py " + file + " -n " + file.split(".sdf")[0] + " --keep-names --long-names --clobber --no-pdb")

			#add the params to the residue_types list
			res_types_file.write(file.split(".sdf")[0] + ".params\n")

#cleanup by deleting sdf and smi files
os.system("rm -drf *smi *sdf")
