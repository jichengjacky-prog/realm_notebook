#script to look at a directory of pdbs of placed ligands using Rosetta's liganddiscoverysearch
#ligands are scored based on values derived from rosetta (ddg, motif counts and real ratio), ligand strain energy with tldr, and attepmter recovery with AutoDock vina
#system scores are summative of all scoring values used to calculate a score, with a larger and more positive value ideally being more favorable

#package imports
import os,sys
import argparse

# Create the parser
parser = argparse.ArgumentParser(
	description="This script looks at a directory of placed ligand pdb files from Rosetta's liganddiscovery search. It scores ligands based on values from Rosetta (ddg, motif counts and real ratio), ligand strain energy using tldr STRAIN, and attempted recovery of the Rosetta prediction with AutoDock Vina."
)

# Add the working location path argument
#optional; if not used, will use current location to look for pdb files
parser.add_argument(
	'-w', '--working_location',
	type=str,
	required=False,
	help='(Optional) Full path to the working location where the pdb files are. Will use the current location if not specified.'
)

# Add the autodock automation script path argument
#run_autodock_on_placed_ligands.py
#optional, will not attempt autodock if flag is not used; also requires the use of the -v flag that indicates where the build of autodock vina is
parser.add_argument(
	'-a', '--autodock_script_path',
	type=str,
	required=False,
	help='(Optional) Full path to the location where the run_autodock_on_placed_ligands.py script is located for attempting to recover Rosetta placements with AutoDock Vina. Recovery will not be attempted if this flag is not used or if the path to the AutoDock Vina Executable (-v, --autodock_vina_path) is not used. This is designed and tested on Linux using AutoDock Vina version 1.1.2, and may not work on other platforms.'
)

# Add the AutoDock Vina executable path argument
#looking for executable called "vina"
#optional, will not attempt autodock if flag is not used; also requires the use of the -a flag that indicates where the autodock automation script is
parser.add_argument(
	'-v', '--autodock_vina_path',
	type=str,
	required=False,
	help='(Optional) Full path to the location where the AutoDock Vina "vina" executable is located for attempting to recover Rosetta placements with AutoDock Vina. Recovery will not be attempted if this flag is not used or if the path to the AutoDock Vina automation script run_autodock_on_placed_ligands.py (-a, --autodock_script_path) is not used. This is designed and tested on Linux using AutoDock Vina version 1.1.2, and may not work on other platforms or versions.'
)

# Add the tldr STRAIN executable path argument
#looking for a script called run_torsion_check_on_placed_ligands.py
#optional, will not attempt strain if  if flag is not used; also requires the use of the -t flag that indicates where the autodock automation script is
parser.add_argument(
	'-r', '--torsion_script_path',
	type=str,
	required=False,
	help='(Optional) Full path to the location where the run_torsion_check_on_placed_ligands.py script is located for evaluating conformer strain energies. Will not be attempted if this flag is not used or if the path to the tldr strain automation script run_autodock_on_placed_ligands.py (-t, --torsion_script_path) is not used. This is designed and tested on Linux, and may not work on other platforms. The script requires a python or conda environment to have openbabel (obabel) to be installed.'
)

# Add the tldr STRAIN executable path argument
#looking for a script called Torsion_Strain.py
#optional, will not attempt strain if  if flag is not used; also requires the use of the -r flag that indicates where the autodock automation script is
parser.add_argument(
	'-t', '--torsion_strain_path',
	type=str,
	required=False,
	help='(Optional) Full path to the location where the tldr strain Torsion_Strain.py script is located for evaluating conformer strain energies. Will not be attempted if this flag is not used or if the path to the tldr strain automation script run_autodock_on_placed_ligands.py (-r, --torsion_strain_path) is not used. This is designed and tested on Linux, and may not work on other platforms. The script requires a python or conda environment to have rdkit to be installed.'
)

# Add the score_weights file path argument
#looking for a file called score_weights.csv
#optional, collects user-set score weights for values used on scoring placed ligands; any value that does not have a weight set is defaulted to 1
parser.add_argument(
	'-s', '--weights_path',
	type=str,
	required=False,
	help='(Optional) Full path to the location where the score weights file score_weights.csv is located. The format within the file is that each line has a score term at index 0 and the score weight at index 1.'
)

# Add the minimum_real_motif_ratio argument
#optional, collects user-set score weights for values used on scoring placed ligands; any value that does not have a weight set is defaulted to 1
parser.add_argument(
	'-m', '--minimum_real_motif_ratio',
	type=float,
	required=False,
	help='(Optional) The minimum ratio of real motifs for a placement to be considered. Any placements with a lower value will not be listed or have autodock/STRAIN be used on them.'
)

# Add the kill argument
#optional, program will not ignore placements if Autodock fails to place the ligand or if STRAIN can not derive an energy score
parser.add_argument(
	'-k', '--kill',
	type=bool,
	required=False,
	help='(Optional) Program will ignore placements if Autodock fails to place the ligand or if STRAIN can not derive an energy score.'
)

#from this option and below, will not use short arguments (i.e. -k), just because there are so many options and only so many letters in the english alphabet
# Add the motif hbond option
#optional, program will read the motifs generated by the placement, and collect the hbond energies to print out
parser.add_argument(
	'--get_motif_hbond_energies',
	type=bool,
	required=False,
	help='(Optional) Program will collect all hbond energies and will derive a count of motifs with nonzero hbond energies.'
)

# Add the motif index translation key option
#optional, program will read a file (can include path) to a csv key file that allows for translation of residue indeces to a desired index
#this is useful because Rosetta can potentially mess up the index of residues between the input and output, and a key file can correct for it.
#line structure is: res_type,original_index,translated_index,difference_in_index (optional)
parser.add_argument(
	'--residue_correction_key_file',
	type=str,
	required=False,
	help='(Optional) Program will read a file (can include path) to a csv key file that allows for translation of residue indeces to a desired index.'
)

# Add option for 
#optional, program requires that a motif (real or not) be found against desired residue indices found in option. Option is comma separated list (with no space) of residues to select for; trailing comma not necessary
#example inputs 86 86,90,112
parser.add_argument(
	'--mandatory_motif_residues',
	type=str,
	required=False,
	help='(Optional) Residues that must have motifs (real or not) identified against them, or else the placed ligand is ignored after checking.'
)

# Add the maximum_autodock_recovery_rmsd
#optional, maximum allowed recovery of autodock rmsd for a placement to be considered. A higher lowest rmsd value and the placement will not be considered
parser.add_argument(
	'--maximum_autodock_recovery_rmsd',
	type=float,
	required=False,
	help='(Optional) Maximum closest rmsd recovery value from Autodock to recover the placed ligand. Any higher closest value and the placement will be ignored.'
)

# Add the maximum_strain_energy
#optional, maximum allowed strain energy on a ligand. A higher strain energy value and the placement will not be considered
parser.add_argument(
	'--maximum_strain_energy',
	type=float,
	required=False,
	help='(Optional) Maximum strain energy for a placed ligand. Any higher strain energy and the placement will be ignored.'
)

#parse arguments
args = parser.parse_args()

#debug print of args
print(args)

#initial setup of working path
working_path = ""

#set up location variable to be either working path or current location,
location = os.getcwd()

#set up working path if -w or --wpath was used
#if 'w' in vars(args) or 'wpath' in vars(args):
if args.working_location != None:
	working_path = args.working_location

	#remove backslash from end of working path if it is there
	if working_path.endswith("/") == True:
		working_path = working_path[:-1]

	location = working_path

	#move to the location
	os.chdir(location)

#determine whether to use autodock and set up autodock variables
autodock_exec_path = ""
autodock_automate_path = ""
#default true, set to false if either relevant flag is not used
attempt_autodock = True
#autodock automate path
#if 'v' in vars(args) or 'autodock_vina_path' in vars(args):
if args.autodock_vina_path != None:
	autodock_exec_path = args.autodock_vina_path
else:
	attempt_autodock = False

#autodock automate path
#if 'a' in vars(args) or 'autodock_script_path' in vars(args):
if args.autodock_script_path != None:
	autodock_automate_path = args.autodock_script_path
	#print(autodock_automate_path)
	#add backslash to end of path if there is not one
	if autodock_automate_path.endswith("/") == False:
		autodock_automate_path = autodock_automate_path + "/"
	#print(autodock_automate_path)
else:
	attempt_autodock = False

#print('-v' in vars(args), '--autodock_vina_path' in vars(args), '-v' in vars(args) or '--autodock_vina_path' in vars(args), '-a' in vars(args), '--autodock_script_path' in vars(args), '-a' in vars(args) or '--autodock_script_path' in vars(args))
#print('v' in vars(args), 'autodock_vina_path' in vars(args), 'v' in vars(args) or 'autodock_vina_path' in vars(args), 'a' in vars(args), 'autodock_script_path' in vars(args), 'a' in vars(args) or 'autodock_script_path' in vars(args))

#determine whether to use tldr strain and set up variables
strain_exec_path = ""
strain_automate_path = ""
#default true, set to false if either relevant flag is not used
attempt_strain = True
#strain torsion path
#if 't' in vars(args) or 'torsion_script_path' in vars(args):
if args.torsion_script_path != None:
	strain_exec_path = args.torsion_script_path

	#add backslash to end of path if there is not one
	if strain_exec_path.endswith("/") == False:
		strain_exec_path = strain_exec_path + "/"
else:
	attempt_strain = False

#strain automate path
#if 'r' in vars(args) or 'torsion_strain_path' in vars(args):
if args.torsion_strain_path != None:
	strain_automate_path = args.torsion_strain_path

	#add backslash to end of path if there is not one
	if strain_automate_path.endswith("/") == False:
		strain_automate_path = strain_automate_path + "/"
else:
	attempt_strain = False

#prepare score weights, set up initial score weight dictionary
#if adjusting weights in score_weights.csv, any changed weight must match one of these entries exactly (including case)
#if a value's weight is set to zero and would otherwise be used, it will not be considered in scoring
score_weights = {"ddg":1,"total_motifs":1,"significant_motifs":1,"real_motif_ratio":1,"hbond_motif_count":1,"hbond_motif_energy_sum":1,"closest_autodock_recovery_rmsd":1,"closest_autodock_recovery_ddg":1,"strain_energy":1}

#if weights_path flag is set, work on adjusting weights
#if 's' in vars(args) or 'weights_path' in vars(args):
if args.weights_path != None:
	weights_path = args.weights_path

	#add backslash to end of path if there is not one
	if weights_path.endswith("/") == False:
		weights_path = weights_path + "/"

	#open the weights file and read the weights data
	weights_file = open(weights_path + "score_weights.csv", "r")

	#read each line and get the score term and weight
	#when splitting line by "," the score term should be index 0 and score be index 1
	#not all score terms need to be in the csv file, and order of terms does not matter
	#there is no header line
	#i.e.:
	"""
	ddg,3
	real_motif_ratio,10
	closest_autodock_recovery_ddg,0
	strain_energy,0.5
	"""
	for line in weights_file.readlines():
		term = line.split(",")[0]
		weight = float(line.split(",")[1].strip())
		
		#adjust the term
		if term in score_weights.keys():
			score_weights[term] = weight
		else:
			print("WARNING: included term '" + term + "' is not valid!")
			print("Valid terms are: " + score_weights.keys())

#minimum motif ratio
minimum_real_motif_ratio = 0
#set ratio if flag is used
if args.minimum_real_motif_ratio != None:
	minimum_real_motif_ratio = args.minimum_real_motif_ratio

#kill failed placements for strain or autodock
kill = False
if args.kill != None:
	kill = args.kill

#set up the key translation file if options provided for it
#mostly concerned with the original and translated indices, but the other 2 can be helpful when translating 
residue_index_dict = {}
use_key = False
if args.residue_correction_key_file != None:
	use_key = True
		#get the key file and read it
	key_file = open(args.residue_correction_key_file,"r")

	for line in key_file.readlines():
				#skip the first line, starts with res_type
		if line.startswith("res_type"):
			continue

				#seed the dictionary with the original and translated indices
		residue_index_dict[line.split(",")[1]] = line.split(",")[2]

#get motif hbond energies
get_hbond_energies = False
if args.get_motif_hbond_energies != None:
	get_hbond_energies = True

#get the mandatory motif residues
#have an empty list, and fille the list if the option was used. the program won't do anything with the lsit if it is empty
mandatory_motifs_list = []
if args.mandatory_motif_residues != None:
	mandatory_motifs_list = args.mandatory_motif_residues.split(",")

#get the maximum closest autodock recovery rmsd
#set as high value of 100
autodock_rmsd_cutoff = 100
if args.maximum_autodock_recovery_rmsd != None:
	autodock_rmsd_cutoff = args.maximum_autodock_recovery_rmsd

#get the maximum strain torsion rmsd
#set as high value of 100
strain_cutoff = 100
if args.maximum_strain_energy != None:
	strain_cutoff = args.maximum_strain_energy

#set up a file to hold the scoring for the placed ligands in csv format both with weighted and non-weighted values
#write the file in the location
#raw score file, generally for debugging and score term weights are not applied
raw_score_file = open(location + "/raw_scores.csv", "w")
#write header line
raw_score_file.write("file,ddg,total_motifs,significant_motifs,real_motif_ratio,hbond_motif_count,hbond_motif_energy_sum,closest_autodock_recovery_rmsd,closest_autodock_recovery_ddg,strain_energy,total\n")
#set up weighted score file; this file should be used for selection
weighted_score_file = open(location + "/weighted_scores.csv", "w")
#write header line
weighted_score_file.write("file,ddg,total_motifs,significant_motifs,real_motif_ratio,hbond_motif_count,hbond_motif_energy_sum,closest_autodock_recovery_rmsd,closest_autodock_recovery_ddg,strain_energy,total\n")

#look at all pdb files in the given directory and get started on scoring
for r,d,f in os.walk(location):
	for file in f:
		#only look in the current location (do not look further down)
		if r != location:
			continue

		#only look at pdb files, and make sure to ignore the minipose file
		if file.endswith(".pdb") and "minipose" not in file:
			#set initial values for score terms when scored and unscored
			#set as 2 entry list initially seeded as ""
			#entry 0 is the raw value and entry 1 is with the weight applied
			ddg = [0,0]
			total_motifs = [0,0]
			significant_motifs = [0,0]
			real_motif_ratio = [0,0]
			hbond_motif_count = [0,0]
			hbond_motif_energy_sum = [0,0]
			#raw starting value set different for selecting closest rmsd value
			closest_autodock_recovery_rmsd = [100,0]
			closest_autodock_recovery_ddg = [0,0]
			strain_energy = [0,0]
			total = [0,0]

			#track if all mandatory residues from mandatory_motif_residues are found
			#start true, turn false when checking
			all_mandatory_residues_found = True
			#list to hold all residues for motifs in placement
			found_motif_residues = []

			#make a directory that contains the starting pdb file and will hold any intermediate data
			dir_name = file.split(".")[0]
			os.system("mkdir " + dir_name)
			os.system("cp " + file + " " + dir_name)

			#extract comment table data from the file to potentially assign values for ddg, total_motifs, significant_motifs, and real_motif_ratio
			#read the file
			placement_file = open(r + "/" + file, "r")

			#get the data
			for line in placement_file.readlines():
				#ddg
				if line.startswith("Scoring: Post-HighResDock system ddG:"):
					ddg = [float(line.split()[len(line.split()) - 1].strip()),float(line.split()[len(line.split()) - 1].strip())*score_weights["ddg"]]
				#total_motifs
				if line.startswith("Placement motifs: Total motifs made:"):
					total_motifs = [float(line.split()[len(line.split()) - 1].strip()),float(line.split()[len(line.split()) - 1].strip())*score_weights["total_motifs"]]
				#significant_motifs
				if line.startswith("Placement motifs: Motifs made against significant residues count:"):
					significant_motifs = [float(line.split()[len(line.split()) - 1].strip()),float(line.split()[len(line.split()) - 1].strip())*score_weights["significant_motifs"]]
				#real_motif_ratio
				if line.startswith("Placement motifs: Real motif ratio:"):
					real_motif_ratio = [float(line.split()[len(line.split()) - 1].strip()),float(line.split()[len(line.split()) - 1].strip())*score_weights["real_motif_ratio"]]

				#pull hbond data 
				#motif data lines
				if ": Placement motif " in line:

					#determine the index of the residue and translate it to add to the motif list(s)
					index = line.split("Hbond_score")[1].split("_")[1][3:]

					#translate the index if in residue_index_dict keys
					if index in residue_index_dict.keys():
						translated_index = residue_index_dict[index]
						index = translated_index
					
					found_motif_residues.append(index)

					#get the hbond score, last piece of data on the line after a :, convert to float
					hbond_score = float(line.split(":")[len(line.split(":")) - 1].strip())

					#adjust non-zero score count
					if hbond_score != 0:
						hbond_motif_count[0] = hbond_motif_count[0] + 1
						hbond_motif_count[1] = hbond_motif_count[1] + 1
					#apply the hbond score to the score sum
					hbond_motif_energy_sum[0] = hbond_motif_energy_sum[0] + hbond_score
					hbond_motif_energy_sum[1] = hbond_motif_energy_sum[1] + hbond_score

			#apply the weights to the hbond motif count and total hbond energy
			hbond_motif_count[1] = hbond_motif_count[1] * score_weights["hbond_motif_count"]
			hbond_motif_energy_sum[1] = hbond_motif_energy_sum[1] * score_weights["hbond_motif_energy_sum"]

			#if the unweighted real motif ratio is below the cutoff, continue to the next ligand and ignore current
			if real_motif_ratio[0] < minimum_real_motif_ratio:
				print("Real motif ratio too low at: " + str(real_motif_ratio[0]))
				continue

			#check if the placement has motifs for all required residues
			#skip check if the mandatory_motifs_list is empty
			if mandatory_motifs_list != []:
				for res_id in mandatory_motifs_list:
					if res_id not in found_motif_residues:
						all_mandatory_residues_found = False

			#if not all mandatory residues have been found, continue and skip this ligand placement
			if all_mandatory_residues_found == False:
				print("Not all mandatory residues had motifs")
				continue

			#attempt strain torsion if able
			if attempt_strain:
				
				print("Attempting tldr STRAIN")

				#call the run_torsion_check_on_placed_ligands.py script
				print("python " + strain_automate_path + "run_torsion_check_on_placed_ligands.py " + r + "/" + dir_name + " " + strain_exec_path)
				os.system("python " + strain_automate_path + "run_torsion_check_on_placed_ligands.py " + r + "/" + dir_name + " " + strain_exec_path)

				#read the torsion data table csv
				#(pdb name)_torsion_data.csv
				strain_file = open(dir_name + "/" + dir_name + "_lig_Torsion_Strain.csv", "r")

				#read the strain file and get the total strain energy if it was calculated
				for line in strain_file.readlines():
					#there is a chance the the energy is inexplicably not calculated
					#behavior for this case
					#split(",")[1] of line is "NA" in this case
					if line.split(",")[1].strip() == "NA":
						#set strain energy to 100
						strain_energy = [100,100*score_weights["strain_energy"]]
					else:
						#get average of high and low energy predictions to use
						high = float(line.split(",")[1])
						low = float(line.split(",")[2])
						average = (high + low) / 2

						strain_energy = [average,average*score_weights["strain_energy"]]

			#if STRAIN was not able to get an energy score, kill and continue
			if strain_energy[0] == 100 and kill:
				print("Unable to derive torsion strain energy")
				continue

			#kill if the strain energy is too high
			if strain_energy[0] > strain_cutoff:
				print("Torsion strain energy too high at " + str(strain_energy[0]))
				continue

			#attempt autodock if able
			#attempt autodock last, as autodock is the slowest step
			if attempt_autodock:

				print("Attempting autodock")

				#call run_autodock_on_placed_ligands.py on the directory for this file
				#path to automation script, path to working file, path to autodock executable
				os.system("python " + autodock_automate_path + "run_autodock_on_placed_ligands.py " + r + "/" + dir_name + " " + autodock_exec_path)

				#read the corresponding autodock_data.csv file and determine the closest placement if there is any (and get corresponding energy score)
				#if there are no placement attempts, set values for closest_autodock_recovery_rmsd to 100 (arbitrary high value that should tank the placement because getting no placement attempts is bad) and closest_autodock_recovery_ddg to 0
				autodock_file = open(dir_name + "/autodock_data.csv", "r")

				for line in autodock_file.readlines():
					#ignore the header line
					if "system,model,rmsd,close_to_rosetta,vina_energy" in line:
						continue
					#handle if there were no placements; if split(,)[1] is 0
					if line.split(",")[1] == "0":
						closest_autodock_recovery_ddg = [0,0*score_weights["closest_autodock_recovery_ddg"]]
						closest_autodock_recovery_rmsd = [100,100*score_weights["closest_autodock_recovery_rmsd"]]
					#otherwise, handle the data extract the closest rmsd with its corresponding energy
					if float(line.split(",")[2]) < closest_autodock_recovery_rmsd[0]:
						#set the new closest rmsd with its ddg
						closest_autodock_recovery_rmsd = [float(line.split(",")[2]), float(line.split(",")[2]) * score_weights["closest_autodock_recovery_rmsd"]]
						closest_autodock_recovery_ddg = [float(line.split(",")[4].strip()), float(line.split(",")[4].strip()) * score_weights["closest_autodock_recovery_ddg"]]

			#if Autodock was not able to place, kill and continue
			if (closest_autodock_recovery_rmsd[0] == 100 or closest_autodock_recovery_ddg[0] == 0 ) and kill:
				print("AutoDock unable to place ligand")
				continue

			#if the autodock closest rmsd was too far, kill it
			if closest_autodock_recovery_rmsd[0] > autodock_rmsd_cutoff:
				print("Closest AutoDock recovery is beyond maximum cutoff at: " + str(closest_autodock_recovery_rmsd[0]))
				continue

			total[0] = ddg[0] + total_motifs[0] + significant_motifs[0] + real_motif_ratio[0] + closest_autodock_recovery_rmsd[0] + closest_autodock_recovery_ddg[0] + strain_energy[0] + hbond_motif_count[0] + hbond_motif_energy_sum[0]
			total[1] = ddg[1] + total_motifs[1] + significant_motifs[1] + real_motif_ratio[1] + closest_autodock_recovery_rmsd[1] + closest_autodock_recovery_ddg[1] + strain_energy[1] + hbond_motif_count[1] + hbond_motif_energy_sum[1]

			#write scores to the score files
			#file,ddg,total_motifs,significant_motifs,real_motif_ratio,hbond_motif_count,hbond_motif_energy_sum,closest_autodock_recovery_rmsd,closest_autodock_recovery_ddg,strain_energy,total
			raw_score_file.write(r + "/" + file + "," + str(ddg[0]) + "," + str(total_motifs[0]) + "," + str(significant_motifs[0]) + "," + str(real_motif_ratio[0]) + "," + str(hbond_motif_count[0]) + "," + str(hbond_motif_energy_sum[0]) + "," + str(closest_autodock_recovery_rmsd[0]) + "," + str(closest_autodock_recovery_ddg[0]) + "," + str(strain_energy[0]) + "," + str(total[0]) + "\n")
			weighted_score_file.write(r + "/" + file + "," + str(ddg[1]) + "," + str(total_motifs[1]) + "," + str(significant_motifs[1]) + "," + str(real_motif_ratio[1]) + "," + str(hbond_motif_count[0]) + "," + str(hbond_motif_energy_sum[0]) + "," + str(closest_autodock_recovery_rmsd[1]) + "," + str(closest_autodock_recovery_ddg[1]) + "," + str(strain_energy[1]) + "," + str(total[1]) + "\n")
