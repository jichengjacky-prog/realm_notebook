#the purpose of this script is to look at the input PDB for Rosetta and a user-specified placements file and create a csv key file that notse the residue index shift between the original and placement
#this is needed because Rosetta changes residue indices of the residues in the placement files
#if we want to visualize reisdues with real vs fake motifs from the motifs library indicated in the placements pdb data with pymol, we need to be able to translate the original vs translated indices with a key file that the pymol generation script reads
#this assumes that the rosetta input pdb and the output placement are aligned (which they should be)

#using the assumption of aligned pdbs, residues indices are figured out by matching the closest residues by rmsd of center of mass (which should be very close to 0 as residues barely shift) and noting the indices 

import os,sys

#args

#rosetta input pdb (which will match the indexing of motif residue indices in the placement file)
starting_file = sys.argv[1]
sf_base = starting_file.split(".pdb")[0].split("/")[-1]

#rosetta placement pdb (which would have the motifs data in it that we translate to, but we do not use in this script)
placement_file = sys.argv[2]
pf_base = placement_file.split(".pdb")[0].split("/")[-1]

#set a output key file name based on the names of the starting files
out_file_name = sf_base + "_vs_" + pf_base + "_key.csv"

#optional, if the user wants a likely cleaner file name, they can input a 3rd arg (this can also be used to indicate outputting to a different path)
if len(sys.argv) == 4:
	out_file_name = sys.argv[3]

	#if the out file name does not end in .csv, make it end in .csv
	if out_file_name.endswith(".csv") == False:
		out_file_name = out_file_name + ".csv"

#start writing the output key file
key_file = open(out_file_name,"w")
#write the header, order is residue type, target index, placement index, difference in indices
key_file.write("res_type," + str(sf_base) + "," + str(pf_base) + ",delta\n")

#now, read through both files, find matching residues, and note their indices in the key
#list of all residues in the starting file
#the lists will consist of the residue type, residue index, and center of mass in x,y,z coordinates
#this will be used to match against an iterated residue from the placement file to determine the closest match of the same type
#the list is likely easier to deal with than a dictionary as we just want to iterate upon (although it could be made a little more efficient by binning by residue type, this simple algo should be fast and only run in a few seconds at worst)
sf_residues_list = []

#read the starting file
sf_read = open(starting_file,"r")
pf_read = open(placement_file,"r")

#list to build the working residue by atoms
working_residue_list = []
working_residue_number = -1

for line in sf_read.readlines():
	#only look at lines starting with "ATOM "
	if line.startswith("ATOM ") == False:
		continue

	#read the atom data
	atom_res_name = line.split()[3]
	atom_res_index = int(line.split()[5])
	atom_x = float(line.split()[6])
	atom_y = float(line.split()[7])
	atom_z = float(line.split()[8])

	#change the workign residue number of it is -1
	if working_residue_number == -1:
		working_residue_number = atom_res_index

	#add the atom to the workign residue list if the working residue is the same as the residie index of the current atom, otherwise cap off the residue and prepare to add it to the residue list
	if atom_res_index == working_residue_number:
		working_residue_list.append([atom_res_name,atom_res_index,atom_x,atom_y,atom_z])
	else:
		#if we are moving on to a new residue, wrap up the previous
		#get the prior residue name and index
		prior_name = working_residue_list[0][0]
		prior_index = working_residue_list[0][1]

		prior_x_center = 0
		prior_y_center = 0
		prior_z_center = 0

		for atom in working_residue_list:
			prior_x_center = prior_x_center + (atom[2] / len(working_residue_list))
			prior_y_center = prior_y_center + (atom[3] / len(working_residue_list))
			prior_z_center = prior_z_center + (atom[4] / len(working_residue_list))

		#we now have the residue centroid and can set the residue to the sf_residues_list
		sf_residues_list.append([prior_name,prior_index,prior_x_center,prior_y_center,prior_z_center])

		#wipe the workign residue and start with the new atom
		working_residue_list = []
		working_residue_number = atom_res_index
		working_residue_list.append([atom_res_name,atom_res_index,atom_x,atom_y,atom_z])

#cap off the last atom from the starting file
prior_name = working_residue_list[0][0]
prior_index = working_residue_list[0][1]

prior_x_center = 0
prior_y_center = 0
prior_z_center = 0

for atom in working_residue_list:
	prior_x_center = prior_x_center + (atom[2] / len(working_residue_list))
	prior_y_center = prior_y_center + (atom[3] / len(working_residue_list))
	prior_z_center = prior_z_center + (atom[4] / len(working_residue_list))

#we now have the residue centroid and can set the residue to the sf_residues_list
sf_residues_list.append([prior_name,prior_index,prior_x_center,prior_y_center,prior_z_center])

#now, look through all atoms in the placement file

#list to build the working residue by atoms
working_residue_list = []
working_residue_number = -1

for line in pf_read.readlines():
	#only look at lines starting with "ATOM "
	if line.startswith("ATOM ") == False:
		continue

	#read the atom data
	atom_res_name = line.split()[3]
	atom_res_index = int(line.split()[5])
	atom_x = float(line.split()[6])
	atom_y = float(line.split()[7])
	atom_z = float(line.split()[8])

	#change the workign residue number of it is -1
	if working_residue_number == -1:
		working_residue_number = atom_res_index

	#add the atom to the workign residue list if the working residue is the same as the residie index of the current atom, otherwise cap off the residue and prepare to add it to the residue list
	if atom_res_index == working_residue_number:
		working_residue_list.append([atom_res_name,atom_res_index,atom_x,atom_y,atom_z])
	else:
		#if we are moving on to a new residue, wrap up the previous
		#get the prior residue name and index
		prior_name = working_residue_list[0][0]
		prior_index = working_residue_list[0][1]

		prior_x_center = 0
		prior_y_center = 0
		prior_z_center = 0

		for atom in working_residue_list:
			prior_x_center = prior_x_center + (atom[2] / len(working_residue_list))
			prior_y_center = prior_y_center + (atom[3] / len(working_residue_list))
			prior_z_center = prior_z_center + (atom[4] / len(working_residue_list))

		#we now have the residue centroid and can match against the sf residues
		#list with notes the residue name,matching index from sf, and the distance
		closest_residue_match = [prior_name,-1,1000000]

		for res in sf_residues_list:
			#determine if the residue types are the same
			if prior_name != res[0]:
				continue

			#otherwise, if the residue names match, get the distance of centroids
			res_centroid_distance = (((prior_x_center - res[2]) ** 2) + ((prior_y_center - res[3]) ** 2) + ((prior_z_center - res[4]) ** 2)) ** 0.5

			#check if this is smaller than the shortest distance recorded so far, and if so, override
			if res_centroid_distance < closest_residue_match[2]:
				closest_residue_match = [prior_name,res[1],res_centroid_distance]

		#now that we have the closest match (if any), write to teh key file
		#skip if the index is -1
		if closest_residue_match[1] != -1:
			#key_file.write("res_type," + str(sf_base) + "," + str(pf_base) + ",delta\n")
			print("WRITING", prior_name, prior_index, "matched_to", closest_residue_match)
			key_file.write(str(closest_residue_match[0]) + "," + str(closest_residue_match[1]) + "," + str(prior_index) + "," + str(prior_index - closest_residue_match[1]) + "\n")


		#wipe the workign residue and start with the new atom
		working_residue_list = []
		working_residue_number = atom_res_index
		working_residue_list.append([atom_res_name,atom_res_index,atom_x,atom_y,atom_z])

#handle the final residue
prior_name = working_residue_list[0][0]
prior_index = working_residue_list[0][1]

prior_x_center = 0
prior_y_center = 0
prior_z_center = 0

for atom in working_residue_list:
	prior_x_center = prior_x_center + (atom[2] / len(working_residue_list))
	prior_y_center = prior_y_center + (atom[3] / len(working_residue_list))
	prior_z_center = prior_z_center + (atom[4] / len(working_residue_list))

#we now have the residue centroid and can match against the sf residues
#list with notes the residue name,matching index from sf, and the distance
closest_residue_match = [prior_name,-1,1000000]

for res in sf_residues_list:
	#determine if the residue types are the same
	if prior_name != res[0]:
		continue

	#otherwise, if the residue names match, get the distance of centroids
	res_centroid_distance = (((prior_x_center - res[2]) ** 2) + ((prior_y_center - res[3]) ** 2) + ((prior_z_center - res[4]) ** 2)) ** 0.5

	#check if this is smaller than the shortest distance recorded so far, and if so, override
	if res_centroid_distance < closest_residue_match[2]:
		closest_residue_match = [prior_name,res[1],res_centroid_distance]

#now that we have the closest match (if any), write to teh key file
#skip if the index is -1
if closest_residue_match[1] != -1:
	#key_file.write("res_type," + str(sf_base) + "," + str(pf_base) + ",delta\n")
	print("WRITING", prior_name, prior_index, "matched_to", closest_residue_match)
	key_file.write(str(closest_residue_match[0]) + "," + str(closest_residue_match[1]) + "," + str(prior_index) + "," + str(prior_index - closest_residue_match[1]) + "\n")