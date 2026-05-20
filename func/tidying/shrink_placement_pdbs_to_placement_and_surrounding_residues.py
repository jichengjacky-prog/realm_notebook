#the purpose of this script is to help reduce memory overhead of placements by taking placement pdbs and only keeping residues that are within 5 angstroms of the placed ligand
#this requires a python environment with pymol2 (actually, maybe just go through and get distances via raw text so we can maintain residue indexing pdb comments)
#this is an individual script, which is meant to be operated in the directory that it is run in or a single command line argument for a path to process can be used
#for the purpose of pipeline cleanup, this will look for a file named placements.tar.gz, unpack the tar, make reduced versions of the placements, and then recompress them
#will also remove extraneous folders from the folder

import os,sys
#import pymol2

working_location = os.getcwd()

#arguments (pre-docked PDB for residue shift reference; optional path to go to operate in)
#reference pdb (this assumes that the placement system and the reference are aligned, which they should be if the placement was derived from a rosetta run with this placement with no other operations performed on it)
reference_pdb = sys.argv[1]

if len(sys.argv) == 3:
	#get the different location
	working_location = sys.argv[2]

#move to the working location
os.chdir(working_location)

#unzip the placements.tar.gz
os.system("tar -xzf placements.tar.gz")

#move into the directory
os.chdir("placements")

placements_location = os.getcwd()

#read the reference pdb and get all residues for center of mass of heavy atoms (no hydrogens)
#doing center of mass shift because atom index names may not align between pdb given to rosetta and rosetta placement pdb
ref_file = open(reference_pdb,"r")

#dictionary to hold the reference pdb residue centers of mass
#key is a tuple of the reference residue index (may not be the same as rosetta output) and residue 3 letter code
ref_res_com = {}

#have a working list of atom coordinates that builds as we find residue atoms
ref_res_atoms = []
cur_pair = ("","")

for line in ref_file.readlines():
	#if it is an ATOM line, we have something to work with
	if line.startswith("ATOM"):
		#if the atom is a hydrogen, skip (last character in line is H)
		if line.strip().endswith("H"):
			continue

		#get atom xyz
		atom_xyz = [float(line[30:38]),float(line[38:46]),float(line[46:54])]

		#get atom residue index and code
		atom_res_index = int(line[22:26])
		atom_res_code = line[17:20]

		#if starting with the first atom, change the current residue to the residue we just hit
		if cur_pair[0] == "" and cur_pair[1] == "":
			#cur_pair[0] = atom_res_index
			#cur_pair[1] = atom_res_code
			cur_pair = (atom_res_index,atom_res_code)

		#if the current atom is in the current residue, add it to the reference residue atoms
		if cur_pair[0] == atom_res_index:
			ref_res_atoms.append(atom_xyz)
		#otherwise, if there is a mismatch, we hit the end of the residue and should derive the center of mass
		else:
			#derive the sum of each dimension
			x_sum = 0
			y_sum = 0
			z_sum = 0
			for atom in ref_res_atoms:
				x_sum = x_sum + atom[0]
				y_sum = y_sum + atom[1]
				z_sum = z_sum + atom[2]

			#get the com coordinate by averaging the sum by number of atoms
			#assign to the dictionary
			ref_res_com[cur_pair] = [x_sum / len(ref_res_atoms), y_sum / len(ref_res_atoms), z_sum / len(ref_res_atoms)]

			#now, set new cur_pair values since we are moving on and reset the ref_res_atoms list for the new residue
			#cur_pair[0] = atom_res_index
			#cur_pair[1] = atom_res_code
			cur_pair = (atom_res_index,atom_res_code)

			ref_res_atoms = []
			ref_res_atoms.append(atom_xyz)

	#at the end when we hit the TER line to indicate end of protein/chain, get the com for the last residue
	if line.startswith("TER") and len(ref_res_atoms) > 0 and cur_pair[0] != "" and cur_pair[1] != "":
		#derive the sum of each dimension
		x_sum = 0
		y_sum = 0
		z_sum = 0
		for atom in ref_res_atoms:
			x_sum = x_sum + atom[0]
			y_sum = y_sum + atom[1]
			z_sum = z_sum + atom[2]

		#get the com coordinate by averaging the sum by number of atoms
		#assign to the dictionary
		ref_res_com[cur_pair] = [x_sum / len(ref_res_atoms), y_sum / len(ref_res_atoms), z_sum / len(ref_res_atoms)]

		ref_res_atoms = []		


#make a dictionary similar to prot_res that will be made, but in global scope that will make a skeleton of residues that do not move at least once, so rosetta palcements can be completely rebuilt
skeleton_res = {}


#iterate over the directory and look for placement pdb files
for r,d,f in os.walk(placements_location):
	#iterate over each pdb file
	for file in f:
		#only work on files in the placements location
		if r == placements_location and file.endswith(".pdb"):

			print(file)

			#if it is the minipose pdb, delete it and continue
			if "minipose" in file:
				os.system("rm -drf " + file)
				continue

			#skip a skeleton if one exists
			if file == "skeleton.pdb":
				continue

			#otherwise, delete the corresponding folder with the file's namesake first
			file_base = file.split(".pdb")[0]

			os.system("rm -drf " + file_base)

			#now, read the file so we can prune away distant residues
			#buffer all lines into a list so we only have to read the file once
			read_file = open(file,"r")

			file_lines = []

			#create a dictionary of residues (line starts with ATOM) with all atoms so we can see if the residue is close enough to any ligand (any res atom in the residue is within 5 angstroms of any lig atom)
			#use index as a key; values will be the atom line, which the xyz coordinates can be derived for distance
			prot_res = {}

			#make a list of ligand atom coordinates (no need to hold the lines since we will write them anyway)
			lig_res_atoms = []

			for line in read_file.readlines():
				file_lines.append(line)

				#if a residue atom, work on adding it
				if line.startswith("ATOM"):
					#add it to the protein residue dictionary
					#check if this is aready in the dicitonary or not and make a new list at the key if it does not exist already
					resnum = line[22:26]
					if resnum not in prot_res.keys():
						prot_res[resnum] = []

					prot_res[resnum].append(line)

				#if a hetatm (ligand) take the coordinates and add to the list
				if line.startswith("HETATM"):
					#append the x,y,z coords to lig_res_Atoms
					lig_res_atoms.append([float(line[30:38]),float(line[38:46]),float(line[46:54])])

			#we have now buffered in the file, determine which residues to keep based on distance
			#list of residues to keep by index
			res_to_keep = []

			
			for res in prot_res.keys():
				#keep any residues with any shift so we retain that data
				#determine the residue heavy atom center of mass
				res_com = [0,0,0]
				nheavy_atoms = 0
				#get residue 3 letter code
				residue_code = ""
				for res_atom in prot_res[res]:
					#skip hydrogens
					if res_atom.strip().endswith("H"):
						continue

					#extract the residue atom coordinates to add to the res com
					#res_com[0] = res_com[0] + float(res_atom.split()[6])
					#res_com[1] = res_com[1] + float(res_atom.split()[7])
					#res_com[2] = res_com[2] + float(res_atom.split()[8])

					res_com[0] = res_com[0] + float(res_atom[30:38])
					res_com[1] = res_com[1] + float(res_atom[38:46])
					res_com[2] = res_com[2] + float(res_atom[46:54])
					nheavy_atoms = nheavy_atoms + 1

					if residue_code == "":
						residue_code = res_atom[17:20]

				#divide the residue com sums by the number of atoms to get a usable average
				res_com[0] = res_com[0] / nheavy_atoms
				res_com[1] = res_com[1] / nheavy_atoms
				res_com[2] = res_com[2] / nheavy_atoms

				shortest_distance = 100

				#run through each residue in the reference and determine which residue is the closest by com and then what the distance is (if greater than 0, we should keep; may have to test and see if there are rounding issues)
				for ref_res in ref_res_com.keys():
					#if the residue code matches the reference residue, consider for comparison
					if residue_code == ref_res[1]:
						#get the distance between the reference and placement residue coms
						distance = (((ref_res_com[ref_res][0] - res_com[0]) ** 2) + ((ref_res_com[ref_res][1] - res_com[1]) ** 2) + ((ref_res_com[ref_res][2] - res_com[2]) ** 2)) ** 0.5

						#check if the distance is shorter than the shortest recorded distance, and if so, replace (no need to track the residue)
						#if the distance is 0, go ahead and break since we don't need to check with any residues
						if distance < shortest_distance:
							shortest_distance = distance

						if shortest_distance == 0:
							break

				#if the shortest distance is greater than 0.00001 (meaning that the residue moved any amount, but accounting for minor rounding errors), add it to residues to keep and continue, skipping proximity keep
				if shortest_distance > 0.00001:
					#debug print of the residue name/index and distance
					print(shortest_distance, res, residue_code)
					res_to_keep.append(res)
					continue
				else:
					#if the residue was close enough to the reference, try to add its data to the skeleton if the residue index is not present
					#keep regardless of proximity to the ligand
					#remember, the data here is a list of file lines for easy writing later
					if res not in skeleton_res.keys():
						skeleton_res[res] = prot_res[res]

				
				#look for residues within 5 angstroms to keep
				for res_atom in prot_res[res]:
					#extract the residue atom coordinates
					x = float(res_atom[30:38])
					y = float(res_atom[38:46])
					z = float(res_atom[46:54])

					#get the coordinate distance of the atom to each atom in the ligand
					for lig_atom in lig_res_atoms:
						distance = (((lig_atom[0] - x) ** 2) + ((lig_atom[1] - y) ** 2) + ((lig_atom[2] - z) ** 2)) ** 0.5

						#if the distance is within 5, note the residue in residues to keep
						if distance < 5:
							res_to_keep.append(res)
							#break
							break

					#if we added the residue to the res_to_keep, break here to move to the next residue
					if res in res_to_keep:
						break

			

			#finally, go back over the buffered file, and when ew get to ATOM lines, determine whether or not to keep the residue based on if the residue is in the res_to_keep
			write_file = open("temp.pdb", "w")

			for line in file_lines:
				#if the line does not start with ATOM, write
				if line.startswith("ATOM") == False:
					write_file.write(line)
				else:
					#if it is a residue atom line, determine the residue index and determine whether we identified if it iwas close enough
					if int(line[22:26]) in res_to_keep:
						write_file.write(line)

			#close the read and write streams
			read_file.close()
			write_file.close()

			#overwrite the original with the trimmed with the exact same name
			os.system("mv temp.pdb " + file)

#finally, write a file of the receptor skeleton containing data on residues that did not move at least once that is indexed with rosetta residue indexing with atoms in reference coordinates
skeleton_file = open("skeleton.pdb","w")
for res in skeleton_res.keys():
	for line in skeleton_res[res]:
		skeleton_file.write(line)

#cap iff file with ter and end
skeleton_file.write("TER\nEND\n")
skeleton_file.close()

#now, all files should be processed, recompress
#move back up to the level where the old compressed palcements file is
os.chdir(working_location)

os.system("tar -czf placements.tar.gz placements")

#remove the palcements_directory
os.system("rm -drf placements")



