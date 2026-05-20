#this script works with a placements directory that was dehydrated with shrink_placement_pdbs_to_placement_and_surrounding_residues.py, using a skeleton pdb to rebuild full pdbs 
#unless an optional argument is used to select a location, this script works in the location where it is called


import os,sys

working_location = os.getcwd()

if len(sys.argv) == 2:
	#get the different location
	working_location = sys.argv[1]

#move to the working location
os.chdir(working_location)

#unzip the placements.tar.gz
os.system("tar -xzf placements.tar.gz")

#move into the directory
os.chdir("placements")

placements_location = os.getcwd()

#read the reference pdb and get all residues for center of mass of heavy atoms (no hydrogens)
#doing center of mass shift because atom index names may not align between pdb given to rosetta and rosetta placement pdb
ref_file = open("skeleton.pdb","r")

#dictionary to hold the reference pdb residue centers of mass
#key is a tuple of the reference residue index (may not be the same as rosetta output) and residue 3 letter code
ref_res_data = {}

#have a working list of atom coordinates that builds as we find residue atoms
ref_res_atoms = []
cur_pair = ("","")

for line in ref_file.readlines():
	#if it is an ATOM line, we have something to work with
	if line.startswith("ATOM"):

		#get atom residue index and code
		atom_res_index = line.split()[5]
		atom_res_code = line.split()[3]

		#if starting with the first atom, change the current residue to the residue we just hit
		if cur_pair[0] == "" and cur_pair[1] == "":
			#cur_pair[0] = atom_res_index
			#cur_pair[1] = atom_res_code
			cur_pair = (atom_res_index,atom_res_code)

		#if the current atom is in the current residue, add it to the reference residue atoms
		if cur_pair[0] == atom_res_index:
			ref_res_atoms.append(line)
		#otherwise, if there is a mismatch, we hit the end of the residue and should derive the center of mass
		else:
			ref_res_data[cur_pair] = ref_res_atoms

			#now, set new cur_pair values since we are moving on and reset the ref_res_atoms list for the new residue
			#cur_pair[0] = atom_res_index
			#cur_pair[1] = atom_res_code
			cur_pair = (atom_res_index,atom_res_code)

			ref_res_atoms = []
			ref_res_atoms.append(line)

	#at the end when we hit the TER line to indicate end of protein/chain, get the com for the last residue
	if line.startswith("TER") and len(ref_res_atoms) > 0 and cur_pair[0] != "" and cur_pair[1] != "":
		ref_res_data[cur_pair] = ref_res_atoms

		ref_res_atoms = []

#now, go through every pdb (besides the skeleton) and rehydrate
for r,d,f in os.walk(placements_location):
	for file in f:
		if r == placements_location and file != "skeleton.pdb" and file.endswith(".pdb"):
			print(file)

			#pull its data the same as we did the skeleton

			read_file = open(file, "r")
			#open a file to write to to splice the skeleton and placement files together as we read the dehydrated placement
			temp_file = open("temp.pdb", "w")

			#dictionary to hold the reference pdb residue centers of mass
			#key is a tuple of the reference residue index (may not be the same as rosetta output) and residue 3 letter code
			placement_res_data = {}

			#have a working list of atom coordinates that builds as we find residue atoms
			ref_res_atoms = []
			cur_pair = ("","")

			#variable to note the last time that we spliced residues between the skeleton and placement, start at 0
			last_index_spliced = 0

			for line in read_file.readlines():
				#if it is an ATOM line, we have something to work with
				if line.startswith("ATOM"):

					#get atom residue index and code
					atom_res_index = line.split()[5]
					atom_res_code = line.split()[3]

					#if starting with the first atom, change the current residue to the residue we just hit
					if cur_pair[0] == "" and cur_pair[1] == "":
						#cur_pair[0] = atom_res_index
						#cur_pair[1] = atom_res_code
						cur_pair = (atom_res_index,atom_res_code)

					#if the current atom is in the current residue, add it to the reference residue atoms
					if cur_pair[0] == atom_res_index:
						ref_res_atoms.append(line)
					#otherwise, if there is a mismatch, we hit the end of the residue and should derive the center of mass
					else:
						placement_res_data[cur_pair] = ref_res_atoms

						#write all atoms of all residues from the skeleton (by index) up to this point and then write the atoms of this residue
						for residues in ref_res_data.keys():
							#read and print all skeleton residues up to the current residue in the placement and from the last time that we printed
							if int(residues[0]) < int(cur_pair[0]) and int(residues[0]) > last_index_spliced:
								for ref_line in ref_res_data[residues]:
									temp_file.write(ref_line)

						#write all of the residue atoms from the placement
						for placement_line in ref_res_atoms:
							temp_file.write(placement_line)

						#reset last_index_spliced to be the current index
						last_index_spliced = int(cur_pair[0])


						#now, set new cur_pair values since we are moving on and reset the ref_res_atoms list for the new residue
						#cur_pair[0] = atom_res_index
						#cur_pair[1] = atom_res_code
						cur_pair = (atom_res_index,atom_res_code)

						ref_res_atoms = []
						ref_res_atoms.append(line)

				#at the end when we hit the TER line to indicate end of protein/chain, get the com for the last residue
				if line.startswith("TER") and len(ref_res_atoms) > 0 and cur_pair[0] != "" and cur_pair[1] != "":
					placement_res_data[cur_pair] = ref_res_atoms

					#write all atoms of all residues from the skeleton (by index) up to this point and then write the atoms of this residue
					for residues in ref_res_data.keys():
						#read and print all skeleton residues up to the current residue in the placement and from the last time that we printed
						#if int(residues[0]) < int(cur_pair[0]) and int(residues[0]) > last_index_spliced:
						if int(residues[0]) > last_index_spliced:
							for ref_line in ref_res_data[residues]:
								temp_file.write(ref_line)

					#write all of the residue atoms from the placement
					for placement_line in ref_res_atoms:
						temp_file.write(placement_line)

					#reset last_index_spliced to be the current index
					#last_index_spliced = int(cur_pair[0])

					ref_res_atoms = []


				if line.startswith("ATOM") == False:
					#write the line if not dealing with atoms
					temp_file.write(line)

			#close the temporary file
			temp_file.close()

			#overwrite the original with the trimmed with the exact same name
			os.system("mv temp.pdb " + file)

#now, all files should be processed, recompress
#move back up to the level where the old compressed palcements file is
os.chdir(working_location)

os.system("tar -czf placements.tar.gz placements")

#remove the palcements_directory
os.system("rm -drf placements")


