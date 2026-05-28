#the purpose of this script is a post-hoc fix (ideally, only used for the pth2 analysis and I fix before we run discovery again) to fix the paths in the refined conformers compiled files
#the existing paths in these files are the tmp locations where discovery was run, which is not helpful
#this will look down the current location at all raw_scores.csv files and fix the pathing
#as a bonus to fixing the old files, this will also write a compiled file at the working location containing all data in the raw_scores.csv files found

#imports
import os,sys

#start writing an all raw scores file
all_scores_file = open("all_raw_scores.csv", "w")

#write a header to the file (hardcoded)
all_scores_file.write("file,ddg,total_motifs,significant_motifs,real_motif_ratio,hbond_motif_count,hbond_motif_energy_sum,closest_autodock_recovery_rmsd,closest_autodock_recovery_ddg,strain_energy,total\n")

#record the starting location
starting_dir = os.getcwd()

#have a counter to print out for progress tracking
file_counter = 0

#look down from the current location and work with all raw_scores.csv
for r,d,f in os.walk(starting_dir):
	for file in f:
		if file == "raw_scores.csv":
			#move to that location
			os.chdir(r)

			#print the working file
			file_counter = file_counter + 1

			print(file_counter, r + "/" + file)

			#rename the raw score sfile to an old raw scores file
			os.system("mv raw_scores.csv old_raw_scores.csv")

			#read the old raw scores file and start writing to a new raw_scores_file
			old_file = open("old_raw_scores.csv", "r")
			new_file = open("raw_scores.csv", "w")

			for line in old_file.readlines():
				#bad tmp path line
				if line.startswith("/tmp/"):
					#split the line on / and take the last entry
					good_data = line.split("/")[-1]

					#write a fixed line based on the root (and add /placement/)
					new_file.write(r + "/placements/" + good_data)

					#also write to the all file
					all_scores_file.write(r + "/placements/" + good_data)

				#non-tmp line (header or otherwise good data to not fix, just write)
				else:
					new_file.write(line)

			#once done, close the files
			old_file.close()
			new_file.close()
