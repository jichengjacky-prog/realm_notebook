import os,sys

#going to try to use a heap data structure to better keep the best scores
import heapq

#value to hold the number of ligands to keep
max_ligands_to_keep = int(sys.argv[1])

#location of shapedb data
shapedb_data_location = sys.argv[2]

#sanitize path to end with a / if it doesn't already
if shapedb_data_location.endswith("/") == False:
	shapedb_data_location = shapedb_data_location + "/"

#run over the whole library from 0 to 53084
min_chunk = 0
min_chunk = int(sys.argv[3])
max_chunk = 10000
max_chunk = int(sys.argv[4])

#create the heap
#each entry will be a tuple of 4 entries: shapedb score, ligand name (with conformer number), chunk, subchunk
conformer_list = []

#heapq.heapify(conformer_list)

#value to hold the currently smallest value of the conformer list so that we don't have to keep calling it every time we want it if the value has not changed
#default to -2, because the value should never get that low
#smallest_value = -2

#create empty heap
conformer_heap = []


#i.e. 0 will have the range 0-99, 1 is 100-199, 10 is 1000-1099, 530 is 53,000-53,099 (this will be specifically cut to 53084 since that is the max)
#make strings for the min and max chunk
min_chunk_str = str(min_chunk)
while len(min_chunk_str) < 5:
	min_chunk_str = "0" + min_chunk_str

max_chunk_str = str(max_chunk)
while len(max_chunk_str) < 5:
	max_chunk_str = "0" + max_chunk_str

#iterate through each chunk and its subchunks of ligands to pull the best scoring ligands for the given shape file prefix
for i in range(min_chunk,max_chunk + 1):
	
	#if the current chunk (i) is greater than 53084 (highest chunk), break the loop
	if i > max_chunk:
		print("We are exceeding the size of the library, breaking the loop.")
		break

	try:
		#construct a 5 digit string to access chunks (lead with zeroes)
		chunk_str = str(i)
		print("On chunk " + chunk_str)

		while len(chunk_str) < 5:
			chunk_str = "0" + chunk_str

		#derive the superchunk that this belongs to for accession
		superchunk_str = str(int(chunk_str[0:3]))

		#iterate through each subchunk (0-9)
		for j in range(0,10):
			
			subchunk_str = str(j)

			#open the text file to read
			read_file_path = shapedb_data_location + superchunk_str + "/" + chunk_str + "/" + chunk_str + "_" + subchunk_str + "_nn_filtered.txt"
			try:
				read_file = open(read_file_path, "r")
			except FileNotFoundError:
				print(read_file_path + " is missing")
				continue
			except Exception as e:
				print("Error opening " + read_file_path + ": " + str(e))
				continue

			#read through the file
			for line in read_file:
				parts = line.split()
				if len(parts) < 2:
					continue

				#get conformer name and score
				#shapedb scores are multiplied by negative 1 so that the worse scoring placements are what gets popped off (since python heaps pop the lowest value)
				try:
					conf_score = float(parts[1].strip()) * -1
				except ValueError:
					continue
				conf_name = parts[0]

				#create entry tuple
				entry = (conf_score, conf_name, chunk_str, subchunk_str)

				# keep adding entries until we are at the max to keep
				if len(conformer_heap) < max_ligands_to_keep:
					conformer_heap.append(entry)
					# once at the max, heapify
					if len(conformer_heap) == max_ligands_to_keep:
						heapq.heapify(conformer_heap)
				else:
					# heap[0] is the worst (smallest) score we are currently keeping
					if conf_score > conformer_heap[0][0]:
						heapq.heapreplace(conformer_heap, entry)

			#close the file
			read_file.close()
	except Exception as e:
		print("Error processing chunk " + str(chunk_str) + ": " + str(e))
	#even smarter approach, only heapify at the end of the file after adding everything from the file into the list, and then pop off the lowest if we are beyond the max (and don't even heapify if we are below the max)
	#even more smart, only do this step when completing a chunk, not upon completion of each subchunk
	#if len(conformer_list) > max_ligands_to_keep:
	if len(conformer_list) > (max_ligands_to_keep * 2):
		conformer_heap = conformer_list

		print("Heapifying list after end of current chunk, " + str(chunk_str))
		heapq.heapify(conformer_heap)

		#use the nlargest command to get max number back as a list to use for the next file
		conformer_list = heapq.nlargest(max_ligands_to_keep,conformer_heap)

		#pop off until the length is at the max
		#while len(conformer_heap) > max_ligands_to_keep:
		#	heapq.heappop()
	
"""
#write to conformer_list and sort
conformer_list = sorted(conformer_heap, reverse=True)

#write the best ligands to a file with accession information to a new file

#make a folder if the source location to store this data (clobber existing data)
os.system("mkdir -p " + shapedb_data_location + "best_" + str(max_ligands_to_keep) + "_chunk_lists")

write_file = open(shapedb_data_location + "best_" + str(max_ligands_to_keep) + "_chunk_lists/best_" + str(max_ligands_to_keep) + "_chunks_" + str(min_chunk_str) + "_" + str(max_chunk_str) + ".txt","w")

#print(list(conformer_list[i]))
#print("-------------")
#print()

for i in conformer_list:
	write_file.write(str(i[0]) + "," + str(i[1]) + "," + str(i[2]) + "," + str(i[3]) + "\n")

write_file.close()
"""