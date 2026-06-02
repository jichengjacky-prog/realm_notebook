import os,sys

compressed_file = sys.argv[1]

readfile = open(compressed_file, "r")
writefile = open("fixed_" + compressed_file, "w")

#read the read file and fix the spacing to write out
for line in readfile.readlines():
	if line.startswith("NAME"):
		#
		#print(line)
		writefile.write(line)
	elif line.startswith("IO_STRING"):
		#
		#print(line)
		writefile.write(line)
	elif line.startswith("TYPE"):
		#
		#print(line)
		writefile.write(line)
	elif line.startswith("AA"):
		#
		#print(line)
		writefile.write(line)
	elif line.startswith("ATOM"):
		#
		#print(line)

		#get line components, remove newline at end
		line_split = line[:-1].split()
		#print("ATOM %-4s %-4s %-4s %.2f\n" % (line_split[1], line_split[2], line_split[3], float(line_split[4])))
		writefile.write("ATOM %-4s %-4s %-4s %.2f\n" % (line_split[1], line_split[2], line_split[3], float(line_split[4])))
	elif line.startswith("BOND_TYPE"):
		#
		#print(line)

		#get line components, remove newline at end
		line_split = line[:-1].split()
		#print("BOND_TYPE %-4s %-4s %-4s\n" % (line_split[1], line_split[2], line_split[3]))
		writefile.write("BOND_TYPE %-4s %-4s %-4s\n" % (line_split[1], line_split[2], line_split[3]))
	
	elif line.startswith("CHI"):
		#
		#print(line)
		
		#get line components, remove newline at end
		line_split = line[:-1].split()
		#print("CHI %i %-4s %-4s %-4s %-4s\n" % (int(line_split[1]), line_split[2], line_split[3], line_split[4], line_split[5]))
		writefile.write("CHI %i %-4s %-4s %-4s %-4s\n" % (int(line_split[1]), line_split[2], line_split[3], line_split[4], line_split[5]))
	
	elif line.startswith("PROTON_CHI"):
		#
		#print(line)
		
		#get line components, remove newline at end
		#line_split = line[:-1].split()
		#print("ATOM %-4s %-4s %-4s %.2f\n" % (line_split[1], line_split[2], line_split[3], line_split[4]))
		#writefile.write("ATOM %-4s %-4s %-4s %.2f\n" % (line_split[1], line_split[2], line_split[3], line_split[4]))
		writefile.write(line)

	elif line.startswith("NBR_ATOM"):
		#
		#print(line)
		
		#get line components, remove newline at end
		#line_split = line[:-1].split()
		#print("ATOM %-4s %-4s %-4s %.2f\n" % (line_split[1], line_split[2], line_split[3], line_split[4]))
		#writefile.write("ATOM %-4s %-4s %-4s %.2f\n" % (line_split[1], line_split[2], line_split[3], line_split[4]))
		writefile.write(line)

	elif line.startswith("NBR_RADIUS"):
		#
		#print(line)
		writefile.write(line)
	elif line.startswith("ICOOR_INTERNAL"):
		#
		#print(line)
		
		#get line components, remove newline at end
		line_split = line[:-1].split()
		#print("ICOOR_INTERNAL   %-4s %11.6f %11.6f %11.6f  %-4s  %-4s  %-4s\n" % (line_split[1], float(line_split[2]), float(line_split[3]), float(line_split[4]), line_split[5], line_split[6], line_split[7]))
		writefile.write("ICOOR_INTERNAL   %-4s %11.6f %11.6f %11.6f  %-4s  %-4s  %-4s\n" % (line_split[1], float(line_split[2]), float(line_split[3]), float(line_split[4]), line_split[5], line_split[6], line_split[7]))
#	elif line.startswith(""):
		#
	else:
		print("Line not accounted for!")
		print(line)
		writefile.write(line)


"""
NAME
IO_STRING 
TYPE 
AA UNK
ATOM
BOND_TYPE
CHI
PROTON_CHI
NBR_ATOM
NBR_RADIUS
ICOOR_INTERNAL
"""
