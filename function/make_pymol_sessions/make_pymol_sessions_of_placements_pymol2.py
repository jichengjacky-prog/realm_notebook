#note, this script was initially written by ChatGPT to write most of the pymol logic, and further modified by me (Ari) to enhance use of inputs

#usage:
#python make_pymol_sessions_of_placements_pymol2.py /path/to/placements/directory/ ligandresidueindex list+string+of+residue+indices
#python make_pymol_sessions_of_placements_pymol2.py /scratch/abgvg9/discovery_results/top_1000_placement/agonist_12M_passing_placements/0 282 227+86+253+63+257

#this uses pymol2, which may be better
import os,sys
import pymol2

# Define the directory containing your files as a command line argument
#directory = 'path/to/your/files'
directory = sys.argv[1]

#have the directory end with a backslash if it doesn't from the input
if directory.endswith("/") == False:
    directory = directory + "/"

#ligand residue index of interest
ligand_residue = sys.argv[2]

# Define the specific residues you want to highlight (e.g., 86, 227, 342)
highlight_residues = sys.argv[3]

#define the path to an optional residue index key csv file
use_key = False

#blank dictionary to hold residue index keys
#this is needed because at the very least, the indexing of residues on motifs collected off of 7l1u (agonists collection) do not match the index of the residue in the produce pdb
#this key helps to translate these mismatched residues (with the help of pdb 4s0v in this case)
#translation was performed by hand before running this script
#line structure is: res_type,original_index,translated_index,difference_in_index (optional)
#mostly concerned with the original and translated indices, but the other 2 can be helpful when translating 
residue_index_dict = {}
if len(sys.argv) == 5:
    use_key = True
    #get the key file and read it
    key_file = open(sys.argv[4],"r")

    for line in key_file.readlines():
        #skip the first line, starts with res_type
        if line.startswith("res_type"):
            continue

        #seed the dictionary with the original and translated indices
        residue_index_dict[line.split(",")[1]] = line.split(",")[2]


# Define the color for the entire protein
protein_color = 'cyan'

# Define the distance threshold (in angstroms)
distance_threshold = 5.0

# Define the selection around which to find neighboring residues (e.g., ligand)
ligand_selection = "resi " + str(ligand_residue)  # Adjust as necessary for your files

# Initialize PyMOL in headless mode
with pymol2.PyMOL() as pymol:
    
    #set the internal gui width
    pymol.cmd.set('internal_gui_width', 600)

    # Set the sphere scale to 0.25 for visualizing motifs
    pymol.cmd.set('sphere_scale', 0.5)

    items = sorted(os.listdir(directory))

    # Loop through each file in the directory
    for filename in items:
        if filename.endswith('.pdb'):  # Ensure only PDB files are processed
            filepath = os.path.join(directory, filename)
            
            # Load the file
            pymol.cmd.load(filepath)
            
            # Get the object name (assumed to be the filename without the extension)
            object_name = os.path.splitext(filename)[0]
            print(object_name)

            #read the file for real motif data and determine which motifs are real and which are not and get their indices

            # Color the entire protein
            pymol.cmd.color(protein_color, object_name)

            
            # Find and show residues within a certain distance of the ligand
            neighboring_residues = f'(byres {ligand_selection} around {distance_threshold}) and {object_name}'
            pymol.cmd.show('sticks', neighboring_residues)
            #pymol.cmd.color('elem', neighboring_residues)
            pymol.cmd.color('yellow', f'{neighboring_residues} and elem C')

            # Show the specific residues as sticks and color them
            pymol.cmd.show('sticks', f'{object_name} and resi {highlight_residues}')
            #pymol.cmd.color('elem', f'{object_name} and resi {highlight_residues}')
            pymol.cmd.color('orange', f'{object_name} and resi {highlight_residues} and elem C')

            #only do this if we should use the key
            if use_key:
                
                #list to hold all motif indices
                all_motifs = []

                #list to hold real motif indices
                real_motifs = []

                #read the pdb file and get the real motif data
                #make spheres for where all motifs are
                #determine which motifs are considered real and color them magenta
                pdb_file = open(filepath,"r")
                for line in pdb_file.readlines():
                    #real motif data lines
                    if "Real motif check" in line:
                            #determine if motif is real based on if "No real match" is in the line
                            is_real = True
                            if "No real match" in line:
                                is_real = False

                            #determine the index of the residue and translate it to add to the motif list(s)
                            index = line.split("Hbond_score")[1].split("_")[1][3:]

                            #translate the index
                            if index in residue_index_dict.keys():

                                translated_index = residue_index_dict[index]

                                all_motifs.append(translated_index)

                                if is_real:
                                    real_motifs.append(translated_index)
                            else:
                                #handling to not crash the program if we hit a residue not in the index dictionary key:
                                print("Warning, index " + index + " not found in the translation key!")

                #once done getting all motifs, make spheres on the residues in all_motifs, and then color the real motifs residues magenta
                #make selection strings for use with pymol
                all_motifs_string = "resi "
                for motif in all_motifs:
                     all_motifs_string = all_motifs_string + motif + "+"

                #remove + at end
                all_motifs_string = all_motifs_string[:-1]

                #repeat for real
                real_motifs_string = "resi "
                for motif in real_motifs:
                    real_motifs_string = real_motifs_string + motif + "+"
                real_motifs_string = real_motifs_string[:-1]

                #show spheres for all motifs and then color the real motifs
                pymol.cmd.show('spheres', f'{object_name} and {all_motifs_string}')
                pymol.cmd.color('brown', f'{object_name} and {all_motifs_string} and elem C')
                pymol.cmd.color('magenta', f'{object_name} and {real_motifs_string} and elem C')

            #color the ligands white
            pymol.cmd.color('white', f'{ligand_selection} and elem C')

            #fixing element coloring:
            pymol.cmd.color('white', 'elem H')
            pymol.cmd.color('red', 'elem O')
            pymol.cmd.color('blue', 'elem N')
            pymol.cmd.color('green', 'elem Cl')
            pymol.cmd.color('cyan', 'elem F')
            pymol.cmd.color('brown', 'elem Br')
            pymol.cmd.color('purple', 'elem I')
            pymol.cmd.color('yellow', 'elem S')

            #

            # Display hydrogen bonds
            pymol.cmd.dist(f'{object_name}_hbonds', f'{object_name} and {ligand_selection}', neighboring_residues, cutoff=3.5, mode=2)
            pymol.cmd.set('dash_color', 'yellow', f'{object_name}_hbonds')
            #pymol.cmd.dist(f'{object_name}_hbonds', f'{object_name} and {ligand_selection}', f'resi {highlight_residues}', cutoff=3.5, mode=2)
            #pymol.cmd.set('dash_color', 'green', f'{object_name}_hbonds')
            pymol.cmd.set('dash_width', 2.0)


             # Print names of all active objects
            active_objects = pymol.cmd.get_names()
            #print(f"Active objects in the session: {active_objects}")

            #attempt to count the number of measurements om the hydrogen bond distance object
            #distance = pymol.cmd.get_distance(f'{object_name}_hbonds')
            #print(distance)

            # Group the molecule and the hydrogen bond distance objects together
            pymol.cmd.group(f'{object_name}_group', f'{object_name} {object_name}_hbonds')


    # Save the session for all proteins
    pymol.cmd.save('all_proteins_session.pse')  # Save the PyMOL session


# Initialize PyMOL in headless mode for hbonds only
with pymol2.PyMOL() as pymol:
    
    #set the internal gui width
    pymol.cmd.set('internal_gui_width', 600)

    # Loop through each file in the directory
    for filename in os.listdir(directory):
        if filename.endswith('.pdb'):  # Ensure only PDB files are processed
            filepath = os.path.join(directory, filename)
            
            # Load the file
            pymol.cmd.load(filepath)
            
            # Get the object name (assumed to be the filename without the extension)
            object_name = os.path.splitext(filename)[0]
            print(object_name)
            
            # Color the entire protein
            pymol.cmd.color(protein_color, object_name)

            
            # Find and show residues within a certain distance of the ligand
            neighboring_residues = f'(byres {ligand_selection} around {distance_threshold}) and {object_name}'
            pymol.cmd.show('sticks', neighboring_residues)
            #pymol.cmd.color('elem', neighboring_residues)
            pymol.cmd.color('yellow', f'{neighboring_residues} and elem C')

            # Show the specific residues as sticks and color them
            pymol.cmd.show('sticks', f'{object_name} and resi {highlight_residues}')
            #pymol.cmd.color('elem', f'{object_name} and resi {highlight_residues}')
            pymol.cmd.color('orange', f'{object_name} and resi {highlight_residues} and elem C')

            #color the ligands magenta
            pymol.cmd.color('white', f'{ligand_selection} and elem C')

            #fixing element coloring:
            pymol.cmd.color('white', 'elem H')
            pymol.cmd.color('red', 'elem O')
            pymol.cmd.color('blue', 'elem N')
            pymol.cmd.color('green', 'elem Cl')
            pymol.cmd.color('cyan', 'elem F')
            pymol.cmd.color('brown', 'elem Br')
            pymol.cmd.color('purple', 'elem I')
            pymol.cmd.color('yellow', 'elem S')

            # Display hydrogen bonds
            pymol.cmd.dist(f'{object_name}_hbonds', f'{object_name} and {ligand_selection}', neighboring_residues, cutoff=3.5, mode=2)
            pymol.cmd.set('dash_color', 'green', f'{object_name}_hbonds')
            pymol.cmd.set('dash_width', 2.0)

             # Print names of all active objects
            active_objects = pymol.cmd.get_names()
            #print(f"Active objects in the session: {active_objects}")

            #attempt to count the number of measurements om the hydrogen bond distance object
            #distance = pymol.cmd.get_distance(f'{object_name}_hbonds')
            #print(distance)

            # Group the molecule and the hydrogen bond distance objects together
            #pymol.cmd.group(f'{object_name}_group', f'{object_name} {object_name}_hbonds')

            #delete the molecule object so we only keep the hbond
            pymol.cmd.delete(f'{object_name}')



    # Save the session for all proteins
    pymol.cmd.save('hbonds_only_session.pse')  # Save the PyMOL session
    