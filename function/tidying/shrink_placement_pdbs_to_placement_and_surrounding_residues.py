#!/usr/bin/env python3
"""
Shrink (dehydrate) placement PDBs to only the placed ligand and surrounding
protein residues.  This dramatically reduces file sizes for storage.

For each placement PDB, keeps:
  - Residues within 5 Å of any ligand atom
  - Residues that moved relative to the reference (pre-docked) structure

Writes a skeleton.pdb of all unmoved residues so placements can be fully
reconstructed later (rehydration).

Usage:
    python shrink_placement_pdbs_to_placement_and_surrounding_residues.py <reference_pdb>
    python shrink_placement_pdbs_to_placement_and_surrounding_residues.py <reference_pdb> <working_dir>

Operates on placements.tar.gz in the working directory (default: cwd).
"""

import os
import shutil
import sys


def build_ref_residue_com(reference_pdb):
    """Build a dictionary of {(res_index, res_code): [cx, cy, cz]} centre-of-mass
    for each residue in the reference PDB (heavy atoms only)."""
    ref_res_com = {}
    cur_pair = ("", "")
    ref_res_atoms = []

    with open(reference_pdb, "r") as ref_file:
        for line in ref_file:
            if not line.startswith("ATOM"):
                if line.startswith("TER") and ref_res_atoms and cur_pair[0] != "":
                    n = len(ref_res_atoms)
                    cx = sum(a[0] for a in ref_res_atoms) / n
                    cy = sum(a[1] for a in ref_res_atoms) / n
                    cz = sum(a[2] for a in ref_res_atoms) / n
                    ref_res_com[cur_pair] = [cx, cy, cz]
                    ref_res_atoms = []
                continue

            if line.strip().endswith("H"):
                continue

            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            res_index = int(line[22:26])
            res_code = line[17:20]

            if cur_pair[0] == "":
                cur_pair = (res_index, res_code)

            if cur_pair[0] == res_index:
                ref_res_atoms.append([x, y, z])
            else:
                if ref_res_atoms:
                    n = len(ref_res_atoms)
                    cx = sum(a[0] for a in ref_res_atoms) / n
                    cy = sum(a[1] for a in ref_res_atoms) / n
                    cz = sum(a[2] for a in ref_res_atoms) / n
                    ref_res_com[cur_pair] = [cx, cy, cz]
                cur_pair = (res_index, res_code)
                ref_res_atoms = [[x, y, z]]

    # Don't forget the last residue
    if ref_res_atoms and cur_pair[0] != "":
        n = len(ref_res_atoms)
        cx = sum(a[0] for a in ref_res_atoms) / n
        cy = sum(a[1] for a in ref_res_atoms) / n
        cz = sum(a[2] for a in ref_res_atoms) / n
        ref_res_com[cur_pair] = [cx, cy, cz]

    print(f"  Built reference COM for {len(ref_res_com)} residues")
    return ref_res_com


def shrink_placement_pdb(pdb_path, ref_res_com, skeleton_res):
    """Shrink a single placement PDB in-place.

    Args:
        pdb_path:     path to the placement PDB file
        ref_res_com:  dict of {(idx, code): [cx, cy, cz]} for reference
        skeleton_res: dict to populate with residues that didn't move
    """
    with open(pdb_path, "r") as fh:
        file_lines = fh.readlines()

    # Collect protein residues and ligand atoms
    prot_res = {}       # resnum -> list of ATOM lines
    lig_res_atoms = []  # list of [x, y, z] for HETATM lines

    for line in file_lines:
        if line.startswith("ATOM"):
            resnum = int(line[22:26])
            prot_res.setdefault(resnum, []).append(line)
        elif line.startswith("HETATM"):
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            lig_res_atoms.append([x, y, z])

    res_to_keep = []

    for resnum, atom_lines in prot_res.items():
        # Compute residue COM (heavy atoms only)
        res_com = [0.0, 0.0, 0.0]
        nheavy = 0
        residue_code = ""
        for atom_line in atom_lines:
            if atom_line.strip().endswith("H"):
                continue
            res_com[0] += float(atom_line[30:38])
            res_com[1] += float(atom_line[38:46])
            res_com[2] += float(atom_line[46:54])
            nheavy += 1
            if residue_code == "":
                residue_code = atom_line[17:20]
        if nheavy == 0:
            continue
        res_com = [c / nheavy for c in res_com]

        # Check if residue moved vs reference
        shortest_distance = 100.0
        for (ref_idx, ref_code), ref_com in ref_res_com.items():
            if residue_code == ref_code:
                dist = (
                    (ref_com[0] - res_com[0]) ** 2
                    + (ref_com[1] - res_com[1]) ** 2
                    + (ref_com[2] - res_com[2]) ** 2
                ) ** 0.5
                if dist < shortest_distance:
                    shortest_distance = dist
                if shortest_distance == 0:
                    break

        # Keep if moved
        if shortest_distance > 0.00001:
            res_to_keep.append(resnum)
            continue
        else:
            # Add to skeleton if not already present
            if resnum not in skeleton_res:
                skeleton_res[resnum] = atom_lines

        # Proximity check: any atom within 5 Å of ligand
        for atom_line in atom_lines:
            x = float(atom_line[30:38])
            y = float(atom_line[38:46])
            z = float(atom_line[46:54])
            for lig_atom in lig_res_atoms:
                dist = (
                    (lig_atom[0] - x) ** 2
                    + (lig_atom[1] - y) ** 2
                    + (lig_atom[2] - z) ** 2
                ) ** 0.5
                if dist < 5.0:
                    res_to_keep.append(resnum)
                    break
            if resnum in res_to_keep:
                break

    # Write shrunken PDB
    temp_path = "temp_shrink.pdb"
    with open(temp_path, "w") as out:
        for line in file_lines:
            if not line.startswith("ATOM"):
                out.write(line)
            elif int(line[22:26]) in res_to_keep:
                out.write(line)
    os.replace(temp_path, pdb_path)


def main():
    if len(sys.argv) < 2:
        print(
            "Usage: python shrink_placement_pdbs_to_placement_and_surrounding_residues.py "
            "<reference_pdb> [working_directory]"
        )
        sys.exit(1)

    reference_pdb = sys.argv[1]
    working_location = sys.argv[2] if len(sys.argv) > 2 else os.getcwd()

    if not os.path.isfile(reference_pdb):
        print(f"ERROR: Reference PDB not found: {reference_pdb}")
        sys.exit(1)

    print(f"Reference PDB: {reference_pdb}")
    print(f"Working directory: {working_location}")

    os.chdir(working_location)

    # Look for placements.tar.gz
    if not os.path.isfile("placements.tar.gz"):
        print("ERROR: placements.tar.gz not found in working directory")
        sys.exit(1)

    # Unpack
    print("Unpacking placements.tar.gz ...")
    os.system("tar -xzf placements.tar.gz")
    os.chdir("placements")
    placements_location = os.getcwd()

    # Build reference COM map
    ref_res_com = build_ref_residue_com(reference_pdb)

    # Process each placement PDB
    skeleton_res = {}
    pdb_count = 0

    for root, dirs, files in os.walk(placements_location):
        for f in files:
            if root != placements_location:
                continue
            if not f.endswith(".pdb"):
                continue
            if f == "skeleton.pdb":
                continue

            # Remove minipose files
            if "minipose" in f:
                os.remove(os.path.join(root, f))
                continue

            print(f"  Shrinking: {f}")
            pdb_path = os.path.join(root, f)

            # Remove any matching folder
            folder_name = f.rsplit(".pdb", 1)[0]
            folder_path = os.path.join(placements_location, folder_name)
            if os.path.isdir(folder_path):
                shutil.rmtree(folder_path, ignore_errors=True)

            shrink_placement_pdb(pdb_path, ref_res_com, skeleton_res)
            pdb_count += 1

    # Write skeleton.pdb
    skeleton_path = os.path.join(placements_location, "skeleton.pdb")
    with open(skeleton_path, "w") as skel:
        for resnum in sorted(skeleton_res.keys()):
            skel.writelines(skeleton_res[resnum])
        skel.write("TER\nEND\n")
    print(f"  Wrote skeleton.pdb with {len(skeleton_res)} residues")

    # Recompress
    os.chdir(working_location)
    print("Recompressing placements.tar.gz ...")
    os.system("tar -czf placements.tar.gz placements")
    os.system("rm -drf placements")

    print(f"Done. Shrunk {pdb_count} placement PDB(s).")


if __name__ == "__main__":
    main()
