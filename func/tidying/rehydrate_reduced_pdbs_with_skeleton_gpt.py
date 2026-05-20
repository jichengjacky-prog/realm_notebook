import os, sys

working_location = os.getcwd()
if len(sys.argv) == 2:
    working_location = sys.argv[1]

os.chdir(working_location)
os.system("tar -xzf placements.tar.gz")
os.chdir("placements")
placements_location = os.getcwd()


def resid_key_from_atom_line(line: str):
    """
    Robust PDB fixed-width parsing:
      resname: cols 18-20  (17:20)
      chain : col  22     (21)
      resseq: cols 23-26  (22:26)
      icode : col  27     (26)
    Returns: (chain, resseq_int, icode, resname)
    """
    resname = line[17:20].strip()
    chain = line[21].strip() or " "          # keep blank chain as " "
    resseq_str = line[22:26].strip()
    icode = line[26].strip() or " "          # insertion code
    # resseq should be integer in Rosetta outputs; if not, this will raise (good to know)
    resseq = int(resseq_str)
    return (chain, resseq, icode, resname)


def parse_protein_blocks_fixed(pdb_path):
    """
    Returns:
      blocks: dict[(chain, resseq, icode)] -> list[str ATOM lines]
      prefix_lines: list[str]     # lines before first ATOM
      nonatom_lines: list[str]    # lines after first ATOM that are not protein ATOM
    """
    blocks = {}
    prefix_lines = []
    nonatom_lines = []

    saw_atom = False
    cur_key = None
    cur_lines = []

    with open(pdb_path, "r") as fh:
        for line in fh:
            if line.startswith("ATOM"):
                saw_atom = True
                key = resid_key_from_atom_line(line)[:3]  # (chain, resseq, icode)

                if cur_key is None:
                    cur_key = key
                    cur_lines = [line]
                elif key == cur_key:
                    cur_lines.append(line)
                else:
                    blocks[cur_key] = cur_lines
                    cur_key = key
                    cur_lines = [line]
            else:
                if not saw_atom:
                    prefix_lines.append(line)
                else:
                    nonatom_lines.append(line)

        # finalize last residue at EOF no matter what came after
        if cur_key is not None and cur_lines:
            blocks[cur_key] = cur_lines

    return blocks, prefix_lines, nonatom_lines


# Parse skeleton once
skeleton_blocks, skeleton_prefix, skeleton_nonatom = parse_protein_blocks_fixed("skeleton.pdb")


def sort_reskeys(keys):
    # keys are (chain, resseq, icode)
    # single-chain => chain constant, but this is stable anyway
    return sorted(keys, key=lambda k: (k[0], k[1], k[2]))


for r, d, f in os.walk(placements_location):
    for file in f:
        if r == placements_location and file.endswith(".pdb") and file != "skeleton.pdb":
            print(file)

            placement_blocks, placement_prefix, placement_nonatom = parse_protein_blocks_fixed(file)

            # DEBUG: detect if SER 376 appears in text but isn't parsed into ATOM blocks
            # (This will quickly tell us whether the issue is parsing vs writing.)
            with open(file, "r") as fh:
                text_has_376 = any((" SER " in ln and ln.startswith("ATOM") and ln[22:26].strip() == "376") for ln in fh)

            if text_has_376:
                key_376 = ("A", 376, " ")  # adjust if your chain is blank; see below
                # If your chain is blank in PDB, key should be (" ", 376, " ")
                if key_376 not in placement_blocks and (" ", 376, " ") in placement_blocks:
                    key_376 = (" ", 376, " ")

                if key_376 not in placement_blocks:
                    tail = sort_reskeys(placement_blocks.keys())[-10:]
                    print(f"WARNING: file contains ATOM SER 376 lines, but parser did not capture residue 376 as a block.")
                    print(f"Last 10 parsed residue keys: {tail}")

            # --- FIX: union of residue keys (handles residues present only in placement or only in skeleton) ---
            all_keys_sorted = sort_reskeys(set(skeleton_blocks.keys()) | set(placement_blocks.keys()))

            out_path = "temp.pdb"
            with open(out_path, "w") as out:
                # Keep placement prefix/header (so your remarks/comments at top survive)
                for line in placement_prefix:
                    out.write(line)

                # Protein: placement overrides skeleton
                for key in all_keys_sorted:
                    if key in placement_blocks:
                        out.writelines(placement_blocks[key])
                    else:
                        out.writelines(skeleton_blocks[key])

                # Then write the rest of the placement file (ligand HETATM + footer/comments/END/etc.)
                out.writelines(placement_nonatom)

            os.system(f"mv {out_path} {file}")


os.chdir(working_location)
os.system("tar -czf placements.tar.gz placements")
os.system("rm -drf placements")
