#!/usr/bin/env python3
"""
Shared utility functions for the Rosetta discovery pipeline.

Consolidates reusable logic from:
  - fix_condensed_param_file_spacing.py
  - extract_single_param_from_condensed_file.py
  - score_placed_ligands_with_filtering.py
  - unified_discovery_controller.py
"""

import os
import re
import shutil
import subprocess
import sys
import tempfile

try:
    import CDPL
    import CDPL.Chem as CDPL_Chem
    import CDPL.ConfGen as ConfGen
    CDPKIT_AVAILABLE = True
except Exception:
    CDPKIT_AVAILABLE = False
    CDPL_Chem = None
    ConfGen = None


# ── CDPKit conformer generation (inline, no cross-module imports) ──────────

def generate_cdpkit_conformers(smiles, num_conformers):
    """Generate conformers for a SMILES string using CDPKit.
    
    Args:
        smiles: SMILES string of the ligand
        num_conformers: maximum number of conformers to generate
        
    Returns:
        list of CDPL.Chem.BasicMolecule objects (one per conformer, with 3D coords),
        or None on failure.
    """
    if not CDPKIT_AVAILABLE:
        print("ERROR: CDPKit not available for conformer generation")
        return None
    
    try:
        mol = CDPL_Chem.BasicMolecule()
        CDPL_Chem.parseSMILES(smiles, mol)
        
        ConfGen.prepareForConformerGeneration(mol)
        
        cg = ConfGen.ConformerGenerator()
        cg.settings.maxNumOutputConformers = num_conformers
        cg.settings.minNumOutputConformers = 1
        cg.generate(mol)
        
        num_generated = cg.numConformers
        if num_generated == 0:
            print("WARNING: No conformers generated")
            return None
        
        cg.setConformers(mol)
        
        conformers = []
        for i in range(num_generated):
            conf_mol = CDPL_Chem.BasicMolecule()
            conf_mol.assign(mol)
            CDPL_Chem.applyConformation(conf_mol, i)
            conformers.append(conf_mol)
        
        return conformers
        
    except Exception as e:
        print(f"ERROR in CDPKit conformer generation: {e}")
        return None



def extract_conformer_params(ligand_name, conf_num, chunk, subchunk, 
							  enamine_path, tmp_root=None, output_dir=None):
	"""Extract and fix a single conformer's params from the Enamine library.

	Args:
		output_dir: if given, writes the fixed .params file there and returns its path.
		            If None, returns the params content as a string (no file written).

	Returns:
		str (params text) if output_dir is None, or str (file path) if output_dir is given.
		Returns None on failure.
	"""
	superchunk_str, chunk_str = chunk_to_path(chunk)
	
	# Pipe tar output directly to memory — no temp files on disk
	tar_cmd = (
		f"tar -xzf {enamine_path}/{superchunk_str}/{chunk_str}/"
		f"condensed_params_and_db_{subchunk}.tar.gz "
		f"condensed_params_and_db_{subchunk}/single_conf_params/"
		f"{ligand_name}_shorthand_params.txt "
		f"--strip-components=2 -O"
	)
	result = subprocess.run(tar_cmd, shell=True, capture_output=True, text=True)
	if result.returncode != 0 or not result.stdout:
		print(f"ERROR: Failed to extract params for {ligand_name}_{conf_num}")
		if result.stderr:
			print(result.stderr, file=sys.stderr, end="")
		return None
	
	# All processing is in-memory from here — no intermediate files
	params_text = extract_single_param_text(result.stdout, conf_num)
	if not params_text:
		print(f"ERROR: Failed to extract conformer {conf_num} for {ligand_name}")
		return None

	params_text = fix_params_spacing_text(params_text)
	
	# Write to output_dir if requested, otherwise return text directly
	if output_dir:
		os.makedirs(output_dir, exist_ok=True)
		dest = os.path.join(output_dir, f"{ligand_name}_{conf_num}.params")
		with open(dest, "w") as fh:
			fh.write(params_text)
		return dest
	
	return params_text


def extract_smiles_from_sdf(sdf_file):
	"""Extract SMILES string from first line of SDF file (assumes obabel format)."""
	try:
		# Use obabel to convert SDF to SMILES
		result = subprocess.run(
			f"obabel {sdf_file} -O - -osmi 2>/dev/null | head -1",
			shell=True, capture_output=True, text=True
		)
		if result.stdout:
			smiles = result.stdout.strip().split()[0]  # Get first field
			return smiles
	except Exception as e:
		print(f"ERROR extracting SMILES: {e}")
	return None




def generate_conformers_and_params_cdpkit(ligand_name, smiles, output_dir,
                                          realm_location, num_conformers=150,
                                          tmp_root=None, work_dir=None):
    """Generate conformers from SMILES using CDPKit and create Rosetta params files.

    Uses CDPKit's ConfGen for conformer generation (no external containers needed).
    If work_dir is not given, a temp directory is created under tmp_root.

    Args:
        ligand_name: base name for output files
        smiles: SMILES string of the ligand
        output_dir: directory to write final .params files into
        realm_location: root of the realm discovery tree
        num_conformers: number of conformers to generate (default 150)
        tmp_root: temp directory root (used if work_dir is None)
        work_dir: explicit working directory (if None, one is created)

    Returns:
        list of params file basenames (e.g. ["lig_1.params", "lig_2.params", ...]),
        or empty list on failure.
    """
    if not CDPKIT_AVAILABLE:
        print(f"ERROR: CDPKit pipeline not available for {ligand_name}")
        return []

    _own_work_dir = None
    if work_dir is None:
        work_dir = tempfile.mkdtemp(dir=tmp_root)
        _own_work_dir = work_dir

    try:
        # ── Step A: Generate conformers in-memory with CDPKit ─────────────
        print(f"  Generating {num_conformers} conformers for {ligand_name} with CDPKit...")
        conformers = generate_cdpkit_conformers(smiles, num_conformers)
        if not conformers:
            print(f"ERROR: CDPKit conformer generation failed for {ligand_name}")
            return []
        print(f"  Generated {len(conformers)} conformers")

        # ── Step B: Write each conformer to SDF, then convert to .params ──
        params_files = []
        for i, conf_mol in enumerate(conformers):
            conf_idx = i + 1
            sdf_path = os.path.join(work_dir, f"{ligand_name}_{conf_idx}.sdf")
            write_cdpkit_conformer_to_sdf(conf_mol, sdf_path)

            pf = convert_conformer_to_params(
                sdf_path, ligand_name, conf_idx, work_dir, realm_location
            )
            if not pf:
                print(f"WARNING: molfile_to_params failed for {ligand_name}_{conf_idx}")
                continue

            # Fix spacing so Rosetta can parse correctly
            fix_params_spacing(pf)

            # Copy to output directory
            dest = os.path.join(output_dir, f"{ligand_name}_{conf_idx}.params")
            shutil.copy2(pf, dest)
            params_files.append(f"{ligand_name}_{conf_idx}.params")

        print(f"Generated {len(params_files)} conformer params for {ligand_name}")
        return params_files

    except Exception as e:
        print(f"ERROR generating conformers for {ligand_name}: {e}")
        return []

    finally:
        if _own_work_dir:
            shutil.rmtree(work_dir, ignore_errors=True)


def create_test_params_dir(ligands_list, batch_dir, realm_location, enamine_path, tmp_root=None, round_name="round1"):
	"""Create per-ligand subdirectories under batch_dir/<round_name>/.
	
	Structure:
	    batch_dir/<round_name>/<ligand_name>/
	        test_params/
	            <ligand_name>_<conf_num>.params
	            residue_types.txt       (lists params in this dir)
	            patches.txt
	            exclude_pdb_component_list.txt
	
	Returns the list of ligand subdirs created.
	"""
	lig_dirs = []
	failed_ligands = []
	round_dir = os.path.join(batch_dir, round_name)
	
	for ligand_name, conf_num, chunk, subchunk in ligands_list:
		# Per-ligand directory: batch_dir/<round_name>/<ligand_name>/
		lig_dir = os.path.join(round_dir, ligand_name)
		tp_dir = os.path.join(lig_dir, "test_params")
		os.makedirs(tp_dir, exist_ok=True)
		
		# Create support files in test_params
		run_cmd("touch exclude_pdb_component_list.txt patches.txt", cwd=tp_dir)

		# Skip if params file already exists and is intact (non-empty)
		expected_params = os.path.join(tp_dir, f"{ligand_name}_{conf_num}.params")
		if os.path.isfile(expected_params) and os.path.getsize(expected_params) > 0:
			print(f"SKIP: {expected_params} already exists and intact")
			res_types_file = os.path.join(tp_dir, "residue_types.txt")
			if not os.path.exists(res_types_file):
				with open(res_types_file, "w") as fh:
					fh.write("## atom_type_set and mm-atom_type_set for Rosetta\n")
					fh.write("TYPE_SET_MODE full_atom\n")
					fh.write("ATOM_TYPE_SET fa_standard\n")
					fh.write("ELEMENT_SET default\n")
					fh.write("MM_ATOM_TYPE_SET fa_standard\n")
					fh.write("ORBITAL_TYPE_SET fa_standard\n")
					fh.write("## Params files\n")
			with open(res_types_file, "a") as fh:
				fh.write(f"{ligand_name}_{conf_num}.params\n")
			lig_dirs.append(lig_dir)
			continue
		
		try:
			params_file = extract_conformer_params(
				ligand_name, conf_num, chunk, subchunk, enamine_path,
				batch_dir, output_dir=tp_dir
			)
			
			if params_file:
				# Write per-ligand residue_types.txt in test_params/
				res_types_file = os.path.join(tp_dir, "residue_types.txt")
				if not os.path.exists(res_types_file):
					with open(res_types_file, "w") as fh:
						fh.write("## atom_type_set and mm-atom_type_set for Rosetta\n")
						fh.write("TYPE_SET_MODE full_atom\n")
						fh.write("ATOM_TYPE_SET fa_standard\n")
						fh.write("ELEMENT_SET default\n")
						fh.write("MM_ATOM_TYPE_SET fa_standard\n")
						fh.write("ORBITAL_TYPE_SET fa_standard\n")
						fh.write("## Params files\n")
				
				with open(res_types_file, "a") as fh:
					fh.write(f"{ligand_name}_{conf_num}.params\n")
				
				lig_dirs.append(lig_dir)
			else:
				failed_ligands.append(f"{ligand_name}_{conf_num}")
		
		except Exception as e:
			print(f"ERROR processing {ligand_name}_{conf_num}: {e}")
			failed_ligands.append(f"{ligand_name}_{conf_num}")
	
	if failed_ligands:
		print(f"WARNING: Failed to process {len(failed_ligands)} ligands")
	
	return lig_dirs


def generate_and_add_conformers_to_test_params(batch_dir, ligands_list, 
											   realm_location, enamine_path,
											   num_conformers=150, tmp_root=None, round_name="round2"):
	"""Generate CDPKit conformers and add to per-ligand test_params subdirs.
	
	Conformers are placed in:
	    batch_dir/<round_name>/<ligand_name>/test_params/
	
	Returns True on success, False on failure."""
	
	processed_ligands = set()
	shared_work_dir = tempfile.mkdtemp(dir=tmp_root)
	round_dir = os.path.join(batch_dir, round_name)
	
	any_success = False
	try:
		for ligand_name, conf_num, chunk, subchunk in ligands_list:
			if ligand_name in processed_ligands:
				continue
			processed_ligands.add(ligand_name)
			
			# Per-ligand directory
			lig_dir = os.path.join(round_dir, ligand_name)
			tp_dir = os.path.join(lig_dir, "test_params")
			os.makedirs(tp_dir, exist_ok=True)
			
			print(f"Generating {num_conformers} CDPKit conformers for {ligand_name}...")
			
			# Extract original conformer params (returns text, no file written)
			params_text = extract_conformer_params(
				ligand_name, conf_num, chunk, subchunk, enamine_path, tmp_root=tmp_root
			)
			
			if not params_text:
				print(f"WARNING: Could not extract params for {ligand_name}, skipping conformer generation")
				continue
			
			# Reconstruct SMILES from the .params content
			smiles = extract_smiles_from_params(params_text)
			if not smiles:
				print(f"WARNING: Could not extract SMILES for {ligand_name}, skipping")
				continue
			print(f"  SMILES: {smiles}")
			
			# Generate conformers and params using CDPKit (writes into tp_dir)
			new_params = generate_conformers_and_params_cdpkit(
				ligand_name, smiles, tp_dir, realm_location,
				num_conformers=num_conformers,
				tmp_root=tmp_root, work_dir=shared_work_dir
			)
			
			if new_params:
				# Write per-ligand residue_types.txt
				res_types_file = os.path.join(tp_dir, "residue_types.txt")
				if not os.path.exists(res_types_file):
					with open(res_types_file, "w") as fh:
						fh.write("## atom_type_set and mm-atom_type_set for Rosetta\n")
						fh.write("TYPE_SET_MODE full_atom\n")
						fh.write("ATOM_TYPE_SET fa_standard\n")
						fh.write("ELEMENT_SET default\n")
						fh.write("MM_ATOM_TYPE_SET fa_standard\n")
						fh.write("ORBITAL_TYPE_SET fa_standard\n")
						fh.write("## Params files\n")
					run_cmd("touch exclude_pdb_component_list.txt patches.txt", cwd=tp_dir)
				
				with open(res_types_file, "a") as fh:
					for params_name in new_params:
						fh.write(f"{params_name}\n")
				any_success = True
			else:
				print(f"WARNING: No conformers generated for {ligand_name}")
		
		return any_success
	
	except Exception as e:
		print(f"ERROR in conformer generation: {e}")
		import traceback; traceback.print_exc()
		return False
	
	finally:
		shutil.rmtree(shared_work_dir, ignore_errors=True)


# ---------------------------------------------------------------------------
# SMILES extraction from Rosetta .params files
# ---------------------------------------------------------------------------

def _rosetta_atom_type_to_element(atom_type):
	"""Map a Rosetta atom type string to an element symbol."""
	if not atom_type:
		return "C"
	type_map = [
		("Ntrp", "N"), ("Nhis", "N"), ("Nlys", "N"), ("Narg", "N"),
		("ONH2", "O"), ("Oaro", "O"), ("COO", "C"), ("aroC", "C"),
		("Haro", "H"), ("Hpol", "H"), ("Hapo", "H"),
		("CH", "C"), ("CA", "C"), ("CB", "C"), ("CD", "C"), ("CE", "C"),
		("CG", "C"), ("CZ", "C"), ("NH", "N"),
		("OH", "O"), ("OC", "O"), ("SH", "S"),
		("HA", "H"), ("HB", "H"), ("HG", "H"), ("HD", "H"), ("HE", "H"), ("HH", "H"),
		("C", "C"), ("N", "N"), ("O", "O"), ("S", "S"), ("H", "H"),
		("F", "F"), ("Cl", "Cl"), ("Br", "Br"), ("I", "I"), ("P", "P"),
	]
	for pattern, elem in type_map:
		if atom_type.startswith(pattern):
			return elem
	for ch in atom_type:
		if ch.isalpha():
			return ch.upper()
	return "C"


def extract_smiles_from_params(params_text):
	"""Reconstruct SMILES from Rosetta .params file content.

	Parses ATOM lines for elements and BOND_TYPE lines for connectivity,
	builds the molecular graph with RDKit, and returns a SMILES string.
	Returns None if reconstruction fails.
	"""
	atoms_info = []   # list of (atom_name, element)
	bond_list = []    # list of (atom1, atom2, bond_order)

	for line in params_text.splitlines():
		line = line.strip()
		if not line:
			continue
		if line.startswith("ATOM "):
			parts = line.split()
			if len(parts) >= 4:
				atoms_info.append((parts[1], _rosetta_atom_type_to_element(parts[2])))
		elif line.startswith("BOND_TYPE "):
			parts = line.split()
			if len(parts) >= 4:
				try:
					bond_list.append((parts[1], parts[2], int(parts[3])))
				except ValueError:
					bond_list.append((parts[1], parts[2], 1))

	if not atoms_info or not bond_list:
		print(f"    Params has {len(atoms_info)} atoms, {len(bond_list)} bonds — insufficient")
		return None

	try:
		from rdkit import Chem as rdChem
		mol = rdChem.RWMol()
		atom_map = {}
		for atom_name, elem in atoms_info:
			atom_map[atom_name] = mol.AddAtom(rdChem.Atom(elem))
		for a1, a2, order in bond_list:
			if a1 in atom_map and a2 in atom_map:
				bt = {1: rdChem.BondType.SINGLE, 2: rdChem.BondType.DOUBLE,
					  3: rdChem.BondType.TRIPLE, 4: rdChem.BondType.AROMATIC}.get(order, rdChem.BondType.SINGLE)
				mol.AddBond(atom_map[a1], atom_map[a2], bt)
		mol = mol.GetMol()
		rdChem.SanitizeMol(mol)
		return rdChem.MolToSmiles(mol)
	except Exception as e:
		print(f"    SMILES extraction failed: {e}")
		return None



# ---------------------------------------------------------------------------
# Shell helpers
# ---------------------------------------------------------------------------

def run_cmd(cmd, cwd=None, description="", stream_output=True):
	"""Run a shell command. Returns returncode.

	When stream_output=True (default), stdout/stderr are printed in real time
	so long-running commands don't appear stuck.  Set to False to capture and
	print only after completion.
	"""
	if description:
		print(f"  {description}")

	if stream_output:
		# Stream output in real time — essential for long-running commands
		# so the user can see progress and diagnose hangs.
		with subprocess.Popen(
			cmd, shell=True, cwd=cwd,
			stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
			text=True, bufsize=1
		) as proc:
			for line in proc.stdout:
				print(line, end="")
			proc.wait()
			return proc.returncode
	else:
		result = subprocess.run(cmd, shell=True, capture_output=True, text=True, cwd=cwd)
		if result.stdout:
			print(result.stdout, end="")
		if result.stderr:
			print(result.stderr, file=sys.stderr, end="")
		return result.returncode


# ---------------------------------------------------------------------------
# Enamine library path helpers
# ---------------------------------------------------------------------------

def chunk_to_path(chunk_str):
	"""Convert a 5-digit zero-padded chunk string to (superchunk_str, chunk_str).
	Example: '00123' -> ('1', '00123')"""
	superchunk_str = str(int(chunk_str[:3]))
	return superchunk_str, chunk_str


def write_cdpkit_conformer_to_sdf(conf_mol, sdf_path):
    """Write a single CDPKit conformer molecule to an SDF file."""
    if not CDPKIT_AVAILABLE:
        print("ERROR: CDPKit not available, cannot write SDF")
        return
    writer = CDPL_Chem.FileSDFMolecularGraphWriter(sdf_path)
    CDPL_Chem.setMultiConfExportParameter(writer, False)
    writer.write(conf_mol)
    writer.close()

def fix_params_spacing_text(params_text):
    """Fix column spacing in Rosetta .params content (in-memory).

    Reformat ATOM, BOND_TYPE, CHI, and ICOOR_INTERNAL lines with consistent
    spacing while preserving the original N→C terminus atom order from
    molfile_to_params.py. All other lines pass through unchanged.

    Args:
        params_text: raw .params file content as a string

    Returns:
        fixed .params content as a string
    """
    output_lines = []
    for line in params_text.splitlines(True):
        s = line.rstrip('\n')
        if s.startswith("ATOM"):
            parts = s.split()
            if len(parts) >= 5:
                output_lines.append(
                    "ATOM %-4s %-4s %-4s %.2f\n"
                    % (parts[1], parts[2], parts[3], float(parts[4]))
                )
            else:
                output_lines.append(line)
        elif s.startswith("BOND_TYPE"):
            parts = s.split()
            if len(parts) >= 4:
                output_lines.append(
                    "BOND_TYPE %-4s %-4s %-4s\n"
                    % (parts[1], parts[2], parts[3])
                )
            else:
                output_lines.append(line)
        elif s.startswith("CHI") and not s.startswith("CHI "):
            # PROTON_CHI lines pass through unchanged
            output_lines.append(line)
        elif s.startswith("CHI"):
            parts = s.split()
            if len(parts) >= 6:
                output_lines.append(
                    "CHI %i %-4s %-4s %-4s %-4s\n"
                    % (int(parts[1]), parts[2], parts[3], parts[4], parts[5])
                )
            else:
                output_lines.append(line)
        elif s.startswith("ICOOR_INTERNAL"):
            parts = s.split()
            if len(parts) >= 8:
                output_lines.append(
                    "ICOOR_INTERNAL   %-4s %11.6f %11.6f %11.6f  %-4s  %-4s  %-4s\n"
                    % (parts[1], float(parts[2]), float(parts[3]),
                       float(parts[4]), parts[5], parts[6], parts[7])
                )
            else:
                output_lines.append(line)
        else:
            output_lines.append(line)
    return "".join(output_lines)


def fix_params_spacing(params_file, output_dir=None):
	"""Fix column spacing in a Rosetta .params file (file-based wrapper).
	Reads *params_file*, writes a fixed version, and replaces the original (or
	writes to *output_dir*).

	Returns the path to the fixed file, or None on failure.
	"""
	if not os.path.isfile(params_file):
		print(f"ERROR: params file not found: {params_file}")
		return None

	try:
		with open(params_file, "r") as fh:
			content = fh.read()
		fixed = fix_params_spacing_text(content)

		if output_dir:
			os.makedirs(output_dir, exist_ok=True)
			dest = os.path.join(output_dir, os.path.basename(params_file))
			with open(dest, "w") as fh:
				fh.write(fixed)
			return dest
		else:
			with open(params_file, "w") as fh:
				fh.write(fixed)
			return params_file
	except Exception as e:
		print(f"ERROR fixing params spacing in {params_file}: {e}")
		return None


def convert_conformer_to_params(sdf_path, ligand_name, conf_idx, output_dir,
                                 realm_location):
    """Convert an SDF conformer to Rosetta .params using molfile_to_params.py."""
    params_name = f"{ligand_name}_{conf_idx}"
    sdf_abs = os.path.abspath(sdf_path)
    output_abs = os.path.abspath(output_dir)
    params_file = os.path.join(output_abs, f"{params_name}.params")

    # Run molfile_to_params inside the conformator container (needs rosetta_py)
    conformator_sif = os.path.join(realm_location, "sif", "conformator_container.sif")
    if os.path.exists(conformator_sif):
        sdf_dir = os.path.dirname(sdf_abs)
        params_cmd = (
            f"singularity exec --bind {sdf_dir} {conformator_sif} python "
            f"/conformator_for_container/molfile_to_params.py {sdf_abs} "
            f"-n {params_name} --keep-names --long-names --clobber --no-pdb"
        )
    else:
        params_cmd = (
            f"molfile_to_params.py {sdf_abs} "
            f"-n {params_name} --keep-names --long-names --clobber --no-pdb"
        )

    rc = run_cmd(params_cmd, cwd=output_abs)
    if rc != 0 or not os.path.exists(params_file):
        print(f"WARNING: molfile_to_params failed for {params_name}")
        return None

    # Fix spacing so Rosetta can parse the file correctly
    fix_params_spacing(params_file)

    return params_file

# ---------------------------------------------------------------------------
# Condensed params extraction
# ---------------------------------------------------------------------------

def extract_single_param_text(condensed_text, conf_identifier):
	"""Extract a single conformer's params from condensed text content (in-memory).
	Takes the full text of a condensed params file, returns the .params content
	as a string for the chosen conformer, or None on failure."""
	condense_dict = {}
	in_keys = False
	in_params = False
	output_lines = []

	for line in condensed_text.splitlines(True):
		line_no_newline = line.rstrip("\n")

		# Section headers
		if len(line_no_newline.split()) == 1 and line_no_newline == "KEYS":
			in_keys = True
			continue
		if len(line_no_newline.split()) == 1 and line_no_newline == "PARAMS":
			in_params = True
			in_keys = False
			continue

		# Read key definitions
		if in_keys:
			if not line_no_newline.startswith("_"):
				print(f"Bad key line (missing underscore): {line_no_newline}")
				continue
			if len(line_no_newline.split()) != 2:
				print(f"Bad key line (expected 2 parts): {line_no_newline}")
				continue

			key, entry = line_no_newline.split()
			if key in condense_dict:
				print(f"Duplicate key {key} in dictionary!")
				continue
			condense_dict[key] = entry

		# Read data lines
		if in_params:
			if line_no_newline.count(":") > 1:
				print(f"Line has too many colons: {line_no_newline}")
				continue

			header_side, data_side = line_no_newline.split(":")

			# Check for conserved marker (*)
			conserved = header_side.endswith("*")
			if conserved:
				header_side = header_side[:-1]

			# Translate header
			if header_side in condense_dict:
				header_side = condense_dict[header_side]

			# Pick the correct conformer's data
			split_data = data_side.split(",")
			if conserved:
				my_data = split_data[0]
			else:
				my_data = split_data[conf_identifier - 1]

			# Translate each token
			my_data_split = my_data.split()
			for i in range(len(my_data_split)):
				if my_data_split[i] in condense_dict:
					my_data_split[i] = condense_dict[my_data_split[i]]

			my_data = " " + " ".join(my_data_split)
			output_lines.append(header_side + my_data + "\n")

	return "".join(output_lines)


def extract_single_param(condensed_file_path, conf_identifier, ligand_name, output_dir=None):
	"""Extract a single conformer's params from a condensed params file and write
	it as a standalone .params file readable by Rosetta (file-based wrapper).

	Args:
		condensed_file_path: path to the condensed params text file
		conf_identifier:    1-based conformer index within the condensed file (1-15)
		ligand_name:        base name for the output .params file
		output_dir:         directory to write the .params file (default: cwd)

	Returns:
		path to the written .params file, or None on failure
	"""
	try:
		with open(condensed_file_path, "r") as fh:
			content = fh.read()
	except IOError as e:
		print(f"ERROR: Cannot open file: {e}")
		return None

	params_text = extract_single_param_text(content, conf_identifier)
	if params_text is None:
		return None

	if output_dir:
		os.makedirs(output_dir, exist_ok=True)
		output_path = os.path.join(output_dir, f"{ligand_name}.params")
	else:
		output_path = f"{ligand_name}.params"

	with open(output_path, "w") as fh:
		fh.write(params_text)
	return output_path


# ---------------------------------------------------------------------------
# Placement PDB scoring
# ---------------------------------------------------------------------------

def parse_placement_scores(pdb_path, residue_index_dict=None):
	"""Parse a Rosetta ligand-discovery placement PDB file and extract scoring
	metrics from its header comments.

	Args:
		pdb_path:            path to the placed ligand PDB file
		residue_index_dict:  optional dict mapping original->translated residue indices

	Returns:
		dict with keys: ddg, total_motifs, significant_motifs, real_motif_ratio,
		hbond_motif_count, hbond_motif_energy_sum, found_motif_residues (list),
		or None if the file cannot be read.
	"""
	if not os.path.isfile(pdb_path):
		print(f"ERROR: PDB file not found: {pdb_path}")
		return None

	scores = {
		"ddg": 0.0,
		"total_motifs": 0.0,
		"significant_motifs": 0.0,
		"real_motif_ratio": 0.0,
		"hbond_motif_count": 0,
		"hbond_motif_energy_sum": 0.0,
		"found_motif_residues": [],
	}

	try:
		with open(pdb_path, "r") as fh:
			for line in fh:
				# ddG
				if line.startswith("Scoring: Post-HighResDock system ddG:"):
					scores["ddg"] = float(line.split()[-1].strip())

				# Total motifs
				elif line.startswith("Placement motifs: Total motifs made:"):
					scores["total_motifs"] = float(line.split()[-1].strip())

				# Significant motifs
				elif line.startswith("Placement motifs: Motifs made against significant residues count:"):
					scores["significant_motifs"] = float(line.split()[-1].strip())

				# Real motif ratio
				elif line.startswith("Placement motifs: Real motif ratio:"):
					scores["real_motif_ratio"] = float(line.split()[-1].strip())

				# Per-motif hbond lines
				elif ": Placement motif " in line:
					# Extract residue index from the motif description
					index_match = re.search(r"Hbond_score.*?_(\d{3})[A-Z]?", line)
					if index_match:
						index = index_match.group(1)
						# Translate if a key dict is provided
						if residue_index_dict and index in residue_index_dict:
							index = residue_index_dict[index]
						scores["found_motif_residues"].append(index)

					# Extract hbond score (last colon-separated field)
					hbond_score_str = line.split(":")[-1].strip()
					try:
						hbond_score = float(hbond_score_str)
					except ValueError:
						continue

					scores["hbond_motif_energy_sum"] += hbond_score
					if hbond_score != 0.0:
						scores["hbond_motif_count"] += 1

	except Exception as e:
		print(f"ERROR reading {pdb_path}: {e}")
		return None

	return scores


def compute_weighted_total(scores_dict, score_weights=None):
	"""Apply score weights to a scores dict (from parse_placement_scores) and
	return (weighted_total, weighted_breakdown_dict).

	Default weights are 1.0 for every term.
	"""
	# score = -ddg  (higher = better candidate)
	defaults = {
		"ddg": -1.0,
		"total_motifs": 0.0,
		"significant_motifs": 0.0,
		"real_motif_ratio": 0.0,
		"hbond_motif_count": 0.0,
		"hbond_motif_energy_sum": 0.0,
		"closest_autodock_recovery_rmsd": 0.0,
		"closest_autodock_recovery_ddg": 0.0,
		"strain_energy": 0.0,
	}
	if score_weights:
		defaults.update(score_weights)

	weighted = {}
	total = 0.0
	for term, weight in defaults.items():
		raw = scores_dict.get(term, 0.0)
		w = raw * weight
		weighted[term] = w
		total += w

	return total, weighted


def score_placements(work_dir, weights_file="", min_motif_ratio=0.25):
	"""Score placed ligands by walking a directory, parsing placement PDBs,
	and computing weighted totals.

	If no placement PDB files are found (e.g. because they were already
	tarred by run_rosetta_discovery_search), falls back to reading
	weighted_scores.csv files from subdirectories.

	Args:
		work_dir:        directory to walk for placement PDB files
		weights_file:    optional path to score weights CSV
		min_motif_ratio: minimum real motif ratio filter (default 0.25)

	Returns:
		dict of {ligand_name: total_weighted_score}
	"""
	scores = {}
	weights = load_score_weights(weights_file) if weights_file else {}

	# ── Pass 1: try to find placement PDB files ─────────────────────────
	for root, dirs, files in os.walk(work_dir):
		for file in files:
			if not (file.endswith(".pdb") and "minipose" not in file):
				continue
			pdb_path = os.path.join(root, file)
			parsed = parse_placement_scores(pdb_path)
			if not parsed:
				continue
			if parsed["real_motif_ratio"] < min_motif_ratio:
				continue
			total, _ = compute_weighted_total(parsed, weights)
			lig_name = file.replace(".pdb", "")
			scores[lig_name] = total

	# ── Pass 2: fallback — aggregate weighted_scores.csv files ──────────
	# (generated by score_placed_ligands_with_filtering.py during Rosetta
	#  post-processing; placement PDBs may have been tarred & deleted)
	if not scores:
		import csv as _csv
		for root, dirs, files in os.walk(work_dir):
			for file in files:
				if file != "weighted_scores.csv":
					continue
				csv_path = os.path.join(root, file)
				try:
					with open(csv_path, "r") as fh:
						reader = _csv.DictReader(fh)
						for row in reader:
							# Ligand name from the PDB filename in the first column
							pdb_name = os.path.basename(row.get("file", ""))
							lig_name = pdb_name.replace(".pdb", "")
							if not lig_name:
								continue
							# Use the weighted total (last column)
							total_str = row.get("total", "0")
							try:
								total = float(total_str)
							except ValueError:
								continue
							scores[lig_name] = total
				except Exception as e:
					print(f"WARNING: could not read {csv_path}: {e}")

	return scores


def load_score_weights(weights_csv_path):
	"""Load score weights from a CSV file (term,weight per line, no header).
	Returns a dict of {term: weight}."""
	weights = {}
	if not weights_csv_path or not os.path.isfile(weights_csv_path):
		return weights

	with open(weights_csv_path, "r") as fh:
		for line in fh:
			line = line.strip()
			if not line:
				continue
			parts = line.split(",")
			if len(parts) >= 2:
				try:
					weights[parts[0].strip()] = float(parts[1].strip())
				except ValueError:
					print(f"WARNING: bad weight line: {line}")
	return weights


def load_residue_index_key(key_file_path):
	"""Load a residue-index translation key CSV.
	Format: res_type,original_index,translated_index[,difference]
	Returns dict mapping original_index -> translated_index."""
	key_dict = {}
	if not key_file_path or not os.path.isfile(key_file_path):
		return key_dict

	with open(key_file_path, "r") as fh:
		for line in fh:
			if line.startswith("res_type"):
				continue
			parts = line.strip().split(",")
			if len(parts) >= 3:
				key_dict[parts[1].strip()] = parts[2].strip()
	return key_dict


# ---------------------------------------------------------------------------
# Rosetta ligand discovery search
# ---------------------------------------------------------------------------

def run_rosetta_discovery_search(target_pdb, anchor_residue_string, motifs_file,
								  test_params_dir,discovery_root, atr, rep, ddg,
								  extra_args_file=None, work_dir=None):
	"""Run a single Rosetta ligand discovery search.

	Writes a Rosetta args file, executes the search via singularity, and
	post-processes the output (placements directory, scoring, compression).

	Args:
		target_pdb:             path to target PDB file
		anchor_residue_string:  comma-separated Rosetta-indexed anchor residues
		motifs_file:            path to motifs file
		test_params_dir:        path to test_params directory (must be named "test_params")
		discovery_root:         root directory of the discovery pipeline
		atr, rep, ddg:          fa_atr, fa_rep, ddg score cutoffs
		extra_args_file:        optional file with additional Rosetta args
		work_dir:               working directory (default: current directory)

	Returns:
		True on success, False on failure.
	"""
	# Resolve all paths to absolute BEFORE changing working directory
	target_pdb       = os.path.abspath(target_pdb)
	motifs_file      = os.path.abspath(motifs_file)
	test_params_dir  = os.path.abspath(test_params_dir)
	discovery_root   = os.path.abspath(discovery_root)
	if extra_args_file:
		extra_args_file = os.path.abspath(extra_args_file)

	if work_dir:
		os.makedirs(work_dir, exist_ok=True)
		orig_cwd = os.getcwd()
		os.chdir(work_dir)
	else:
		orig_cwd = None

	try:
		# Ensure test_params_dir ends with /
		if not test_params_dir.endswith("/"):
			test_params_dir = test_params_dir + "/"

		# Write Rosetta args file to the current working directory
		# (so the singularity --bind {cwd}/args:/input/args always finds it)
		with open(os.path.join(os.getcwd(), "args"), "w") as args_fh:
			args_fh.write("#keep seed constant\n")
			args_fh.write("-constant_seed 1\n")
			args_fh.write("#ignore unrecognized residues to help mitigate crashes\n")
			args_fh.write("-ignore_unrecognized_res\n")
			args_fh.write("#handle ligand repeats if using multiple anchor residues\n")
			args_fh.write("-in::file::override_database_params true\n")
			args_fh.write("#constrain coordinates\n")
			args_fh.write("-constrain_relax_to_start_coords\n")
			args_fh.write("#keep all placements\n")
			args_fh.write("-best_pdbs_to_keep 0\n")
			
			# User-input dependent
			args_fh.write("#mapped protein system\n")
			args_fh.write("-s /input/" + os.path.basename(target_pdb) + "\n")
			args_fh.write("#mapped motifs file\n")
			args_fh.write("-motif_filename /input/" + os.path.basename(motifs_file) + "\n")
			args_fh.write("#mapped test_params directory\n")
			args_fh.write("-params_directory_path /input/" + os.path.basename(test_params_dir.rstrip("/")) + "/\n")
			args_fh.write("#rosetta-indexed anchor residue index/indices\n")
			args_fh.write("-protein_discovery_locus " + anchor_residue_string + "\n")
			args_fh.write("#fa_atr cutoff\n")
			args_fh.write("-fa_atr_cutoff " + str(atr) + "\n")
			args_fh.write("#fa_rep cutoff\n")
			args_fh.write("-fa_rep_cutoff " + str(rep) + "\n")
			args_fh.write("#ddg cutoff\n")
			args_fh.write("-ddg_cutoff " + str(ddg) + "\n")

			# Extra user args
			if extra_args_file and os.path.isfile(extra_args_file):
				args_fh.write("###################################################\n")
				args_fh.write("#extra user args from: " + extra_args_file + "\n")
				with open(extra_args_file, "r") as extra_fh:
					args_fh.write(extra_fh.read())

		# Build singularity command
		rosetta_sif = os.path.join(discovery_root, "sif", "rosetta_condensed_6_25_2024.sif")
		rosetta_cmd = " ".join([
			"singularity exec",
			"--bind " + test_params_dir + ":/input/test_params/",
			"--bind " + os.getcwd() + "/args:/input/args",
			"--bind " + target_pdb + ":/input/" + os.path.basename(target_pdb),
			"--bind " + motifs_file + ":/input/" + os.path.basename(motifs_file),
			rosetta_sif,
			"/rosetta/source/bin/ligand_discovery_search_protocol.linuxgccrelease @/input/args",
		])
		##debug print the command
		print(rosetta_cmd)

		print("Running Rosetta discovery search...")
		rosetta_result  = run_cmd(rosetta_cmd)

		if rosetta_result != 0:
			print("ERROR: Rosetta discovery search failed")
			return False

		# Post-processing: organize placements
		run_cmd("mkdir -p placements")
		run_cmd("mv *pdb placements 2>/dev/null; true")

		os.chdir("placements")
		# Rename each PDB with the anchor residue prefix
		for r, d, f in os.walk(os.getcwd()):
			for file in f:
				if file.endswith(".pdb") and r == os.getcwd():
					os.rename(file, "res" + anchor_residue_string + "_" + file)

		# Run placement scoring
		score_script = os.path.join(discovery_root, "function", "rosetta",
									"score_placed_ligands_with_filtering.py")
		run_cmd("python " + score_script)

		# Copy score CSVs up one level
		run_cmd("cp *csv .. 2>/dev/null; true")

		os.chdir("..")

		# Compress placements
		run_cmd("tar -czf placements.tar.gz placements")
		run_cmd("rm -drf placements")

		# Dehydrate to minimize overhead
		dehydrate_script = os.path.join(discovery_root, "function", "tidying",
										"shrink_placement_pdbs_to_placement_and_surrounding_residues.py")
		run_cmd("python " + dehydrate_script + " " + target_pdb)

		print("Rosetta discovery search complete.")
		return True

	except Exception as e:
		print(f"ERROR in Rosetta discovery search: {e}")
		return False

	finally:
		if orig_cwd:
			os.chdir(orig_cwd)


# ═══════════════════════════════════════════════════════════════════════════════
# Tidying / cleanup utilities
# ═══════════════════════════════════════════════════════════════════════════════

def condense_test_params_directories(root_dir, throttle_limit=300, queue="long"):
    """Walk root_dir and submit bsub jobs to tar+gzip then delete every test_params directory found.

    Each job:  tar -czf test_params.tar.gz test_params && rm -drf test_params

    Args:
        root_dir:         top-level directory to walk
        throttle_limit:   max concurrent bsub jobs before sleeping (default 300)
        queue:            LSF queue name (default "long")
    """
    import time
    root_dir = os.path.abspath(root_dir)
    orig_cwd = os.getcwd()
    count = 0

    for r, d, _ in os.walk(root_dir):
        for dire in d:
            if dire != "test_params":
                continue
            count += 1
            print(count)

            os.chdir(r)
            cmd = (
                f'bsub -q {queue} -W 1:00 -u "" -R "rusage[mem=5000]" '
                f'"tar -czf test_params.tar.gz test_params && rm -drf test_params"'
            )
            print(cmd)
            os.system(cmd)

            # throttle
            job_count = int(subprocess.run(
                "bjobs | wc -l", shell=True, capture_output=True, text=True
            ).stdout.strip() or 0)
            while job_count > throttle_limit:
                time.sleep(1)
                job_count = int(subprocess.run(
                    "bjobs | wc -l", shell=True, capture_output=True, text=True
                ).stdout.strip() or 0)

    os.chdir(orig_cwd)
    print(f"Submitted {count} condense jobs.")


def copy_placements_from_list(list_file, target_destination, throttle_limit=500):
    """Copy placement PDBs listed in a CSV file to a review directory.

    Each listed placement is extracted from its parent placements.tar.gz,
    rehydrated, and moved to the target destination.  Submitted as individual
    bsub jobs (short queue).

    Args:
        list_file:          path to CSV (first column = full path to placement PDB)
        target_destination: directory to copy extracted placements into
        throttle_limit:     max concurrent bsub short-queue jobs (default 500)
    """
    import shlex, time

    target_destination = target_destination.rstrip("/") + "/"
    os.makedirs(target_destination, exist_ok=True)
    os.chdir(target_destination)

    file_counter = 0
    dir_counter = 0
    os.makedirs(str(dir_counter), exist_ok=True)
    os.chdir(str(dir_counter))

    rehydrate_script = "/pi/summer.thyme-umw/enamine-REAL-2.6billion/umass_chan_REAL-M_platform/tidying/rehydrate_reduced_pdbs_with_skeleton_gpt.py"

    with open(list_file, "r") as fh:
        for line in fh:
            if line.startswith("file,ddg,"):
                continue
            full_file = line.split(",")[0]
            placements_file = full_file.split("/placement/")[0] + "/placements.tar.gz"
            file_root = full_file.split("/placement/")[1].split(".pdb")[0]

            file_counter += 1
            if file_counter % 100 == 0:
                dir_counter += 1
                os.chdir(target_destination)
                os.makedirs(str(dir_counter), exist_ok=True)
                os.chdir(str(dir_counter))

            os.makedirs(file_root, exist_ok=True)
            os.chdir(file_root)

            job_cmd = (
                f"cp {shlex.quote(placements_file)} . && "
                f"python {rehydrate_script} && "
                f"tar -xzf placements.tar.gz --strip-components=1 "
                f"{shlex.quote('placements/' + file_root + '.pdb')} && "
                f"mv {shlex.quote(file_root + '.pdb')} .. && "
                f"rm -rf placements* && "
                f"rm -drf ../{file_root}"
            )

            cmd = [
                "bsub", "-q", "short", "-W", "2:00", "-u", "",
                "-R", "rusage[mem=5000]",
                "bash", "-lc", job_cmd,
            ]

            # throttle
            job_count = int(subprocess.run(
                "bjobs | grep short | wc -l", shell=True, capture_output=True, text=True
            ).stdout.strip() or 0)
            while job_count > throttle_limit:
                time.sleep(1)
                job_count = int(subprocess.run(
                    "bjobs | grep short | wc -l", shell=True, capture_output=True, text=True
                ).stdout.strip() or 0)

            print(" ".join(cmd))
            subprocess.run(cmd, check=True)
            os.chdir("..")


def dehydrate_placements(cleaner_root, reference_pdb, throttle_limit=100):
    """Submit bsub dehydration jobs for every directory containing both
    placements.tar.gz and raw_scores.csv.

    Args:
        cleaner_root:   top-level directory to search
        reference_pdb:  path to reference PDB for skeleton generation
        throttle_limit: max concurrent short-queue bsub jobs (default 100)
    """
    import time
    shrink_script = "/pi/summer.thyme-umw/enamine-REAL-2.6billion/umass_chan_REAL-M_platform/tidying/shrink_placement_pdbs_to_placement_and_surrounding_residues.py"

    os.chdir(cleaner_root)
    for r, _, files in os.walk(cleaner_root):
        for f in files:
            if f != "placements.tar.gz":
                continue
            if not os.path.exists(os.path.join(r, "raw_scores.csv")):
                continue

            print(r + "/raw_scores.csv " + r + "/placements.tar.gz")
            cmd = (
                f'bsub -q short -W 1:00 -u "" -R "rusage[mem=5000]" '
                f'"python {shrink_script} {reference_pdb} {r}"'
            )
            print(cmd)
            os.system(cmd)

            job_count = int(subprocess.run(
                "bjobs | grep short | wc -l", shell=True, capture_output=True, text=True
            ).stdout.strip() or 0)
            while job_count > throttle_limit:
                time.sleep(1)
                job_count = int(subprocess.run(
                    "bjobs | grep short | wc -l", shell=True, capture_output=True, text=True
                ).stdout.strip() or 0)


def rehydrate_placements(working_location=None):
    """Rehydrate reduced placement PDBs in-place using a skeleton.pdb.

    Expects placements.tar.gz in the working directory.  Unpacks it,
    replaces missing protein residues in each placement PDB with the
    corresponding atoms from skeleton.pdb, then recompresses.

    Args:
        working_location: directory containing placements.tar.gz and
                          skeleton.pdb (default: os.getcwd())
    """
    if working_location is None:
        working_location = os.getcwd()

    os.chdir(working_location)
    os.system("tar -xzf placements.tar.gz")
    os.chdir("placements")
    placements_location = os.getcwd()

    # ── helpers ─────────────────────────────────────────────────────────
    def _resid_key_from_atom_line(line):
        resname = line[17:20].strip()
        chain = line[21].strip() or " "
        resseq = int(line[22:26].strip())
        icode = line[26].strip() or " "
        return (chain, resseq, icode, resname)

    def _parse_protein_blocks(pdb_path):
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
                    key = _resid_key_from_atom_line(line)[:3]
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
            if cur_key is not None and cur_lines:
                blocks[cur_key] = cur_lines
        return blocks, prefix_lines, nonatom_lines

    def _sort_reskeys(keys):
        return sorted(keys, key=lambda k: (k[0], k[1], k[2]))

    # ── main ────────────────────────────────────────────────────────────
    skeleton_blocks, skeleton_prefix, skeleton_nonatom = _parse_protein_blocks("skeleton.pdb")

    for r, d, f in os.walk(placements_location):
        for file in f:
            if r != placements_location or not file.endswith(".pdb") or file == "skeleton.pdb":
                continue
            print(file)

            placement_blocks, placement_prefix, placement_nonatom = _parse_protein_blocks(file)
            all_keys_sorted = _sort_reskeys(set(skeleton_blocks.keys()) | set(placement_blocks.keys()))

            out_path = "temp.pdb"
            with open(out_path, "w") as out:
                for line in placement_prefix:
                    out.write(line)
                for key in all_keys_sorted:
                    if key in placement_blocks:
                        out.writelines(placement_blocks[key])
                    else:
                        out.writelines(skeleton_blocks[key])
                out.writelines(placement_nonatom)

            os.system(f"mv {out_path} {file}")

    os.chdir(working_location)
    os.system("tar -czf placements.tar.gz placements")
    os.system("rm -drf placements")


def shrink_placements(reference_pdb, working_location=None):
    """Shrink placement PDBs to only residues within 5 Å of the ligand plus
    residues that moved relative to the reference.

    Operates on placements.tar.gz in working_location: unpacks, shrinks
    every placement PDB, writes skeleton.pdb of unmoved residues, and
    recompresses.

    Args:
        reference_pdb:    path to reference (pre-docked) PDB
        working_location: directory containing placements.tar.gz
                          (default: os.getcwd())
    """
    if working_location is None:
        working_location = os.getcwd()

    os.chdir(working_location)
    os.system("tar -xzf placements.tar.gz")
    os.chdir("placements")
    placements_location = os.getcwd()

    # ── 1. Build reference residue centres of mass ──────────────────────
    ref_res_com = {}
    cur_pair = ("", "")
    ref_res_atoms = []

    with open(reference_pdb, "r") as ref_file:
        for line in ref_file:
            if not line.startswith("ATOM"):
                if line.startswith("TER") and ref_res_atoms and cur_pair[0] != "":
                    x_sum = sum(a[0] for a in ref_res_atoms)
                    y_sum = sum(a[1] for a in ref_res_atoms)
                    z_sum = sum(a[2] for a in ref_res_atoms)
                    n = len(ref_res_atoms)
                    ref_res_com[cur_pair] = [x_sum / n, y_sum / n, z_sum / n]
                    ref_res_atoms = []
                continue

            if line.strip().endswith("H"):
                continue

            xyz = [float(line[30:38]), float(line[38:46]), float(line[46:54])]
            res_index = int(line[22:26])
            res_code = line[17:20]

            if cur_pair[0] == "":
                cur_pair = (res_index, res_code)

            if cur_pair[0] == res_index:
                ref_res_atoms.append(xyz)
            else:
                if ref_res_atoms:
                    x_sum = sum(a[0] for a in ref_res_atoms)
                    y_sum = sum(a[1] for a in ref_res_atoms)
                    z_sum = sum(a[2] for a in ref_res_atoms)
                    n = len(ref_res_atoms)
                    ref_res_com[cur_pair] = [x_sum / n, y_sum / n, z_sum / n]
                cur_pair = (res_index, res_code)
                ref_res_atoms = [xyz]

    # ── 2. Process each placement PDB ───────────────────────────────────
    skeleton_res = {}

    for r, d, f in os.walk(placements_location):
        for file in f:
            if r != placements_location or not file.endswith(".pdb"):
                continue
            print(file)

            if "minipose" in file:
                os.remove(os.path.join(r, file))
                continue
            if file == "skeleton.pdb":
                continue

            # delete matching folder
            file_base = file.rsplit(".pdb", 1)[0]
            folder_path = os.path.join(placements_location, file_base)
            if os.path.isdir(folder_path):
                shutil.rmtree(folder_path, ignore_errors=True)

            # buffer file
            with open(file, "r") as fh:
                file_lines = fh.readlines()

            prot_res = {}
            lig_res_atoms = []

            for line in file_lines:
                if line.startswith("ATOM"):
                    resnum = int(line[22:26])
                    prot_res.setdefault(resnum, []).append(line)
                elif line.startswith("HETATM"):
                    lig_res_atoms.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])

            res_to_keep = []

            for res in prot_res:
                # compute residue COM (heavy atoms only)
                res_com = [0.0, 0.0, 0.0]
                nheavy = 0
                residue_code = ""
                for res_atom in prot_res[res]:
                    if res_atom.strip().endswith("H"):
                        continue
                    res_com[0] += float(res_atom[30:38])
                    res_com[1] += float(res_atom[38:46])
                    res_com[2] += float(res_atom[46:54])
                    nheavy += 1
                    if residue_code == "":
                        residue_code = res_atom[17:20]
                if nheavy == 0:
                    continue
                res_com = [c / nheavy for c in res_com]

                # check if residue moved vs reference
                shortest_distance = 100.0
                for (ref_idx, ref_code), ref_com in ref_res_com.items():
                    if residue_code == ref_code:
                        dist = ((ref_com[0] - res_com[0]) ** 2 +
                                (ref_com[1] - res_com[1]) ** 2 +
                                (ref_com[2] - res_com[2]) ** 2) ** 0.5
                        if dist < shortest_distance:
                            shortest_distance = dist
                        if shortest_distance == 0:
                            break

                if shortest_distance > 0.00001:
                    res_to_keep.append(res)
                    continue
                else:
                    # add to skeleton if not already present
                    if res not in skeleton_res:
                        skeleton_res[res] = prot_res[res]

                # proximity check: any atom within 5 Å of ligand
                for res_atom in prot_res[res]:
                    x = float(res_atom[30:38])
                    y = float(res_atom[38:46])
                    z = float(res_atom[46:54])
                    for lig_atom in lig_res_atoms:
                        dist = ((lig_atom[0] - x) ** 2 +
                                (lig_atom[1] - y) ** 2 +
                                (lig_atom[2] - z) ** 2) ** 0.5
                        if dist < 5:
                            res_to_keep.append(res)
                            break
                    if res in res_to_keep:
                        break

            # write shrunken PDB
            with open("temp.pdb", "w") as out:
                for line in file_lines:
                    if not line.startswith("ATOM"):
                        out.write(line)
                    elif int(line[22:26]) in res_to_keep:
                        out.write(line)
            os.system("mv temp.pdb " + file)

    # ── 3. Write skeleton ──────────────────────────────────────────────
    with open("skeleton.pdb", "w") as skel:
        for res in sorted(skeleton_res.keys()):
            skel.writelines(skeleton_res[res])
        skel.write("TER\nEND\n")

    # ── 4. Recompress ──────────────────────────────────────────────────
    os.chdir(working_location)
    os.system("tar -czf placements.tar.gz placements")
    os.system("rm -drf placements")


if __name__ == "__main__":
	import argparse

	parser = argparse.ArgumentParser(description="Discovery pipeline utilities")
	sub = parser.add_subparsers(dest="command")

	# fix-params
	p_fix = sub.add_parser("fix-params", help="Fix spacing in a condensed params file")
	p_fix.add_argument("params_file", help="Path to the .params file")

	# extract-param
	p_ext = sub.add_parser("extract-param", help="Extract single conformer from condensed file")
	p_ext.add_argument("condensed_file", help="Path to condensed params file")
	p_ext.add_argument("conf_id", type=int, help="Conformer identifier (1-based)")
	p_ext.add_argument("ligand_name", help="Base ligand name")
	p_ext.add_argument("--output-dir", "-o", default=None, help="Output directory")

	# parse-placement
	p_place = sub.add_parser("parse-placement", help="Parse scores from a placement PDB")
	p_place.add_argument("pdb_path", help="Path to the placed ligand PDB")
	p_place.add_argument("--residue-key", "-k", default=None, help="Residue index translation key CSV")

	args = parser.parse_args()

	if args.command == "fix-params":
		result = fix_params_spacing(args.params_file)
		print(f"Fixed: {result}" if result else "Failed.")

	elif args.command == "extract-param":
		result = extract_single_param(args.condensed_file, args.conf_id,
									  args.ligand_name, args.output_dir)
		print(f"Written: {result}" if result else "Failed.")

	elif args.command == "parse-placement":
		key_dict = load_residue_index_key(args.residue_key) if args.residue_key else None
		scores = parse_placement_scores(args.pdb_path, key_dict)
		if scores:
			total, weighted = compute_weighted_total(scores)
			print(f"ddg={scores['ddg']} motifs={scores['total_motifs']} "
				  f"sig_motifs={scores['significant_motifs']} ratio={scores['real_motif_ratio']} "
				  f"hbond_count={scores['hbond_motif_count']} hbond_sum={scores['hbond_motif_energy_sum']}")
			print(f"Found motif residues: {scores['found_motif_residues']}")
		else:
			print("Failed to parse placement PDB.")
			sys.exit(1)

	else:
		parser.print_help()
