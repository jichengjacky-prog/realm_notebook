#!/usr/bin/env python3
"""
CDPKit Conformer Pipeline

Pipeline to generate 150 conformers per ligand using CDPKit (instead of the
conformator container), then rank each conformer with Rosetta energy scores.

Flow:
  1. Read shapedb input list  (score, ligand_conf, chunk, subchunk)
  2. Extract .params from Enamine library
  3. Extract/reconstruct SMILES for each ligand
  4. Generate 150 conformers with CDPKit
  5. Convert each conformer to Rosetta .params format
  6. Run Rosetta scoring on each conformer
  7. Rank conformers and output top hits

Usage:
    python cdpkit_conformer_pipeline.py \
        --input-list test/conformer_input.text \
        --target-pdb input_pdb/processed/target.pdb \
        --anchor-residues 79 \
        --motifs-file motifs/FINAL_motifs_list_filtered_2_3_2023.motifs \
        --output-dir output/cdpkit_run \
        [--num-conformers 150] \
        [--top-n 10]
"""

import argparse
import heapq
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
# ── CDPKit imports ──────────────────────────────────────────────────────────
try:
    import CDPL
    import CDPL.Chem as Chem
    import CDPL.ConfGen as ConfGen
    CDPKIT_AVAILABLE = True
    print("CDPKit is available.")
except Exception:
    # ImportError or linter/static analysis resolution failures
    CDPKIT_AVAILABLE = False
    Chem = None
    ConfGen = None
    print("WARNING: CDPKit not found. Conformer generation will be skipped.",
          file=sys.stderr)

# ── Rosetta utils (from sibling discovery modules) ──────────────────────────
# Add parent dir to path so we can import utils
_this_dir = os.path.dirname(os.path.abspath(__file__))
if _this_dir not in sys.path:
    sys.path.insert(0, _this_dir)

try:
    from utils import (
        run_cmd, chunk_to_path, extract_single_param_text,
        fix_params_spacing_text, parse_placement_scores,
        compute_weighted_total, load_score_weights,
    )
except ImportError:
    print("WARNING: Cannot import discovery/utils.py. Scoring disabled.",
          file=sys.stderr)
    # Minimal fallbacks
    def run_cmd(cmd, cwd=None, description=""):
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True, cwd=cwd)
        if result.stdout:
            print(result.stdout, end="")
        if result.stderr:
            print(result.stderr, file=sys.stderr, end="")
        return result.returncode

    def chunk_to_path(chunk_str):
        chunk_str = str(chunk_str).zfill(5)
        superchunk_str = str(int(chunk_str[:3]))
        return superchunk_str, chunk_str

    parse_placement_scores = None
    compute_weighted_total = None
    load_score_weights = None

    def extract_single_param_text(condensed_text, conf_identifier):
        # Minimal extraction — see utils.py for full implementation
        lines = []
        in_block = False
        for line in condensed_text.splitlines(True):
            if line.startswith(f"CONFORMER {conf_identifier}"):
                in_block = True
                continue
            if in_block:
                if line.startswith("CONFORMER") and not line.startswith(f"CONFORMER {conf_identifier}"):
                    break
                lines.append(line)
        return "".join(lines) if lines else None

    def fix_params_spacing_text(params_text):
        return params_text  # no-op fallback


# ═══════════════════════════════════════════════════════════════════════════════
# Step 1 — Extract params from Enamine library
# ═══════════════════════════════════════════════════════════════════════════════

def extract_params_from_enamine(ligand_name, conf_num, chunk, subchunk,
                                 enamine_path, output_dir):
    """Extract and fix a single conformer's .params from the Enamine library.
    Returns path to the written .params file, or None on failure."""
    superchunk_str, chunk_str = chunk_to_path(chunk)

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
        return None

    params_text = extract_single_param_text(result.stdout, conf_num)
    if not params_text:
        print(f"ERROR: Failed to extract conformer {conf_num} for {ligand_name}")
        return None

    params_text = fix_params_spacing_text(params_text)

    os.makedirs(output_dir, exist_ok=True)
    dest = os.path.join(output_dir, f"{ligand_name}_{conf_num}.params")
    with open(dest, "w") as fh:
        fh.write(params_text)
    return dest


# ═══════════════════════════════════════════════════════════════════════════════
# Step 2 — Extract SMILES for CDPKit input
# ═══════════════════════════════════════════════════════════════════════════════

def _rosetta_atom_type_to_element(atom_type):
    """Map a Rosetta atom type string to an element symbol.
    Examples: CH3→C, OH1→O, Ntrp→N, aroC→C, Haro→H, COO→C, ONH2→O"""
    if not atom_type:
        return "C"
    # Map specific Rosetta atom types to elements
    # Format: (pattern, element) — checked in order
    type_map = [
        # Multi-char exact/prefix matches
        ("Ntrp", "N"), ("Nhis", "N"), ("Nlys", "N"), ("Narg", "N"),
        ("ONH2", "O"), ("Oaro", "O"),
        ("COO", "C"),
        ("aroC", "C"),
        ("Haro", "H"), ("Hpol", "H"), ("Hapo", "H"),
        # Standard atom type prefixes
        ("CH", "C"), ("CA", "C"), ("CB", "C"), ("CD", "C"), ("CE", "C"),
        ("CG", "C"), ("CZ", "C"),
        ("NH", "N"),
        ("OH", "O"), ("OC", "O"),
        ("SH", "S"),
        ("HA", "H"), ("HB", "H"), ("HG", "H"), ("HD", "H"), ("HE", "H"),
        ("HH", "H"),
        # Single-letter elements
        ("C", "C"), ("N", "N"), ("O", "O"), ("S", "S"), ("H", "H"),
        ("F", "F"), ("Cl", "Cl"), ("Br", "Br"), ("I", "I"), ("P", "P"),
    ]
    for pattern, elem in type_map:
        if atom_type.startswith(pattern):
            return elem
    # Fallback: first alphabetic character
    for ch in atom_type:
        if ch.isalpha():
            return ch.upper()
    return "C"


def _bond_order_to_rdkit(order_int):
    """Convert Rosetta BOND_TYPE integer to RDKit bond type."""
    try:
        from rdkit import Chem as rdChem
    
        mapping = {1: rdChem.BondType.SINGLE, 2: rdChem.BondType.DOUBLE,
                   3: rdChem.BondType.TRIPLE, 4: rdChem.BondType.AROMATIC}
        return mapping.get(int(order_int), rdChem.BondType.SINGLE)
    except ImportError:
        return int(order_int)  # return raw int for CDPKit fallback


def _build_mol_from_params(atoms_info, bond_list):
    """Build an RDKit molecule from atom elements and (a1, a2, bond_order) bonds.
    Returns SMILES string or None."""
    try:
        from rdkit import Chem as rdChem
        mol = rdChem.RWMol()
        atom_map = {}
        for atom_name, elem in atoms_info:
            atom_idx = mol.AddAtom(rdChem.Atom(elem))
            atom_map[atom_name] = atom_idx
        for a1, a2, order in bond_list:
            if a1 in atom_map and a2 in atom_map:
                bt = _bond_order_to_rdkit(order)
                mol.AddBond(atom_map[a1], atom_map[a2], bt)
        mol = mol.GetMol()
        rdChem.SanitizeMol(mol)
        return rdChem.MolToSmiles(mol)
    except Exception as e:
        print(f"    RDKit molecule build failed: {e}")
        return None


def _build_mol_cdpkit(atoms_info, bond_list):
    """Build a CDPKit molecule and return SMILES, or None."""
    try:
        import CDPL.Chem as CDPL_Chem
        mol = CDPL_Chem.BasicMolecule()
        atom_map = {}
        for atom_name, elem in atoms_info:
            atom = mol.addAtom()
            CDPL_Chem.setSymbol(atom, elem)
            atom_map[atom_name] = atom.getIndex()
        for a1, a2, order in bond_list:
            if a1 in atom_map and a2 in atom_map:
                bt = int(order)
                mol.addBond(atom_map[a1], atom_map[a2], bt)
        return CDPL_Chem.generateSMILES(mol)
    except Exception as e:
        print(f"    CDPKit molecule build failed: {e}")
        return None


def _build_smiles_pure_python(atoms_info, bond_list):
    """Pure-Python SMILES generator from atom elements and bonds.
    No external dependencies. Handles single/double/triple/aromatic bonds.

    Uses a DFS traversal with ring-closure numbering to produce a
    syntactically valid (though not strictly canonical) SMILES string.
    """
    n = len(atoms_info)
    if n == 0:
        return None

    elements = [elem for _, elem in atoms_info]
    name_to_idx = {name: i for i, (name, _) in enumerate(atoms_info)}

    # Build adjacency list: adj[i] = list of (j, bond_order)
    adj = [[] for _ in range(n)]
    for a1, a2, order in bond_list:
        i, j = name_to_idx.get(a1), name_to_idx.get(a2)
        if i is not None and j is not None and i != j:
            order_int = int(order)
            adj[i].append((j, order_int))
            adj[j].append((i, order_int))

    if not any(adj):
        return None

    # Bond order → SMILES symbol
    bond_sym = {1: '', 2: '=', 3: '#', 4: ':'}

    visited = [False] * n
    ring_closure = {}   # ring_id → (src_atom_idx, bond_order)
    next_ring_id = 1
    result_parts = []

    def dfs(u, parent, incoming_order):
        nonlocal next_ring_id
        visited[u] = True

        # Emit bond symbol before atom (except for root)
        if parent is not None:
            # Check if this bond closes a ring
            for ring_id, (src_idx, b_order) in list(ring_closure.items()):
                if src_idx == u:
                    if b_order == 1:
                        result_parts.append('')
                    else:
                        result_parts.append(bond_sym.get(b_order, ''))
                    result_parts.append(f'%{ring_id}' if ring_id < 10 else f'%{ring_id}')
                    del ring_closure[ring_id]
                    # Don't continue traversal through a ring-closure bond
            else:
                result_parts.append(bond_sym.get(incoming_order, ''))

        # Emit atom
        # Handle two-letter elements: lowercase second letter
        elem = elements[u]
        result_parts.append(elem)

        # Collect unvisited neighbors
        unvisited = [(v, order) for v, order in adj[u] if not visited[v]]
        # Sort by bond order (heaviest first) and then by atomic number for consistency
        unvisited.sort(key=lambda x: (-x[1], elements[x[0]]))

        branch_count = 0
        for v, order in unvisited[:-1] if len(unvisited) > 1 else []:
            # Open branch for all but the last neighbor
            if branch_count == 0:
                result_parts.append('(')
            dfs(v, u, order)
            branch_count += 1

        # Process last unvisited neighbor (or only one) without branching
        if unvisited:
            dfs(unvisited[-1][0], u, unvisited[-1][1])
            if branch_count > 0:
                result_parts.append(')')
                branch_count = 0

        # Close any branches that were opened
        while branch_count > 0:
            result_parts.append(')')
            branch_count -= 1

        # Handle ring closures — mark visited neighbors that form rings
        for v, order in adj[u]:
            if v != parent and visited[v] and all(v != r_src for r_id, (r_src, _) in ring_closure.items()):
                # Found a ring closure bond — assign a ring ID
                rid = next_ring_id
                next_ring_id += 1
                ring_closure[rid] = (v, order)

    # Start DFS from the first atom
    dfs(0, None, 0)

    return ''.join(result_parts)


def extract_smiles_from_params(params_text):
    """Reconstruct SMILES from Rosetta .params file content.

    Parses ATOM lines for elements and BOND_TYPE lines for connectivity
    with proper bond orders (1=single, 2=double, 3=triple, 4=aromatic).
    Builds the molecular graph and returns SMILES using RDKit, CDPKit,
    or a pure-Python fallback.

    Returns the SMILES string, or None if reconstruction fails.
    """
    atoms_info = []   # list of (atom_name, element)
    bond_list = []    # list of (atom1, atom2, bond_order)

    for line in params_text.splitlines():
        line = line.strip()
        if not line:
            continue

        # Parse ATOM lines:  ATOM <name> <type> <mm_type> <charge>
        if line.startswith("ATOM "):
            parts = line.split()
            if len(parts) >= 4:
                atom_name = parts[1]
                ros_type = parts[2]
                elem = _rosetta_atom_type_to_element(ros_type)
                atoms_info.append((atom_name, elem))

        # Parse BOND_TYPE lines:  BOND_TYPE <atom1> <atom2> <order>
        elif line.startswith("BOND_TYPE "):
            parts = line.split()
            if len(parts) >= 4:
                try:
                    order = int(parts[3])
                    bond_list.append((parts[1], parts[2], order))
                except ValueError:
                    bond_list.append((parts[1], parts[2], 1))

    if not atoms_info or not bond_list:
        print(f"    Params has {len(atoms_info)} atoms, {len(bond_list)} bonds — insufficient to reconstruct SMILES")
        return None

    print(f"    Parsed {len(atoms_info)} atoms, {len(bond_list)} bonds from params")

    # Try RDKit first, then CDPKit, then pure Python
    smiles = _build_mol_from_params(atoms_info, bond_list)
    if smiles:
        print(f"    Generated SMILES via RDKit: {smiles}")
        return smiles
    smiles = _build_mol_cdpkit(atoms_info, bond_list)
    if smiles:
        print(f"    Generated SMILES via CDPKit: {smiles}")
        return smiles
    smiles = _build_smiles_pure_python(atoms_info, bond_list)
    if smiles:
        print(f"    Generated SMILES via pure Python: {smiles}")
        return smiles
    return None


def extract_smiles_via_obabel(enamine_path, ligand_name, chunk, subchunk,
                               tmp_root=None):
    """Extract the SDF from Enamine library and convert to SMILES with obabel."""
    superchunk_str, chunk_str = chunk_to_path(chunk)
    work_dir = tempfile.mkdtemp(dir=tmp_root) if tmp_root else tempfile.mkdtemp()

    try:
        # Extract the SDF for this ligand from the library
        # The SDF is typically in a separate .sdf.gz per subchunk
        sdf_gz = os.path.join(
            enamine_path, superchunk_str, chunk_str,
            f"condensed_params_and_db_{subchunk}.tar.gz"
        )
        # Alternative: try extracting db.db and using shapedb to get SMILES
        # For now, extract the SDF directly if available
        sdf_path = os.path.join(
            enamine_path, superchunk_str, chunk_str,
            f"{subchunk}.sdf.gz"
        )
        if not os.path.exists(sdf_path):
            # Try without .gz
            sdf_path = os.path.join(
                enamine_path, superchunk_str, chunk_str,
                f"{subchunk}.sdf"
            )
        if not os.path.exists(sdf_path):
            print(f"WARNING: SDF not found at {sdf_path}, cannot extract SMILES")
            return None

        # Extract only the ligand of interest from the SDF
        # Use zgrep to find the ligand, then obabel to convert
        extract_cmd = (
            f"zcat {sdf_path} 2>/dev/null | "
            f"obabel -isdf - -O - --separate -f 1 -l 1 "
            f"-osmi --title '{ligand_name}' 2>/dev/null | head -1"
        )
        if not sdf_path.endswith(".gz"):
            extract_cmd = extract_cmd.replace("zcat", "cat")

        result = subprocess.run(extract_cmd, shell=True, capture_output=True, text=True)
        if result.stdout.strip():
            smiles = result.stdout.strip().split()[0]
            return smiles

        # Fallback: grep for the molecule name in the SDF
        grep_cmd = (
            f"zgrep -A 100 '{ligand_name}' {sdf_path} 2>/dev/null | "
            f"obabel -isdf - -osmi 2>/dev/null | head -1"
        )
        if not sdf_path.endswith(".gz"):
            grep_cmd = grep_cmd.replace("zgrep", "grep")

        result = subprocess.run(grep_cmd, shell=True, capture_output=True, text=True)
        if result.stdout.strip():
            smiles = result.stdout.strip().split()[0]
            return smiles

        return None
    finally:
        shutil.rmtree(work_dir, ignore_errors=True)


def rebuild_smiles_from_params(params_text):
    """Alias for extract_smiles_from_params — kept for backward compatibility."""
    return extract_smiles_from_params(params_text)


# ═══════════════════════════════════════════════════════════════════════════════
# Step 3 — Generate conformers with CDPKit
# ═══════════════════════════════════════════════════════════════════════════════

def generate_cdpkit_conformers(smiles, num_conformers=150):
    """Generate conformers from SMILES using CDPKit.
    Returns a list of CDPKit Molecule objects (one per conformer)."""
    if not CDPKIT_AVAILABLE:
        print("ERROR: CDPKit not available")
        return []

    try:
        mol = Chem.parseSMILES(smiles)
        if mol is None:
            print(f"ERROR: CDPKit could not parse SMILES: {smiles}")
            return []

        # Add hydrogens for Rosetta compatibility
        Chem.makeHydrogenComplete(mol)

        # Perceive molecular properties required for conformer generation
        Chem.perceiveHybridizationStates(mol, True)
        Chem.calcImplicitHydrogenCounts(mol, True)
        ConfGen.prepareForConformerGeneration(mol)

        # Generate conformers into the molecule
        conf_gen = ConfGen.ConformerGenerator()
        settings = conf_gen.getSettings()
        settings.setMaxNumOutputConformers(num_conformers)
        settings.setMinRMSD(0.5)  # diversity threshold in Angstrom
        conf_gen.generate(mol)
        conf_gen.setConformers(mol)

        # Extract each conformer as a standalone molecule for SDF output.
        # Follows the gen_mol_ph4s.py pattern: iterate conformer indices and
        # use AtomConformer3DCoordinatesFunctor to access 3D coordinates.
        num_confs = Chem.getNumConformations(mol)
        conformers = []
        for conf_idx in range(num_confs):
            conf_mol = mol.clone()
            # set the 3D coordinate source to this conformer's coords
            Chem.applyConformation(conf_mol, conf_idx)
            conformers.append(conf_mol)

        print(f"  CDPKit generated {len(conformers)} conformers from SMILES")
        return conformers

    except Exception as e:
        print(f"ERROR: CDPKit conformer generation failed: {e}")
        return []


# ═══════════════════════════════════════════════════════════════════════════════
# Step 4 — Convert CDPKit conformers to Rosetta .params files
# ═══════════════════════════════════════════════════════════════════════════════

def write_cdpkit_conformer_to_sdf(conf_mol, sdf_path):
    """Write a single CDPKit conformer molecule to an SDF file."""
    if not CDPKIT_AVAILABLE:
        print("ERROR: CDPKit not available, cannot write SDF")
        return
    writer = Chem.FileSDFMolecularGraphWriter(sdf_path)
    writer.write(conf_mol)
    writer.close()


def convert_conformer_to_params(sdf_path, ligand_name, conf_idx, output_dir,
                                 realm_location):
    """Convert an SDF conformer to Rosetta .params using molfile_to_params.py."""
    params_name = f"{ligand_name}_{conf_idx}"
    # Resolve to absolute paths so the singularity container can see them
    sdf_abs = os.path.abspath(sdf_path)
    output_abs = os.path.abspath(output_dir)
    params_file = os.path.join(output_abs, f"{params_name}.params")

    # Use molfile_to_params via the Rosetta/singularity container
    conformator_sif = os.path.join(realm_location, "sif", "conformator_container.sif")
    if os.path.exists(conformator_sif):
        # Bind the directories containing the SDF and output so singularity
        # can access them (auto-mount only covers $HOME, /tmp, and the CWD).
        sdf_dir = os.path.dirname(sdf_abs)
        bind_paths = f"{realm_location},{sdf_dir},{output_abs}"
        params_cmd = (
            f"singularity exec --bind {sdf_dir} {conformator_sif} python "
            f"/conformator_for_container/molfile_to_params.py {sdf_abs} "
            f"-n {params_name} --keep-names --long-names --clobber --no-pdb"
        )
    else:
        # Fallback: try molfile_to_params directly (if Rosetta is in PATH)
        params_cmd = (
            f"molfile_to_params.py {sdf_abs} "
            f"-n {params_name} --keep-names --long-names --clobber --no-pdb"
        )

    rc = run_cmd(params_cmd, cwd=output_abs)
    if rc != 0 or not os.path.exists(params_file):
        print(f"WARNING: molfile_to_params failed for {params_name}")
        return None

    return params_file


# ═══════════════════════════════════════════════════════════════════════════════
# Step 5 — Score conformers with Rosetta
# ═══════════════════════════════════════════════════════════════════════════════

def score_single_conformer_with_rosetta(params_file, target_pdb, residue,
                                         motifs_file, discovery_root,
                                         atr, rep, ddg, work_dir,
                                         extra_params=None):
    """Run Rosetta scoring for a single conformer and return the score.
    Uses the Rosetta ligand docking script in score-only mode."""
    rosetta_script = os.path.join(
        discovery_root, "func", "discovery", "run_ligand_discovery_search.py"
    )
    extra_flag = f"--extra-args-file {extra_params}" if extra_params else ""

    # Create a minimal test_params dir with just this one ligand
    test_dir = os.path.join(work_dir, "test_params")
    os.makedirs(test_dir, exist_ok=True)
    run_cmd("touch exclude_pdb_component_list.txt patches.txt", cwd=test_dir)

    # Copy params file and create residue_types.txt
    shutil.copy2(params_file, test_dir)
    res_types_file = os.path.join(test_dir, "residue_types.txt")
    with open(res_types_file, "w") as fh:
        fh.write("## atom_type_set and mm-atom_type_set for Rosetta\n")
        fh.write("TYPE_SET_MODE full_atom\n")
        fh.write("ATOM_TYPE_SET fa_standard\n")
        fh.write("ELEMENT_SET default\n")
        fh.write("MM_ATOM_TYPE_SET fa_standard\n")
        fh.write("ORBITAL_TYPE_SET fa_standard\n")
        fh.write("## Params files\n")
        fh.write(f"{os.path.basename(params_file)}\n")

    # Run Rosetta (not via bsub — direct execution for scoring)
    cmd = (
        f"python {rosetta_script} {target_pdb} {residue} {motifs_file} "
        f"{test_dir} {discovery_root} {atr} {rep} {ddg} {extra_flag}"
    )
    print(f"  Running: {cmd[:120]}...")
    rc = run_cmd(cmd, cwd=work_dir)

    # Parse scores from output PDBs
    if parse_placement_scores is None or compute_weighted_total is None:
        print("WARNING: Scoring functions unavailable (utils.py not imported).")
        return {}

    scores = {}
    for root, dirs, files in os.walk(work_dir):
        for f in files:
            if f.endswith(".pdb") and "placed" in f and "minipose" not in f:
                pdb_path = os.path.join(root, f)
                parsed = parse_placement_scores(pdb_path)
                if parsed:
                    total, _ = compute_weighted_total(parsed, {})
                    lig_name = f.replace(".pdb", "")
                    scores[lig_name] = total

    return scores


# ═══════════════════════════════════════════════════════════════════════════════
# Main pipeline
# ═══════════════════════════════════════════════════════════════════════════════

def process_single_ligand(score, ligand_name, conf_num, chunk, subchunk,
                           target_pdb, anchor_residues, motifs_file,
                           output_dir, realm_location, enamine_path,
                           atr, rep, ddg, num_conformers, extra_params):
    """Full pipeline for one ligand entry."""
    print(f"\n{'='*60}")
    print(f"Processing: {ligand_name} (input conf {conf_num}, chunk {chunk})")
    print(f"{'='*60}")

    work_dir = tempfile.mkdtemp(dir=os.path.join(realm_location, "tmp"))
    params_dir = os.path.join(work_dir, "original_params")
    conf_dir = os.path.join(work_dir, "cdpkit_conformers")

    try:
        # ── Step A: Extract original params from Enamine library ─────────────
        print("  [Step A] Extracting params from Enamine library...")
        orig_params = extract_params_from_enamine(
            ligand_name, conf_num, chunk, subchunk,
            enamine_path, params_dir
        )
        if not orig_params:
            print(f"  ERROR: Failed to extract params for {ligand_name}")
            return []
        print(f"  Extracted: {orig_params}")

        # ── Step B: Extract SMILES ──────────────────────────────────────────
        print("  [Step B] Extracting SMILES...")
        with open(orig_params, "r") as fh:
            params_text = fh.read()

        smiles = extract_smiles_from_params(params_text)
        if not smiles:
            smiles = rebuild_smiles_from_params(params_text)
        if not smiles:
            smiles = extract_smiles_via_obabel(
                enamine_path, ligand_name, chunk, subchunk,
                tmp_root=os.path.join(realm_location, "tmp")
            )
        if not smiles:
            print(f"  ERROR: Could not extract SMILES for {ligand_name}")
            return []
        print(f"  SMILES: {smiles}")

        # ── Step C: Generate conformers with CDPKit ─────────────────────────
        print(f"  [Step C] Generating {num_conformers} conformers with CDPKit...")
        conformers = generate_cdpkit_conformers(smiles, num_conformers)
        if not conformers:
            return []
        print(f"  Generated {len(conformers)} conformers")

        # ── Step D: Convert each conformer to Rosetta .params ────────────────
        print("  [Step D] Converting conformers to Rosetta .params...")
        os.makedirs(conf_dir, exist_ok=True)
        params_files = []
        for i, conf_mol in enumerate(conformers):
            sdf_path = os.path.join(conf_dir, f"{ligand_name}_{i+1}.sdf")
            write_cdpkit_conformer_to_sdf(conf_mol, sdf_path)

            pf = convert_conformer_to_params(
                sdf_path, ligand_name, i+1, conf_dir, realm_location
            )
            if pf:
                params_files.append(pf)

        print(f"  Created {len(params_files)} Rosetta .params files")

        # ── Step E: Score each conformer with Rosetta ────────────────────────
        print("  [Step E] Scoring conformers with Rosetta...")
        all_scores = {}  # conformer_name -> rosetta_score
        for pf in params_files:
            conf_name = os.path.splitext(os.path.basename(pf))[0]
            score_dir = os.path.join(work_dir, f"score_{conf_name}")
            os.makedirs(score_dir, exist_ok=True)

            residue = anchor_residues.split(",")[0].strip()
            conf_scores = score_single_conformer_with_rosetta(
                pf, target_pdb, residue, motifs_file,
                realm_location, atr, rep, ddg, score_dir,
                extra_params=extra_params
            )
            all_scores.update(conf_scores)

        # ── Step F: Rank and return top conformers ──────────────────────────
        print(f"  [Step F] Ranking {len(all_scores)} scored conformers...")
        ranked = sorted(all_scores.items(), key=lambda x: x[1], reverse=True)
        return ranked

    finally:
        shutil.rmtree(work_dir, ignore_errors=True)


def main():
    parser = argparse.ArgumentParser(
        description="CDPKit conformer pipeline: pull params, generate "
                    "conformers with CDPKit, rank with Rosetta"
    )
    parser.add_argument("--input-list", required=True,
                        help="Shapedb results list (score,ligand_conf,chunk,subchunk)")
    parser.add_argument("--target-pdb", required=True, help="Target PDB file")
    parser.add_argument("--anchor-residues", required=True,
                        help="Anchor residues (e.g., 79 or 11,79,55)")
    parser.add_argument("--motifs-file", required=True, help="Motifs file path")
    parser.add_argument("--output-dir", required=True, help="Output directory")
    parser.add_argument("--realm-location",
                        default="/pi/summer.thyme-umw/Ji_rosetta_discovery",
                        help="Realm root directory")
    parser.add_argument("--enamine-path",
                        default="/pi/summer.thyme-umw/enamine-REAL-2.6billion",
                        help="Path to Enamine library")
    parser.add_argument("--num-conformers", type=int, default=150,
                        help="Conformers to generate per ligand (default: 150)")
    parser.add_argument("--top-n", type=int, default=20,
                        help="Number of top conformers to report (default: 20)")
    parser.add_argument("--atr", type=float, default=-2.0, help="fa_atr cutoff")
    parser.add_argument("--rep", type=float, default=150.0, help="fa_rep cutoff")
    parser.add_argument("--ddg", type=float, default=-9.0, help="ddg cutoff")
    parser.add_argument("--extra-params", default="",
                        help="Path to extra Rosetta arguments file")

    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    # ── Read input list ─────────────────────────────────────────────────────
    entries = []
    with open(args.input_list, "r") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            fields = line.split(",")
            if len(fields) < 4:
                print(f"Skipping malformed line: {line}")
                continue
            score, ligand_conf, chunk, subchunk = fields[:4]
            ligand_name = "_".join(ligand_conf.split("_")[:-1])
            conf_num = int(ligand_conf.split("_")[-1])
            entries.append((float(score), ligand_name, conf_num, chunk, subchunk))

    print(f"Loaded {len(entries)} entries from {args.input_list}")

    # ── Process each ligand ─────────────────────────────────────────────────
    global_ranked = []  # min-heap of (-rosetta_score, ligand, conf_idx)

    for entry_idx, (score, ligand_name, conf_num, chunk, subchunk) in enumerate(entries):
        print(f"\n{'#'*60}")
        print(f"Entry {entry_idx+1}/{len(entries)}: {ligand_name}")
        print(f"{'#'*60}")

        ranked = process_single_ligand(
            score, ligand_name, conf_num, chunk, subchunk,
            args.target_pdb, args.anchor_residues, args.motifs_file,
            args.output_dir, args.realm_location, args.enamine_path,
            args.atr, args.rep, args.ddg, args.num_conformers,
            args.extra_params
        )

        if not ranked:
            print(f"  WARNING: No scores for {ligand_name}")
            continue

        # Merge into global ranking
        for conf_name, conf_score in ranked:
            entry = (-conf_score, ligand_name, conf_name, conf_score)
            if len(global_ranked) < args.top_n:
                heapq.heappush(global_ranked, entry)
            elif conf_score > -global_ranked[0][0]:
                heapq.heapreplace(global_ranked, entry)

    # ── Write final results ─────────────────────────────────────────────────
    output_file = os.path.join(args.output_dir, "cdpkit_ranked_conformers.csv")
    sorted_results = sorted(global_ranked, reverse=True)

    with open(output_file, "w") as fh:
        fh.write("rank,ligand,conformer,rosetta_score,normalized_score\n")
        max_score = max(s[3] for s in sorted_results) if sorted_results else 1.0
        for rank, (neg_score, ligand, conf_name, conf_score) in enumerate(sorted_results, 1):
            norm = conf_score / max_score if max_score != 0 else 0.0
            fh.write(f"{rank},{ligand},{conf_name},{conf_score:.4f},{norm:.4f}\n")

    print(f"\n{'='*60}")
    print(f"Pipeline complete!")
    print(f"Results written to: {output_file}")
    print(f"Top {len(sorted_results)} conformers saved")
    if sorted_results:
        print(f"Best: {sorted_results[0][1]} {sorted_results[0][2]} "
              f"(score={sorted_results[0][3]:.4f})")
    print(f"{'='*60}")


if __name__ == "__main__":
    main()
