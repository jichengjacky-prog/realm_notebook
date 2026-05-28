#!/usr/bin/env python3
"""
HPC‑optimized rewrite of run_ligand_discovery_search.py

Key features:
- All heavy I/O performed inside node‑local scratch
- Automatic scratch directory detection with LSF‑safe fallback
- Clean subprocess usage
- No wildcard mv, no PANFS directory churn
- Deterministic, maintainable structure
"""

import os
import sys
import shutil
import subprocess
from pathlib import Path


############################################################
# Utility helpers
############################################################

def run(cmd: str):
    """Run a shell command safely with error checking."""
    subprocess.run(cmd, shell=True, check=True)


def copytree(src: Path, dst: Path):
    """Copy a directory tree with overwrite."""
    if dst.exists():
        shutil.rmtree(dst)
    shutil.copytree(src, dst)


############################################################
# Determine scratch directory
############################################################

scratch = os.environ.get("SCRATCH_DIR")

if not scratch:
    user = os.environ.get("USER", "unknown")
    jobid = os.environ.get("LSB_JOBID", "nojobid")
    jobindex = os.environ.get("LSB_JOBINDEX")

    if jobindex:
        scratch = f"/tmp/{user}/{jobid}_{jobindex}.scratch"
    else:
        scratch = f"/tmp/{user}/{jobid}.scratch"

SCRATCH_DIR = Path(scratch)
SCRATCH_DIR.mkdir(parents=True, exist_ok=True)

os.chdir(SCRATCH_DIR)

print("Operating out of directory: " + str(SCRATCH_DIR))

############################################################
# Parse arguments
############################################################

if len(sys.argv) < 8:
    print("Usage: run_ligand_discovery_search.py <pdb> <anchor> <motifs> <test_params> <atr> <rep> <ddg> [extra_args]")
    sys.exit(1)

target_pdb = Path(sys.argv[1])
anchor = sys.argv[2]
motifs_file = Path(sys.argv[3])
test_params_dir = Path(sys.argv[4])
atr, rep, ddg = sys.argv[5:8]
extra_args_file = Path(sys.argv[8]) if len(sys.argv) >= 9 else None

############################################################
# Prepare final destination and move to it
############################################################

final_dest = test_params_dir.parent / str(anchor)
final_dest.mkdir(exist_ok=True)
#os.chdir(final_dest)

############################################################
# Copy inputs into scratch
############################################################

local_pdb = SCRATCH_DIR / target_pdb.name
local_motifs = SCRATCH_DIR / motifs_file.name
local_params = SCRATCH_DIR / "test_params"

#test print of local_params residue_types
print("test print of local_params residue_types file")
os.system("ls " + str(local_params) + "/residue_types.txt")
os.system("cat " + str(local_params) + "/residue_types.txt")

shutil.copy(target_pdb, local_pdb)
shutil.copy(motifs_file, local_motifs)
copytree(test_params_dir, local_params)

if extra_args_file:
    local_extra = SCRATCH_DIR / "extra_args"
    shutil.copy(extra_args_file, local_extra)
else:
    local_extra = None


############################################################
# Write Rosetta args file
############################################################

args_path = SCRATCH_DIR / "args"

with args_path.open("w") as f:
    f.write(
        f"""
# Rosetta ligand discovery args
-constant_seed 1
-ignore_unrecognized_res
-in::file::override_database_params true
-constrain_relax_to_start_coords
-best_pdbs_to_keep 0

# Input files mapped into container
-s /input/{local_pdb.name}
-motif_filename /input/{local_motifs.name}
-params_directory_path /input/test_params/

# Anchor residue(s)
-protein_discovery_locus {anchor}

# Score cutoffs
-fa_atr_cutoff = {atr}
-fa_rep_cutoff = {rep}
-ddg_cutoff = {ddg}
"""
    )

    if local_extra:
        f.write("\n# Extra user args\n")
        f.write(local_extra.read_text())


############################################################
# Run Rosetta inside Singularity
############################################################
#sanity check of printing the arts to ensure successful mount
run(
    f"singularity exec "
    f"--bind {local_params}:/input/test_params "
    f"--bind {args_path}:/input/args "
    f"--bind {local_pdb}:/input/{local_pdb.name} "
    f"--bind {local_motifs}:/input/{local_motifs.name} "
    f"/pi/summer.thyme-umw/enamine-REAL-2.6billion/rosetta_condensed_6_25_2024.sif "
    f"cat /input/args"
)
run(
    f"singularity exec "
    f"--bind {local_params}:/input/test_params "
    f"--bind {args_path}:/input/args "
    f"--bind {local_pdb}:/input/{local_pdb.name} "
    f"--bind {local_motifs}:/input/{local_motifs.name} "
    f"/pi/summer.thyme-umw/enamine-REAL-2.6billion/rosetta_condensed_6_25_2024.sif "
    f"ls /input/test_params"
)
run(
    f"singularity exec "
    f"--bind {local_params}:/input/test_params "
    f"--bind {args_path}:/input/args "
    f"--bind {local_pdb}:/input/{local_pdb.name} "
    f"--bind {local_motifs}:/input/{local_motifs.name} "
    f"/pi/summer.thyme-umw/enamine-REAL-2.6billion/rosetta_condensed_6_25_2024.sif "
    f"cat /input/test_params/residue_types.txt"
)

run(
    f"singularity exec "
    f"--bind {local_params}:/input/test_params "
    f"--bind {args_path}:/input/args "
    f"--bind {local_pdb}:/input/{local_pdb.name} "
    f"--bind {local_motifs}:/input/{local_motifs.name} "
    f"/pi/summer.thyme-umw/enamine-REAL-2.6billion/rosetta_condensed_6_25_2024.sif "
    f"/rosetta/source/bin/ligand_discovery_search_protocol.linuxgccrelease @/input/args"
)


############################################################
# Process Rosetta outputs
############################################################

placements = SCRATCH_DIR / "placements"
placements.mkdir(exist_ok=True)

# Move and rename PDBs deterministically
for pdb in SCRATCH_DIR.glob("*.pdb"):
    newname = placements / f"res{anchor}_{pdb.name}"
    pdb.rename(newname)

#move into placements
os.chdir("placements")

# Score placements
run(
    "python /pi/summer.thyme-umw/enamine-REAL-2.6billion/umass_chan_REAL-M_platform/rosetta/score_placed_ligands_with_filtering.py",
)

#move back up
os.chdir("..")

# Copy CSVs up one level
for csv in placements.glob("*.csv"):
    shutil.copy(csv, SCRATCH_DIR)

# Tar placements
run(f"tar -czf {SCRATCH_DIR}/placements.tar.gz -C {SCRATCH_DIR} placements")

############################################################
# Dehydrate the placements to reduce overhead
############################################################

#run the dehydrate script to minimize overhead (until it is time to process the discovery results)
run("python /pi/summer.thyme-umw/enamine-REAL-2.6billion/umass_chan_REAL-M_platform/tidying/shrink_placement_pdbs_to_placement_and_surrounding_residues.py " + str(target_pdb))

############################################################
# Copy final results back to submission directory
############################################################

#final_dest = Path.cwd()
final_dest = test_params_dir.parent / str(anchor)
final_dest.mkdir(exist_ok=True)

shutil.copy(SCRATCH_DIR / "placements.tar.gz", final_dest)

for csv in SCRATCH_DIR.glob("*.csv"):
    shutil.copy(csv, final_dest)

