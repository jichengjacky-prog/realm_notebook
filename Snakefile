"""
Snakemake pipeline for the Rosetta Ligand Discovery workflow.

Converts the unified_discovery_controller.py into a reproducible,
parallelisable Snakemake pipeline with checkpoint/resume support.

Usage:
    snakemake --cores 16 --configfile config.yaml

    # Or override config on the command line:
    snakemake --cores 16 \
        --config shapedb_list=results.txt target_pdb=receptor.pdb \
        anchor_residues="239,242" motifs_file=motifs.motifs \
        output_dir=output/

Notes:
Before running, make sure you have CDPKit installed and properly licensed, and that the Enamine library is
downloaded and accessible at the specified path. Also ensure that Rosetta is set up for command-line use and 
that the discovery scripts are in place. The config.yaml file should be populated with the appropriate parameters 
for your run, or you can override them via --config on the command line as shown above. 
"""

import os
import sys
import subprocess
import time
import csv
from pathlib import Path

# ═══════════════════════════════════════════════════════════════════════════════
# Configuration — set in config.yaml or via --config
# ═══════════════════════════════════════════════════════════════════════════════

configfile: "config.yaml"

# Defaults (overridden by config.yaml / --config)
SNAKEFILE_DIR  = workflow.basedir if hasattr(workflow, 'basedir') else os.getcwd()
SHAPEDB_LIST    = config.get("shapedb_list",    "")
TARGET_PDB      = config.get("target_pdb",      "")
ANCHOR_RESIDUES = config.get("anchor_residues", "")
MOTIFS_FILE     = config.get("motifs_file",     "")
OUTPUT_DIR      = config.get("output_dir",      "output")
REALM_LOCATION  = config.get("realm_location",  SNAKEFILE_DIR)
ENAMINE_PATH    = config.get("enamine_path",    "/pi/summer.thyme-umw/enamine-REAL-2.6billion")
NUM_CONFORMERS  = config.get("num_conformers",  400)
TOP_N           = config.get("top_n",           10)
BATCH_SIZE      = config.get("batch_size",      10)
MAX_LIGANDS     = config.get("max_ligands",     10)
TOP_HITS        = config.get("top_hits",        0)
ATR             = config.get("atr",             -2.0)
REP             = config.get("rep",             150.0)
DDG             = config.get("ddg",             -9.0)
EXTRA_PARAMS    = config.get("extra_params",    "extra_arg/sample_extra_arg")
LICENSE_KEY     = config.get("license_key",     "")
MIN_MOTIF_RATIO = config.get("min_motif_ratio", 0.25)
WEIGHTS_FILE    = config.get("weights_file") or os.path.join(REALM_LOCATION, "function", "discovery", "score_weights.json")
TMP_ROOT        = config.get("tmp_root")     or os.path.join(OUTPUT_DIR, "tmp")

os.makedirs(OUTPUT_DIR, exist_ok=True)
os.makedirs(TMP_ROOT, exist_ok=True)

# ── LSF cluster resource defaults ───────────────────────────────────────────
# Override in config.yaml or via --config
LSF_QUEUE_ROSETTA   = config.get("lsf_queue_rosetta",   "long")
LSF_QUEUE_DEFAULT   = config.get("lsf_queue_default",   "short")
LSF_WALLTIME_ROSETTA = config.get("lsf_walltime_rosetta", "12:00")
LSF_WALLTIME_DEFAULT = config.get("lsf_walltime_default", "2:00")

# ═══════════════════════════════════════════════════════════════════════════════
# Parse shapedb list → batch files
# ═══════════════════════════════════════════════════════════════════════════════

BATCHES_DIR = os.path.join(OUTPUT_DIR, "batches")
BATCH_SHARD_SIZE = 1000   # files per subdirectory — keeps NFS directories small
PARSE_DONE_FLAG = os.path.join(OUTPUT_DIR, ".parse_shapedb.done")
os.makedirs(BATCHES_DIR, exist_ok=True)

def batch_file_path(batch_id):
    """Return the sharded path for a batch file.
    Sharding (1000 files per subdir) keeps NFS directories fast even with
    hundreds of thousands of batches."""
    shard = int(batch_id) // BATCH_SHARD_SIZE
    shard_dir = os.path.join(BATCHES_DIR, f"{shard:04d}")
    return os.path.join(shard_dir, f"batch_{batch_id}.txt")

def _discover_batch_ids():
    """Walk BATCHES_DIR and return a sorted list of existing batch IDs."""
    ids = []
    for root, dirs, files in os.walk(BATCHES_DIR):
        for f in files:
            if f.startswith("batch_") and f.endswith(".txt"):
                try:
                    ids.append(int(f[len("batch_"):-len(".txt")]))
                except ValueError:
                    pass
    return sorted(ids)

if os.path.exists(PARSE_DONE_FLAG):
    # ── Skip parsing; reconstruct BATCH_IDS from existing batch files ──
    batch_ids = _discover_batch_ids()
    print(f"Skipped shapedb parse ({PARSE_DONE_FLAG} exists) — "
          f"found {len(batch_ids)} existing batches "
          f"({(len(batch_ids) + BATCH_SHARD_SIZE - 1) // BATCH_SHARD_SIZE} shards)")
else:
    # ── Read the shapedb list ───────────────────────────────────────────
    entries = []
    with open(SHAPEDB_LIST, "r") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if TOP_HITS > 0 and len(entries) >= TOP_HITS:
                break
            fields = [f.strip() for f in line.split(",")]
            if len(fields) < 4:
                continue
            score, ligand_conf, chunk, subchunk = fields[:4]
            ligand_name = "_".join(ligand_conf.split("_")[:-1])
            conf_num = int(ligand_conf.split("_")[-1])
            entries.append((float(score), ligand_name, conf_num, chunk, subchunk))

    # ── Split into batches, write batch files (sharded into subdirectories) ─
    batch_ids = []
    batch = []
    for i, (score, lig, cn, ch, sc) in enumerate(entries):
        batch.append((lig, cn, ch, sc))
        if len(batch) >= BATCH_SIZE or i == len(entries) - 1:
            batch_id = len(batch_ids)
            bf_path = batch_file_path(batch_id)
            os.makedirs(os.path.dirname(bf_path), exist_ok=True)
            with open(bf_path, "w") as fh:
                for l, c, h, s in batch:
                    fh.write(f"{l},{c},{h},{s}\n")
            batch_ids.append(batch_id)
            batch = []

    print(f"Split {len(entries)} entries into {len(batch_ids)} batches "
          f"({(len(batch_ids) + BATCH_SHARD_SIZE - 1) // BATCH_SHARD_SIZE} shards)")

    # ── Write done flag ─────────────────────────────────────────────────
    with open(PARSE_DONE_FLAG, "w") as f:
        f.write(f"{len(batch_ids)} batches, {len(entries)} entries\n")

# Wildcard: batch ID
BATCH_IDS = batch_ids

# ── Round 2 sub-batch splitting ─────────────────────────────────────────────
#   max_confs_per_sub — hard ceiling: conformers per sub-batch (default 10).
#                       Number of sub-batches scales up automatically to
#                       ceil(total_confs / max_confs_per_sub).  Sub-batches
#                       are discovered dynamically via a Snakemake checkpoint
#                       so only the actually-needed Rosetta jobs are submitted.
MAX_CONFS_PER_SUB = config.get("max_confs_per_sub", 10)

# ═══════════════════════════════════════════════════════════════════════════════
# Rules
# ═══════════════════════════════════════════════════════════════════════════════

rule all:
    input:
        os.path.join(OUTPUT_DIR, "cleanup.done"),

# ── Step 1: Extract library params → test_params directory ──────────────

rule extract_params:
    """Extract conformer params from the Enamine library for one batch."""
    input:
        batch_file = lambda wildcards: batch_file_path(wildcards.batch_id),
    output:
        done_flag        = os.path.join(TMP_ROOT, "batch_{batch_id}", "extract_params.done"),
        test_params_done = os.path.join(TMP_ROOT, "batch_{batch_id}", "test_params", ".round1_ready"),
    params:
        realm_location = REALM_LOCATION,
        enamine_path   = ENAMINE_PATH,
        tmp_root       = TMP_ROOT,
    resources:
        mem_mb=2000,
        cpus=1,
        queue=LSF_QUEUE_DEFAULT,
        walltime=LSF_WALLTIME_DEFAULT,
    log:
        os.path.join(TMP_ROOT, "batch_{batch_id}", "extract_params.log"),
    shell:
        """
        set -e
        BATCH_DIR=$(dirname {output.done_flag})
        TEST_PARAMS_DIR=$BATCH_DIR/test_params
        mkdir -p $BATCH_DIR

        /home/ji.cheng4-umw/miniforge3/envs/realm_env/bin/python3.11 -c "
import sys; sys.path.insert(0, '{params.realm_location}/function/discovery')
from utils import create_test_params_dir

ligands = []
with open('{input.batch_file}') as fh:
    for line in fh:
        parts = line.strip().split(',')
        if len(parts) == 4:
            ligands.append((parts[0], int(parts[1]), parts[2], parts[3]))

create_test_params_dir(ligands, '$BATCH_DIR', '{params.realm_location}',
                       '{params.enamine_path}', tmp_root='{params.tmp_root}')
print('extract_params complete')
" > {log} 2>&1
        touch {output.test_params_done}
        touch {output.done_flag}
        """


# ── Step 2: Run Rosetta discovery (bsub) ────────────────────────────────

rule rosetta_discovery_round1:
    """Run Rosetta ligand discovery search for one batch (round 1 — library conformers)."""
    input:
        params_done      = os.path.join(TMP_ROOT, "batch_{batch_id}", "extract_params.done"),
        test_params_done = os.path.join(TMP_ROOT, "batch_{batch_id}", "test_params", ".round1_ready"),
    output:
        done_flag   = os.path.join(TMP_ROOT, "batch_{batch_id}", "rosetta_round1.done"),
    params:
        target_pdb      = TARGET_PDB,
        anchor_residues = ANCHOR_RESIDUES,
        motifs_file     = MOTIFS_FILE,
        discovery_root  = REALM_LOCATION,
        atr             = str(ATR),
        rep             = str(REP),
        ddg             = str(DDG),
        extra_params    = EXTRA_PARAMS,
        tmp_root        = TMP_ROOT,
    resources:
        mem_mb=8000,
        cpus=1,
        queue=LSF_QUEUE_ROSETTA,
        walltime=LSF_WALLTIME_ROSETTA,
    log:
        os.path.join(TMP_ROOT, "batch_{batch_id}", "rosetta_round1.log"),
    shell:
        """
        set -e
        BATCH_DIR=$(dirname {output.done_flag})
        mkdir -p $BATCH_DIR

        /home/ji.cheng4-umw/miniforge3/envs/realm_env/bin/python3.11 -c "
import sys, os, glob
sys.path.insert(0, '{params.discovery_root}/function/discovery')
from utils import run_rosetta_discovery_search

any_placements = False
for residue in '{params.anchor_residues}'.split(','):
    res = residue.strip()
    res_dir = os.path.join('$BATCH_DIR', 'round1', res)
    os.makedirs(res_dir, exist_ok=True)
    success = run_rosetta_discovery_search(
        '{params.target_pdb}', res, '{params.motifs_file}',
        '$BATCH_DIR/test_params/', '{params.discovery_root}',
        '{params.atr}', '{params.rep}', '{params.ddg}',
        extra_args_file='{params.extra_params}' if '{params.extra_params}' else None,
        work_dir=res_dir
    )
    if not success:
        # Rosetta may return non-zero even when placements were produced.
        # Check if any placement PDBs or tarballs exist before giving up.
        pdbs = glob.glob(os.path.join(res_dir, '*.pdb'))
        tars = glob.glob(os.path.join(res_dir, '*.tar.gz'))
        if pdbs or tars:
            print(f'Rosetta exited non-zero but placements exist for residue {{res}} ({{len(pdbs)}} PDBs, {{len(tars)}} tarballs)')
            any_placements = True
        else:
            print(f'ERROR: Rosetta discovery failed for residue {{res}} with no placements')
    else:
        any_placements = True
        

if not any_placements:
    raise SystemExit(1)
open('{output.done_flag}', 'w').close()

print('Rosetta round 1 complete')
" > {log} 2>&1
        """


# ── Step 3: Score placements from round 1 ───────────────────────────────

rule score_round1:
    """Score placements from round 1 and write per-batch scores CSV."""
    input:
        done_flag   = os.path.join(TMP_ROOT, "batch_{batch_id}", "rosetta_round1.done"),
    output:
        scores_csv  = os.path.join(TMP_ROOT, "batch_{batch_id}", "scores_round1.csv"),
    params:
        weights_file = WEIGHTS_FILE,
        realm        = REALM_LOCATION,
        batches_dir  = BATCHES_DIR,
    resources:
        mem_mb=2000,
        cpus=1,
        queue=LSF_QUEUE_DEFAULT,
        walltime=LSF_WALLTIME_DEFAULT,
    log:
        os.path.join(TMP_ROOT, "batch_{batch_id}", "score_round1.log"),
    shell:
        """
        set -e
        BATCH_DIR=$(dirname {output.scores_csv})

        /home/ji.cheng4-umw/miniforge3/envs/realm_env/bin/python3.11 -c "
import sys, os, csv, glob
sys.path.insert(0, '{params.realm}/function/discovery')
from utils import score_placements

# ── Load batch ligand list for dedup matching ────────────────────────
batch_id = int('{wildcards.batch_id}')
batch_file = os.path.join('{params.batches_dir}', f'{{batch_id // 1000:04d}}', f'batch_{{batch_id}}.txt')
batch_ligands = set()
with open(batch_file, 'r') as fh:
    for line in fh:
        parts = line.strip().split(',')
        if len(parts) >= 1:
            batch_ligands.add(parts[0])

def get_base_ligand(placement_name):
    '''Match placement_name against known batch ligands (longest match first).'''
    for lig in sorted(batch_ligands, key=len, reverse=True):
        if lig in placement_name:
            return lig
    # Fallback: strip last underscore-separated field
    return '_'.join(placement_name.split('_')[:-1])

scores = score_placements(os.path.join('$BATCH_DIR', 'round1'), weights_file='{params.weights_file}')
# Also extract real_motif_ratio from weighted_scores.csv files
motif_ratios = {{}}
for csv_path in glob.glob(os.path.join('$BATCH_DIR', 'round1', '*', 'weighted_scores.csv')):
    with open(csv_path, 'r') as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            pdb_name = os.path.basename(row.get('file', ''))
            lig_name = pdb_name.replace('.pdb', '')
            if lig_name:
                try:
                    motif_ratios[lig_name] = float(row.get('real_motif_ratio', 0))
                except ValueError:
                    pass

# ── Deduplicate: keep only best-scoring conformer per base ligand ─
best = {{}}  # base_ligand -> (score, full_ligand_name, motif_ratio)
for lig, sc in scores.items():
    base = get_base_ligand(lig)
    if base not in best or sc > best[base][0]:
        best[base] = (sc, lig, motif_ratios.get(lig, 0.0))

with open('{output.scores_csv}', 'w', newline='') as fh:
    writer = csv.writer(fh)
    writer.writerow(['ligand', 'score', 'real_motif_ratio'])
    for base, (sc, lig, mr) in sorted(best.items(), key=lambda x: -x[1][0]):
        writer.writerow([lig, f'{{sc:.6f}}', f'{{mr:.6f}}'])
print(f'Scored {{len(best)}} unique ligands (from {{len(scores)}} placements)')
" > {log} 2>&1

        # ── Delete round1 placement files immediately after scoring ─────────
        # Scores are saved; placement PDBs & tarballs are no longer needed.
        rm -rf "$BATCH_DIR"/round1
        rm -f "$BATCH_DIR"/rosetta_round1.log
        """


# ── Step 3.5: Filter top-N hits from round 1 & cleanup intermediates ────

rule filter_top_round1:
    """Aggregate round-1 scores from all batches, keep top-N ligands via heapq,
    write a filter list, and clean up round-1 residue directories for ALL
    batches (placements already scored in scores_round1.csv) to prevent
    mixing with round-2 outputs.  Preserves .done and .csv files."""
    localrule: True
    input:
        expand(os.path.join(TMP_ROOT, "batch_{batch_id}", "scores_round1.csv"),
               batch_id=BATCH_IDS),
    output:
        top_list = os.path.join(OUTPUT_DIR, "top_ligands_round1.txt"),
    params:
        max_ligands     = MAX_LIGANDS,
        min_motif_ratio = MIN_MOTIF_RATIO,
        tmp_root        = TMP_ROOT,
    resources:
        mem_mb=4000,
        cpus=1,
    log:
        os.path.join(OUTPUT_DIR, "filter_top_round1.log"),
    shell:
        """
        set -e
        /home/ji.cheng4-umw/miniforge3/envs/realm_env/bin/python3.11 -c "
import heapq, csv, shutil, glob, os

# ── Collect all scores, deduplicate by base ligand (keep best conformer) ─
# Build a global set of batch ligand names for matching placement names.
# tmp_root is {{output_dir}}/tmp, so output_dir = tmp_root/../..
output_dir = os.path.dirname(os.path.dirname('{params.tmp_root}'))
batches_dir = os.path.join(output_dir, 'batches')
batch_ligands_global = set()
for batch_dir in glob.glob(os.path.join('{params.tmp_root}', 'batch_*')):
    batch_name = os.path.basename(batch_dir)  # e.g. "batch_0"
    batch_id = int(batch_name.split('_')[-1])
    batch_file = os.path.join(batches_dir, f'{{batch_id // 1000:04d}}', f'batch_{{batch_id}}.txt')
    if os.path.exists(batch_file):
        with open(batch_file, 'r') as fh:
            for line in fh:
                parts = line.strip().split(',')
                if len(parts) >= 1:
                    batch_ligands_global.add(parts[0])

def get_base_ligand_global(placement_name):
    '''Match placement_name against all known batch ligands (longest match first).'''
    for lig in sorted(batch_ligands_global, key=len, reverse=True):
        if lig in placement_name:
            return lig
    return '_'.join(placement_name.split('_')[:-1])

seen = {{}}  # base_ligand -> (best_score, full_ligand_name, real_motif_ratio)
for sf in glob.glob(os.path.join('{params.tmp_root}', 'batch_*', 'scores_round1.csv')):
    if not os.path.exists(sf):
        continue
    with open(sf, 'r') as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            lig = row['ligand']
            score = float(row['score'])
            mr = float(row['real_motif_ratio'])
            if mr <= {params.min_motif_ratio}:
                continue
            base = get_base_ligand_global(lig)
            if base not in seen or score > seen[base][0]:
                seen[base] = (score, lig, mr)

# ── Build top-N heap from unique ligands ─────────────────────────
heap = []
for base, (score, lig, mr) in seen.items():
    entry = (-score, lig, score, mr)
    if len(heap) < {params.max_ligands}:
        heapq.heappush(heap, entry)
    elif score > -heap[0][0]:
        heapq.heapreplace(heap, entry)

top_ligands = {{lig for _, lig, _, _ in heap}}
print(f'  Top {{len(top_ligands)}} unique ligands kept from round 1')

# ── Clean up round-1 residue directories for ALL batches ────────
for batch_dir in glob.glob(os.path.join('{params.tmp_root}', 'batch_*')):
    scores_csv = os.path.join(batch_dir, 'scores_round1.csv')
    if not os.path.exists(scores_csv):
        continue
    keep = False
    with open(scores_csv, 'r') as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            if row['ligand'] in top_ligands:
                keep = True
                break
    round1_dir = os.path.join(batch_dir, 'round1')
    if os.path.isdir(round1_dir):
        shutil.rmtree(round1_dir, ignore_errors=True)
    if not keep:
        print(f'  Cleaned intermediate files for {{os.path.basename(batch_dir)}} (non-top)')
    else:
        print(f'  Cleaned round-1 placements for {{os.path.basename(batch_dir)}} (top batch, kept scores)')

# ── Write filter list (ligand,score,real_motif_ratio) ─────────
# Build a lookup: full_ligand_name -> (score, real_motif_ratio)
top_info = {{}}
for _, lig, score, mr in heap:
    top_info[lig] = (score, mr)
with open('{output.top_list}', 'w') as fh:
    fh.write('ligand,score,real_motif_ratio\\n')
    for lig in sorted(top_info.keys()):
        score, mr = top_info[lig]
        fh.write(f'{{lig}},{{score:.6f}},{{mr:.6f}}\\n')
print(f'  Wrote {{len(top_info)}} top ligands to {output.top_list}')
" > {log} 2>&1
        """


# ── Step 4: Generate CDPKit conformers & add to test_params ─────────────

rule generate_conformers:
    """Re-extract top-N ligand params from Enamine library, generate CDPKit
    conformers, and create a fresh test_params directory for round 2.

    Disregards the round-1 extract_params.done signal — re-extracts from
    scratch for only the top-N ligands."""
    input:
        batch_file  = lambda wildcards: batch_file_path(wildcards.batch_id),
        top_list    = os.path.join(OUTPUT_DIR, "top_ligands_round1.txt"),
    output:
        done_flag        = os.path.join(TMP_ROOT, "batch_{batch_id}", "conformers_generated.done"),
        test_params_done = os.path.join(TMP_ROOT, "batch_{batch_id}", "test_params", ".round2_ready"),
    resources:
        mem_mb=4000,
        cpus=1,
        queue=LSF_QUEUE_DEFAULT,
        walltime="2:00",
    params:
        realm_location = REALM_LOCATION,
        enamine_path   = ENAMINE_PATH,
        num_conformers = NUM_CONFORMERS,
        license_key    = LICENSE_KEY,
        tmp_root       = TMP_ROOT,
    log:
        os.path.join(TMP_ROOT, "batch_{batch_id}", "generate_conformers.log"),
    shell:
        """
        set -e
        BATCH_DIR=$(dirname {output.done_flag})
        TEST_PARAMS_DIR=$BATCH_DIR/test_params
        mkdir -p $BATCH_DIR
        # Clean test_params for round 2 but keep .round1_ready sentinel
        # so Snakemake still sees extract_params as complete.
        mkdir -p "$TEST_PARAMS_DIR"
        find "$TEST_PARAMS_DIR" -mindepth 1 ! -name '.round1_ready' -exec rm -rf {{}} +
        # Also clean any lingering round-1 residue directories so round-2
        # placements never mix with round-1 placements (belt & suspenders).
        rm -rf "$BATCH_DIR"/round1

        /home/ji.cheng4-umw/miniforge3/envs/realm_env/bin/python3.11 -c "
import sys, os
sys.path.insert(0, '{params.realm_location}/function/discovery')
from utils import (
    create_test_params_dir,
    generate_and_add_conformers_to_test_params,
)

# ── Load top-N filter list (CSV: ligand,score,real_motif_ratio) ───
top_set = set()
with open('{input.top_list}') as fh:
    header = fh.readline()  # skip header
    for line in fh:
        line = line.strip()
        if line:
            name = line.split(',')[0]
            top_set.add(name)
print(f'Top-N filter: {{len(top_set)}} ligands')

# ── Read batch file and keep only top-N ─────────────────────────────
# top_set contains full placement names (e.g. res242_..._PV-003529332293_3_91)
# Batch file has base ligand names (e.g. PV-003529332293) — use substring match.
ligands = []
with open('{input.batch_file}') as fh:
    for line in fh:
        parts = line.strip().split(',')
        if len(parts) == 4:
            name = parts[0]
            if any(name in entry for entry in top_set):
                ligands.append((name, int(parts[1]), parts[2], parts[3]))
print(f'Top-N ligands in this batch: {{len(ligands)}}')

if not ligands:
    print('No top-N ligands in this batch, skipping')
    os.makedirs('$TEST_PARAMS_DIR', exist_ok=True)
    open('{output.test_params_done}', 'w').close()
    open('{output.done_flag}', 'w').close()
    sys.exit(0)

# ── Step A: Re-extract library params for top ligands ───────────────
print('Re-extracting params from Enamine library for top ligands...')
create_test_params_dir(
    ligands, '$BATCH_DIR', '{params.realm_location}',
    '{params.enamine_path}', tmp_root='{params.tmp_root}'
)

# ── Step B: Generate CDPKit conformers and add to test_params ───────
print('Generating CDPKit conformers for top ligands...')
ok = generate_and_add_conformers_to_test_params(
    '$TEST_PARAMS_DIR', ligands,
    '{params.realm_location}', '{params.enamine_path}',
    num_conformers={params.num_conformers},
    tmp_root='{params.tmp_root}'
)
if ok:
    open('{output.test_params_done}', 'w').close()
    open('{output.done_flag}', 'w').close()
    print('Conformer generation complete.')
else:
    print('ERROR: Conformer generation failed')
    import sys; sys.exit(1)
" > {log} 2>&1
        """


# ── Step 4.5: Split round-2 conformers into sub-batches (checkpoint) ───

checkpoint split_conformers_round2:
    """Split round-2 .params files into round2/<sub>/test_params/ sub-dirs.
    This is a Snakemake checkpoint — downstream rules discover the actual
    sub-batch count from its output, so only needed Rosetta jobs are submitted."""
    input:
        confs_done       = os.path.join(TMP_ROOT, "batch_{batch_id}", "conformers_generated.done"),
        test_params_done = os.path.join(TMP_ROOT, "batch_{batch_id}", "test_params", ".round2_ready"),
    output:
        done_flag  = os.path.join(TMP_ROOT, "batch_{batch_id}", "split_round2.done"),
        num_subs   = os.path.join(TMP_ROOT, "batch_{batch_id}", "num_subs.txt"),
    params:
        max_per_sub = MAX_CONFS_PER_SUB,
    resources:
        mem_mb=500,
        cpus=1,
        queue=LSF_QUEUE_DEFAULT,
        walltime="0:30",
    log:
        os.path.join(TMP_ROOT, "batch_{batch_id}", "split_round2.log"),
    shell:
        """
        set -e
        BATCH_DIR=$(dirname {output.done_flag})
        SRC_DIR="$BATCH_DIR/test_params"

        # Count .params files (excluding sentinel / non-params)
        params_files=($(find "$SRC_DIR" -maxdepth 1 -name '*.params' | sort))
        total_params=${{#params_files[@]}}
        echo "Found $total_params .params files in test_params/"

        if [ "$total_params" -eq 0 ]; then
            echo "No params files — nothing to split"
            echo "0" > {output.num_subs}
            touch {output.done_flag}
            exit 0
        fi

        # ── Scale sub-batches to fit: keep per_sub fixed, grow num_subs ─
        per_sub={params.max_per_sub}           # hard ceiling per sub-batch
        num_subs=$(( (total_params + per_sub - 1) / per_sub ))
        echo "$total_params conformers → $num_subs sub-batch(es) × ≤$per_sub each"

        # ── Capture residue_types.txt header (first 7 lines) ──────────
        RT_HEADER=$(head -7 "$SRC_DIR/residue_types.txt" 2>/dev/null)

        # ── Distribute into round2/<sub>/test_params/ ────────────────────
        sub=0
        count=0
        for pf in "${{params_files[@]}}"; do
            SUB_TEST_PARAMS="$BATCH_DIR/round2/$sub/test_params"
            mkdir -p "$SUB_TEST_PARAMS"
            cp "$pf" "$SUB_TEST_PARAMS/"
            count=$((count + 1))
            if [ $count -ge $per_sub ]; then
                # Copy static support files
                for sf in patches.txt exclude_pdb_component_list.txt; do
                    [ -f "$SRC_DIR/$sf" ] && cp "$SRC_DIR/$sf" "$SUB_TEST_PARAMS/" || true
                done
                # Generate filtered residue_types.txt for this sub-batch
                echo "$RT_HEADER" > "$SUB_TEST_PARAMS/residue_types.txt"
                (cd "$SUB_TEST_PARAMS" && ls -1 *.params) >> "$SUB_TEST_PARAMS/residue_types.txt"
                touch "$SUB_TEST_PARAMS/.round2_ready"
                sub=$((sub + 1))
                count=0
            fi
        done

        # Last (possibly partial) sub-batch
        if [ $count -gt 0 ]; then
            SUB_TEST_PARAMS="$BATCH_DIR/round2/$sub/test_params"
            for sf in patches.txt exclude_pdb_component_list.txt; do
                [ -f "$SRC_DIR/$sf" ] && cp "$SRC_DIR/$sf" "$SUB_TEST_PARAMS/" || true
            done
            echo "$RT_HEADER" > "$SUB_TEST_PARAMS/residue_types.txt"
            (cd "$SUB_TEST_PARAMS" && ls -1 *.params) >> "$SUB_TEST_PARAMS/residue_types.txt"
            touch "$SUB_TEST_PARAMS/.round2_ready"
            sub=$((sub + 1))
        fi

        echo "Split $total_params params into $sub sub-batch(es) under round2/"
        echo "$sub" > {output.num_subs}

        # Clean up original .params files from test_params/ — they're now
        # distributed into round2/<sub>/test_params/ and no longer needed.
        find "$SRC_DIR" -maxdepth 1 -name '*.params' -exec rm {{}} +
        echo "Cleaned original .params files from $SRC_DIR"

        touch {output.done_flag}
        """


# ── Dynamic sub-batch discovery (reads checkpoint output) ───────────────

def get_round2_sub_ids(wildcards):
    """Return the list of sub-batch IDs that were actually created.
    Called by Snakemake after the split checkpoint completes."""
    chk = checkpoints.split_conformers_round2.get(batch_id=wildcards.batch_id)
    num_subs_path = chk.output.num_subs
    with open(num_subs_path) as f:
        n = int(f.read().strip())
    return list(range(n))


# ── Step 5: Rerun Rosetta discovery with new conformers (sub-batched) ──

rule rosetta_discovery_round2:
    """Run Rosetta discovery again with CDPKit-generated conformers,
    one job per sub-batch (sub_id).  Sub-batches are discovered
    dynamically — only actually-populated sub-IDs are submitted."""
    input:
        split_done = os.path.join(TMP_ROOT, "batch_{batch_id}", "split_round2.done"),
    output:
        done_flag = os.path.join(TMP_ROOT, "batch_{batch_id}", "rosetta_round2_{sub_id}.done"),
    resources:
        mem_mb=8000,
        cpus=1,
        queue=LSF_QUEUE_ROSETTA,
        walltime=LSF_WALLTIME_ROSETTA,
    params:
        target_pdb      = TARGET_PDB,
        anchor_residues = ANCHOR_RESIDUES,
        motifs_file     = MOTIFS_FILE,
        discovery_root  = REALM_LOCATION,
        atr             = str(ATR),
        rep             = str(REP),
        ddg             = str(DDG),
        extra_params    = EXTRA_PARAMS,
    log:
        os.path.join(TMP_ROOT, "batch_{batch_id}", "rosetta_round2_{sub_id}.log"),
    shell:
        """
        set -e
        BATCH_DIR=$(dirname {output.done_flag})
        SUB_ID={wildcards.sub_id}
        PARAMS_DIR="$BATCH_DIR/round2/$SUB_ID/test_params"

        # Skip sub-batch if no .params files exist in this sub-directory
        if [ ! -d "$PARAMS_DIR" ] || [ -z "$(find "$PARAMS_DIR" -maxdepth 1 -name '*.params' 2>/dev/null)" ]; then
            echo "Sub-batch $SUB_ID has no .params files — skipping"
            touch {output.done_flag}
            exit 0
        fi

        param_count=$(find "$PARAMS_DIR" -maxdepth 1 -name '*.params' | wc -l)
        echo "Sub-batch $SUB_ID: $param_count params file(s) in round2/$SUB_ID/test_params/"

        /home/ji.cheng4-umw/miniforge3/envs/realm_env/bin/python3.11 -c "
import sys, os, glob
sys.path.insert(0, '{params.discovery_root}/function/discovery')
from utils import run_rosetta_discovery_search

any_placements = False
for residue in '{params.anchor_residues}'.split(','):
    res = residue.strip()
    res_dir = os.path.join('$BATCH_DIR', 'round2', '$SUB_ID', res)
    os.makedirs(res_dir, exist_ok=True)
    success = run_rosetta_discovery_search(
        '{params.target_pdb}', res, '{params.motifs_file}',
        '$PARAMS_DIR/', '{params.discovery_root}',
        '{params.atr}', '{params.rep}', '{params.ddg}',
        extra_args_file='{params.extra_params}' if '{params.extra_params}' else None,
        work_dir=res_dir
    )
    if not success:
        # Rosetta may return non-zero even when placements were produced.
        pdbs = glob.glob(os.path.join(res_dir, '*.pdb'))
        tars = glob.glob(os.path.join(res_dir, '*.tar.gz'))
        if pdbs or tars:
            print(f'Rosetta exited non-zero but placements exist for residue {{res}} ({{len(pdbs)}} PDBs, {{len(tars)}} tarballs)')
            any_placements = True
        else:
            print(f'ERROR: Round 2 Rosetta discovery failed for residue {{res}} with no placements')
    else:
        any_placements = True

if not any_placements:
    raise SystemExit(1)
open('{output.done_flag}', 'w').close()
print('Rosetta round 2 sub-batch $SUB_ID complete')
" > {log} 2>&1
        """


# ── Step 6: Score placements from round 2 ───────────────────────────────

rule score_round2:
    """Score placements from round 2 (all sub-batches) and write per-batch scores CSV.
    Sub-batch count is read from the split checkpoint — only actual sub-batches
    are waited on."""
    input:
        lambda wildcards: expand(
            os.path.join(TMP_ROOT, "batch_{batch_id}", "rosetta_round2_{sub_id}.done"),
            batch_id=wildcards.batch_id, sub_id=get_round2_sub_ids(wildcards)),
    output:
        scores_csv  = os.path.join(TMP_ROOT, "batch_{batch_id}", "scores_round2.csv"),
    resources:
        mem_mb=2000,
        cpus=1,
        queue=LSF_QUEUE_DEFAULT,
        walltime=LSF_WALLTIME_DEFAULT,
    params:
        weights_file = WEIGHTS_FILE,
        realm        = REALM_LOCATION,
        batches_dir  = BATCHES_DIR,
    log:
        os.path.join(TMP_ROOT, "batch_{batch_id}", "score_round2.log"),
    shell:
        """
        set -e
        BATCH_DIR=$(dirname {output.scores_csv})

        /home/ji.cheng4-umw/miniforge3/envs/realm_env/bin/python3.11 -c "
import sys, os, csv, glob
sys.path.insert(0, '{params.realm}/function/discovery')
from utils import score_placements

# ── Load batch ligand list for dedup matching ────────────────────────
batch_id = int('{wildcards.batch_id}')
batch_file = os.path.join('{params.batches_dir}', f'{{batch_id // 1000:04d}}', f'batch_{{batch_id}}.txt')
batch_ligands = set()
with open(batch_file, 'r') as fh:
    for line in fh:
        parts = line.strip().split(',')
        if len(parts) >= 1:
            batch_ligands.add(parts[0])

def get_base_ligand(placement_name):
    '''Match placement_name against known batch ligands (longest match first).'''
    for lig in sorted(batch_ligands, key=len, reverse=True):
        if lig in placement_name:
            return lig
    return '_'.join(placement_name.split('_')[:-1])

scores = score_placements(os.path.join('$BATCH_DIR', 'round2'), weights_file='{params.weights_file}')
# Also extract real_motif_ratio from weighted_scores.csv files
# (now two levels deep: round2/<sub_id>/<residue>/weighted_scores.csv)
motif_ratios = {{}}
for csv_path in glob.glob(os.path.join('$BATCH_DIR', 'round2', '*', '*', 'weighted_scores.csv')):
    with open(csv_path, 'r') as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            pdb_name = os.path.basename(row.get('file', ''))
            lig_name = pdb_name.replace('.pdb', '')
            if lig_name:
                try:
                    motif_ratios[lig_name] = float(row.get('real_motif_ratio', 0))
                except ValueError:
                    pass

# ── Deduplicate: keep only best-scoring conformer per base ligand ─
best = {{}}  # base_ligand -> (score, full_ligand_name, motif_ratio)
for lig, sc in scores.items():
    base = get_base_ligand(lig)
    if base not in best or sc > best[base][0]:
        best[base] = (sc, lig, motif_ratios.get(lig, 0.0))

with open('{output.scores_csv}', 'w', newline='') as fh:
    writer = csv.writer(fh)
    writer.writerow(['ligand', 'score', 'real_motif_ratio'])
    for base, (sc, lig, mr) in sorted(best.items(), key=lambda x: -x[1][0]):
        writer.writerow([lig, f'{{sc:.6f}}', f'{{mr:.6f}}'])
print(f'Scored {{len(best)}} unique ligands (from {{len(scores)}} placements) (round 2)')
" > {log} 2>&1

        # ── Delete round2 placement files immediately after scoring ─────────
        # Scores are saved; placement PDBs & tarballs are no longer needed.
        rm -rf "$BATCH_DIR"/round2
        rm -f "$BATCH_DIR"/rosetta_round2_*.log
        """


# ── Step 6.5: Final cleanup — remove all intermediate batch directories ─

rule cleanup_intermediates:
    """Remove heavy intermediate files (placements, PDBs, tarballs, residue dirs)
    while preserving .done checkpoints, score CSVs, test_params/ and _ready
    sentinel files so Snakemake can still track completed steps on re-runs."""
    localrule: True
    input:
        top_ligands = os.path.join(OUTPUT_DIR, "top_ligands.csv"),
    output:
        os.path.join(OUTPUT_DIR, "cleanup.done"),
    params:
        tmp_root = TMP_ROOT,
    resources:
        mem_mb=500,
        cpus=1,
    log:
        os.path.join(OUTPUT_DIR, "cleanup_intermediates.log"),
    shell:
        """
        set -e
        /home/ji.cheng4-umw/miniforge3/envs/realm_env/bin/python3.11 -c "
import shutil, glob, os

for batch_dir in glob.glob(os.path.join('{params.tmp_root}', 'batch_*')):
    for entry in os.listdir(batch_dir):
        epath = os.path.join(batch_dir, entry)
        if os.path.isdir(epath):
            if entry not in ('test_params', 'round1', 'round2'):
                shutil.rmtree(epath, ignore_errors=True)
            elif entry in ('round1', 'round2'):
                # round dirs should already be cleaned by score rules;
                # remove if still present
                shutil.rmtree(epath, ignore_errors=True)
        elif os.path.isfile(epath):
            if not (entry.endswith('.done') or
                    entry.endswith('.csv') or
                    '_ready' in entry):
                os.remove(epath)
    print(f'  Cleaned heavy intermediates in {{os.path.basename(batch_dir)}}')
# Touch done flag
with open('{output}', 'w') as fh:
    fh.write('done\\n')
print('  All heavy intermediate files cleaned; .done, .csv & _ready preserved.')
" > {log} 2>&1
        """


# ── Step 7: Aggregate → top_ligands.csv ─────────────────────────────────

rule aggregate_results:
    """Aggregate all per-batch scores, deduplicate by base ligand, and write top-N CSV."""
    localrule: True
    input:
        # Depend on all score files being ready before aggregation
        expand(os.path.join(TMP_ROOT, "batch_{batch_id}", "scores_round2.csv"),
               batch_id=BATCH_IDS),
    output:
        os.path.join(OUTPUT_DIR, "top_ligands.csv"),
    params:
        max_ligands     = MAX_LIGANDS,
        min_motif_ratio = MIN_MOTIF_RATIO,
        tmp_root        = TMP_ROOT,
    resources:
        mem_mb=16000,
        cpus=1,
    log:
        os.path.join(OUTPUT_DIR, "aggregate_results.log"),
    shell:
        """
        set -e
        /home/ji.cheng4-umw/miniforge3/envs/realm_env/bin/python3.11 -c "
import heapq, csv, glob, os

# Find all round-2 score files produced by upstream rules
scores_files = glob.glob(os.path.join('{params.tmp_root}', 'batch_*', 'scores_round2.csv'))

# ── Build global batch ligand set for dedup matching, plus chunk lookup ─
output_dir = os.path.dirname(os.path.dirname('{params.tmp_root}'))
batches_dir = os.path.join(output_dir, 'batches')
batch_ligands_global = set()
ligand_chunk_info = {{}}  # base_ligand_name -> (chunk, subchunk)
for batch_dir in glob.glob(os.path.join('{params.tmp_root}', 'batch_*')):
    batch_name = os.path.basename(batch_dir)
    batch_id = int(batch_name.split('_')[-1])
    batch_file = os.path.join(batches_dir, f'{{batch_id // 1000:04d}}', f'batch_{{batch_id}}.txt')
    if os.path.exists(batch_file):
        with open(batch_file, 'r') as fh:
            for line in fh:
                parts = line.strip().split(',')
                if len(parts) >= 4:
                    batch_ligands_global.add(parts[0])
                    ligand_chunk_info[parts[0]] = (parts[2], parts[3])

def get_base_ligand_global(placement_name):
    '''Match placement_name against all known batch ligands (longest match first).'''
    for lig in sorted(batch_ligands_global, key=len, reverse=True):
        if lig in placement_name:
            return lig
    return '_'.join(placement_name.split('_')[:-1])

# Collect all scores, deduplicate by base ligand name
seen = {{}}  # base_ligand -> (best_score, full_name, mr, chunk, subchunk)
for sf in scores_files:
    if not os.path.exists(sf):
        continue
    with open(sf, 'r') as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            lig = row['ligand']
            score = float(row['score'])
            mr = float(row['real_motif_ratio'])
            if mr <= {params.min_motif_ratio}:
                continue
            base = get_base_ligand_global(lig)
            chunk, subchunk = ligand_chunk_info.get(base, ('', ''))
            if base not in seen or score > seen[base][0]:
                seen[base] = (score, lig, mr, chunk, subchunk)

# Build top-N heap
heap = []
for base, (score, lig, mr, chunk, subchunk) in seen.items():
    entry = (-score, lig, score, mr, chunk, subchunk)
    if len(heap) < {params.max_ligands}:
        heapq.heappush(heap, entry)
    elif score > -heap[0][0]:
        heapq.heapreplace(heap, entry)

sorted_results = sorted(heap, reverse=True)

with open('{output}', 'w', newline='') as fh:
    writer = csv.writer(fh)
    writer.writerow(['rank', 'ligand_name', 'score', 'real_motif_ratio', 'chunk', 'subchunk'])
    for rank, (neg_score, lig, score, mr, chunk, subchunk) in enumerate(sorted_results, 1):
        writer.writerow([rank, lig, f'{{score:.6f}}', f'{{mr:.6f}}', chunk, subchunk])

print(f'Wrote {{len(sorted_results)}} top ligands to {output}')
" > {log} 2>&1
        """
