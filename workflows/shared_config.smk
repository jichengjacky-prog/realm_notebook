# ═══════════════════════════════════════════════════════════════════════════════
# Shared configuration for all Rosetta Discovery workflow steps.
# Include this AFTER configfile: in each step-specific Snakefile:
#   configfile: "config.yaml"
#   include: "workflows/shared_config.smk"
# ═══════════════════════════════════════════════════════════════════════════════

import os

# ── Paths & parameters (read from config dict, set by configfile: directive) ──
SNAKEFILE_DIR  = workflow.basedir if hasattr(workflow, 'basedir') else os.getcwd()
SHAPEDB_LIST    = config.get("shapedb_list",    "")
TARGET_PDB      = config.get("target_pdb",      "")
ANCHOR_RESIDUES = config.get("anchor_residues", "")
MOTIFS_FILE     = config.get("motifs_file",     "")
OUTPUT_DIR      = config.get("output_dir",      "output")
REALM_LOCATION  = config.get("realm_location",  SNAKEFILE_DIR)
ENAMINE_PATH    = config.get("enamine_path",    "")
PYTHON_BIN      = config.get("python_bin",      "python3.11")
SNAKEMAKE_BIN   = config.get("snakemake_bin",   "snakemake")

# ── Algorithm parameters ─────────────────────────────────────────────────
NUM_CONFORMERS  = config.get("num_conformers",  400)
TOP_N           = config.get("top_n",           10)
BATCH_SIZE      = config.get("batch_size",      10)
MAX_LIGANDS     = config.get("max_ligands",     10)
TOP_HITS        = config.get("top_hits",        0)
ATR             = config.get("atr",             -2.0)
REP             = config.get("rep",             150.0)
DDG             = config.get("ddg",             -9.0)
MIN_MOTIF_RATIO = config.get("min_motif_ratio", 0.25)
MAX_CONFS_PER_SUB = config.get("max_confs_per_sub", 10)

# ── Optional paths ───────────────────────────────────────────────────────
EXTRA_PARAMS    = config.get("extra_params",    "extra_arg/sample_extra_arg")
LICENSE_KEY     = config.get("license_key",     "")
WEIGHTS_FILE    = config.get("weights_file") or os.path.join(REALM_LOCATION, "function", "discovery", "score_weights.json")
TMP_ROOT        = config.get("tmp_root")     or os.path.join(OUTPUT_DIR, "tmp")

os.makedirs(OUTPUT_DIR, exist_ok=True)
os.makedirs(TMP_ROOT, exist_ok=True)

# ── Derived paths ────────────────────────────────────────────────────────
BATCHES_DIR = os.path.join(OUTPUT_DIR, "batches")
BATCH_SHARD_SIZE = 1000
PARSE_DONE_FLAG = os.path.join(OUTPUT_DIR, ".parse_shapedb.done")
os.makedirs(BATCHES_DIR, exist_ok=True)

# ── LSF cluster resource defaults ────────────────────────────────────────
LSF_QUEUE_ROSETTA    = config.get("lsf_queue_rosetta",    "long")
LSF_QUEUE_DEFAULT    = config.get("lsf_queue_default",    "short")
LSF_WALLTIME_ROSETTA  = config.get("lsf_walltime_rosetta",  "24:00")
LSF_WALLTIME_DEFAULT  = config.get("lsf_walltime_default",  "4:00")

# ── Helper: batch file path (sharded into subdirectories) ────────────────
def batch_file_path(batch_id):
    shard = int(batch_id) // BATCH_SHARD_SIZE
    shard_dir = os.path.join(BATCHES_DIR, f"{shard:04d}")
    return os.path.join(shard_dir, f"batch_{batch_id}.txt")

# ── Batch slicing support (for large runs) ───────────────────────────────
# Set via --config batch_start=0 batch_end=500 to process a slice.
# When set, BATCH_IDS is generated from start/end without glob_wildcards.
BATCH_START = int(config.get("batch_start", -1))
BATCH_END   = int(config.get("batch_end",   -1))

# ── Helper: discover or generate batch IDs ────────────────────────────────
def discover_batch_ids(pattern_subdir, pattern_file):
    """
    If batch_start/batch_end are set, generate IDs from the range.
    Otherwise, use glob_wildcards to discover from filesystem.
    This avoids expensive NFS globbing for sliced runs.
    """
    if BATCH_START >= 0 and BATCH_END > BATCH_START:
        return [str(i) for i in range(BATCH_START, BATCH_END)]
    # Fallback: glob from filesystem (slow for 30K dirs on NFS)
    ids, = glob_wildcards(os.path.join(TMP_ROOT, pattern_subdir, pattern_file))
    return sorted(ids, key=int)
