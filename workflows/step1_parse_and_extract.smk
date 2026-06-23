# ═══════════════════════════════════════════════════════════════════════════════
# STEP 1: Parse shapedb list → batch files + extract library params → test_params
#
# Run:  snakemake -s workflows/step1_parse_and_extract.smk \
#           --configfile yaml/config_state3.yaml --profile profile/lsf
# ═══════════════════════════════════════════════════════════════════════════════

import os

configfile: "config.yaml"
include: "shared_config.smk"


# ── Parse shapedb list → batch files (idempotent, uses done flag) ────────

def _read_batch_count_from_flag():
    with open(PARSE_DONE_FLAG) as f:
        first_line = f.readline().strip()
    if not first_line:
        return None
    tokens = first_line.split()
    if not tokens:
        return None
    try:
        return int(tokens[0])
    except (ValueError, IndexError):
        return None


should_parse = True
if os.path.exists(PARSE_DONE_FLAG):
    n = _read_batch_count_from_flag()
    if n is not None and n > 0:
        batch_ids = list(range(n))
        should_parse = False
    else:
        os.remove(PARSE_DONE_FLAG)

if should_parse:
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

    with open(PARSE_DONE_FLAG, "w") as f:
        f.write(f"{len(batch_ids)} batches, {len(entries)} entries\n")

BATCH_IDS = discover_batch_ids("batch_{batch_id}", "extract_params.done")

# ── Step sentinel ────────────────────────────────────────────────────────
STEP1_DONE = os.path.join(OUTPUT_DIR, ".step1_extract.done")


# ── Rules ─────────────────────────────────────────────────────────────────

rule all:
    input: STEP1_DONE


rule extract_params:
    """Extract conformer params from the Enamine library for one batch."""
    input:
        batch_file = lambda wildcards: batch_file_path(wildcards.batch_id),
    output:
        done_flag        = os.path.join(TMP_ROOT, "batch_{batch_id}", "extract_params.done"),
        test_params_done = os.path.join(TMP_ROOT, "batch_{batch_id}", ".round1_ready"),
    params:
        realm_location = REALM_LOCATION,
        enamine_path   = ENAMINE_PATH,
        tmp_root       = TMP_ROOT,
    resources:
        load=1,
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

lig_dirs = create_test_params_dir(ligands, '$BATCH_DIR', '{params.realm_location}',
                       '{params.enamine_path}', tmp_root='{params.tmp_root}', round_name='round1')
print(f'extract_params complete — {{len(lig_dirs)}} ligand dirs created')
" > {log} 2>&1
        touch {output.test_params_done}
        touch {output.done_flag}
        """


rule step1_sentinel:
    """Aggregate sentinel: depends on all extract_params completing."""
    localrule: True
    input:
        expand(os.path.join(TMP_ROOT, "batch_{batch_id}", "extract_params.done"),
               batch_id=BATCH_IDS),
    output:
        STEP1_DONE,
    resources:
        mem_mb=500,
        cpus=1,
    shell:
        """
        touch {output}
        echo "Step 1 (extract_params) complete"
        """
