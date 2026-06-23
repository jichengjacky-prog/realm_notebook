# ═══════════════════════════════════════════════════════════════════════════════
# STEP 5: Generate CDPKit conformers for top-N ligands → round2 test_params
#
# Run:  snakemake -s workflows/step5_conformers.smk \
#           --configfile yaml/config_state3.yaml --profile profile/lsf
# ═══════════════════════════════════════════════════════════════════════════════

import os

configfile: "config.yaml"
include: "shared_config.smk"

# ── Discover batches from step 1 output (all batches, but only top-N ligands processed) ─
BATCH_IDS = discover_batch_ids("batch_{batch_id}", "extract_params.done")

TOP_LIST = os.path.join(OUTPUT_DIR, "top_ligands_round1.txt")
STEP5_DONE = os.path.join(OUTPUT_DIR, ".step5_conformers.done")


# ── Rules ─────────────────────────────────────────────────────────────────

rule all:
    input: STEP5_DONE


rule generate_conformers:
    """Re-extract top-N ligand params and generate CDPKit conformers.
    Submits one LSF job per ligand internally via bsub, then waits for all."""
    input:
        batch_file  = lambda wildcards: batch_file_path(wildcards.batch_id),
        top_list    = TOP_LIST,
    output:
        done_flag        = os.path.join(TMP_ROOT, "batch_{batch_id}", "conformers_generated.done"),
        test_params_done = os.path.join(TMP_ROOT, "batch_{batch_id}", ".round2_ready"),
    resources:
        load=1,
        mem_mb=2000,
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
        mkdir -p $BATCH_DIR
        rm -rf "$BATCH_DIR"/round1

        LIGANDS_FILE="$BATCH_DIR/.r2_ligands.txt"
        /home/ji.cheng4-umw/miniforge3/envs/realm_env/bin/python3.11 -c "
import sys, os
sys.path.insert(0, '{params.realm_location}/function/discovery')

top_set = set()
with open('{input.top_list}') as fh:
    fh.readline()
    for line in fh:
        line = line.strip()
        if line:
            top_set.add(line.split(',')[0])

ligands = []
with open('{input.batch_file}') as fh:
    for line in fh:
        parts = line.strip().split(',')
        if len(parts) == 4:
            name = parts[0]
            if any(name in entry for entry in top_set):
                ligands.append((name, int(parts[1]), parts[2], parts[3]))

with open('$LIGANDS_FILE', 'w') as f:
    for l in ligands:
        f.write(','.join(str(x) for x in l) + '\\n')
print(f'Top-N ligands in this batch: {{len(ligands)}}')
" 2>&1 | tee -a {log}

        if [ ! -s "$LIGANDS_FILE" ]; then
            echo 'No top-N ligands in this batch, skipping'
            mkdir -p "$BATCH_DIR"
            touch {output.test_params_done}
            touch {output.done_flag}
            exit 0
        fi

        JOB_IDS=""
        while IFS=, read -r lig_name conf_num chunk subchunk || [ -n "$lig_name" ]; do
            [ -n "$lig_name" ] || continue
            JOB_NAME="smk_genconf_${{lig_name:0:20}}"
            SCRIPT="$BATCH_DIR/genconf_$lig_name.py"
            cat > "$SCRIPT" << 'PYEOF'
import sys, os
lig_name = sys.argv[1]
conf_num = int(sys.argv[2])
chunk = sys.argv[3]
subchunk = sys.argv[4]
batch_dir = sys.argv[5]
realm = sys.argv[6]
enamine = sys.argv[7]
num_conf = int(sys.argv[8])

sys.path.insert(0, os.path.join(realm, 'function/discovery'))
from utils import create_test_params_dir, generate_and_add_conformers_to_test_params

create_test_params_dir(
    [(lig_name, conf_num, chunk, subchunk)], batch_dir, realm, enamine,
    round_name='round2'
)
generate_and_add_conformers_to_test_params(
    batch_dir, [(lig_name, conf_num, chunk, subchunk)], realm, enamine,
    num_conformers=num_conf, round_name='round2'
)
print(f'Conformers generated for {{lig_name}}')
PYEOF
            JOB_ID=$(bsub -q {resources.queue} \
                -W {resources.walltime} \
                -M 4000 \
                -n 1 \
                -R 'span[hosts=1] rusage[mem=4000]' \
                -J "$JOB_NAME" \
                -o "$BATCH_DIR/genconf_$lig_name.out" \
                -e "$BATCH_DIR/genconf_$lig_name.err" \
                /home/ji.cheng4-umw/miniforge3/envs/realm_env/bin/python3.11 "$SCRIPT" \
                "$lig_name" "$conf_num" "$chunk" "$subchunk" \
                "$BATCH_DIR" '{params.realm_location}' '{params.enamine_path}' \
                '{params.num_conformers}' 2>&1 | grep -oP '<\d+>' | tr -d '<>' || true)
            [ -n "$JOB_ID" ] && JOB_IDS="$JOB_IDS $JOB_ID"
            echo "Submitted $JOB_NAME ($JOB_ID)"
            sleep 0.1
        done < "$LIGANDS_FILE"

        if [ -z "$JOB_IDS" ]; then
            echo "No conformer jobs submitted — skipping"
            touch {output.test_params_done}
            touch {output.done_flag}
            exit 0
        fi

        echo "Submitted jobs:$JOB_IDS — waiting for .params files..."
        EXPECTED_LIGS=$(wc -l < "$LIGANDS_FILE")
        echo "Expecting params for $EXPECTED_LIGS ligand(s)"

        while true; do
            PARAMS_COUNT=$(find "$BATCH_DIR"/round2 -mindepth 3 -name '*.params' 2>/dev/null | wc -l)
            echo "  Params progress: $PARAMS_COUNT .params files (need at least $EXPECTED_LIGS)"
            if [ "$PARAMS_COUNT" -ge "$EXPECTED_LIGS" ]; then
                echo "All .params files present — conformer generation complete"
                break
            fi
            STILL_RUNNING=0
            for JID in $JOB_IDS; do
                STAT=$(bjobs -o stat -noheader "$JID" 2>/dev/null | tr -d ' ')
                if [ -n "$STAT" ]; then
                    STILL_RUNNING=1
                    break
                fi
            done
            if [ "$STILL_RUNNING" -eq 0 ] && [ "$PARAMS_COUNT" -lt "$EXPECTED_LIGS" ]; then
                echo "ERROR: All LSF jobs finished but only $PARAMS_COUNT params found"
                exit 1
            fi
            sleep 30
        done

        touch {output.test_params_done}
        touch {output.done_flag}
        echo 'Conformer generation complete'
        """


rule step5_sentinel:
    """Aggregate sentinel: depends on all conformer_generated completing."""
    localrule: True
    input:
        expand(os.path.join(TMP_ROOT, "batch_{batch_id}", "conformers_generated.done"),
               batch_id=BATCH_IDS),
    output:
        STEP5_DONE,
    resources:
        mem_mb=500,
        cpus=1,
    shell:
        """
        touch {output}
        echo "Step 5 (generate_conformers) complete"
        """
