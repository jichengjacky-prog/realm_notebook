# ═══════════════════════════════════════════════════════════════════════════════
# STEP 2: Rosetta discovery round 1 — per-ligand LSF jobs via bsub
#
# Run:  snakemake -s workflows/step2_rosetta1.smk \
#           --configfile yaml/config_state3.yaml --profile profile/lsf \
#           --resources load=150
# ═══════════════════════════════════════════════════════════════════════════════

import os

configfile: "config.yaml"
include: "shared_config.smk"

# ── Discover batches from step 1 output ──────────────────────────────────
BATCH_IDS = discover_batch_ids("batch_{batch_id}", "extract_params.done")

# ── Step sentinel (slice-aware) ──────────────────────────────────────────
SLICE_TAG = f"slice_{BATCH_START}_{BATCH_END}" if BATCH_START >= 0 else "all"
STEP2_DONE = os.path.join(OUTPUT_DIR, f".step2_rosetta1_{SLICE_TAG}.done")


# ── Rules ─────────────────────────────────────────────────────────────────

rule all:
    input: STEP2_DONE


rule rosetta_discovery_round1:
    """Run Rosetta ligand discovery search for one batch.
    Submits one LSF job per ligand internally via bsub, then waits for all."""
    input:
        params_done      = os.path.join(TMP_ROOT, "batch_{batch_id}", "extract_params.done"),
        test_params_done = os.path.join(TMP_ROOT, "batch_{batch_id}", ".round1_ready"),
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
        rosetta_queue   = LSF_QUEUE_ROSETTA,
        default_queue   = LSF_QUEUE_DEFAULT,
        rosetta_walltime = LSF_WALLTIME_ROSETTA,
        python_bin      = PYTHON_BIN,
    resources:
        load=10,         
        mem_mb=8000,
        cpus=1,
        queue=LSF_QUEUE_ROSETTA,
        walltime=LSF_WALLTIME_ROSETTA,
    log:
        os.path.join(TMP_ROOT, "batch_{batch_id}", "rosetta_round1.log"),
    shell:
        """
        # NOTE: set +e (not set -e) because some per-ligand Rosetta jobs
        # may exit non-zero (no placements found) which is normal.
        # Explicit error checks are used where needed.
        BATCH_DIR=$(dirname {output.done_flag})
        mkdir -p $BATCH_DIR
        echo "Starting Rosetta round 1 for batch in $BATCH_DIR"

        # ── Submit one LSF job per ligand in parallel ──────────────────
        JOB_IDS=()
        LIG_NAMES=()
        ANCHOR_COUNT=$(echo '{params.anchor_residues}' | tr ',' '\n' | sed '/^[[:space:]]*$/d' | wc -l)
        for TP_DIR in "$BATCH_DIR"/round1/*/test_params/; do
            [ -d "$TP_DIR" ] || continue
            LIG_NAME=$(basename $(dirname "$TP_DIR"))
            # Skip if ligand already has all expected CSVs (idempotent across restarts)
            EXISTING_CSV_COUNT=$(find "$BATCH_DIR/round1/$LIG_NAME" -name 'weighted_scores.csv' 2>/dev/null | wc -l)
            if [ "$EXISTING_CSV_COUNT" -ge "$ANCHOR_COUNT" ] && [ "$ANCHOR_COUNT" -gt 0 ]; then
                echo "Ligand $LIG_NAME: already complete ($EXISTING_CSV_COUNT CSVs), skipping submission"
                rm -f "$BATCH_DIR/round1/ros1_$LIG_NAME.out" \
                      "$BATCH_DIR/round1/ros1_$LIG_NAME.err" \
                      "$BATCH_DIR/round1/ros1_$LIG_NAME.py"
                touch "$BATCH_DIR/round1/$LIG_NAME.done"
                continue
            fi
            LIG_NAMES+=("$LIG_NAME")
            JOB_NAME="smk_ros1_${{LIG_NAME:0:20}}"
            # Write per-ligand Python script to round1/
            SCRIPT="$BATCH_DIR/round1/ros1_$LIG_NAME.py"
            cat > "$SCRIPT" << 'PYEOF'
import sys, os, glob
sys.path.insert(0, '{params.discovery_root}/function/discovery')
from utils import run_rosetta_discovery_search

lig_name = sys.argv[1]
tp_dir = sys.argv[2]
batch_dir = sys.argv[3]
any_placements = False
for residue in '{params.anchor_residues}'.split(','):
    res = residue.strip()
    res_dir = os.path.join(batch_dir, 'round1', lig_name, res)
    os.makedirs(res_dir, exist_ok=True)
    # Idempotency guard: skip if weighted_scores.csv already exists for this residue
    if os.path.isfile(os.path.join(res_dir, 'weighted_scores.csv')):
        print(f'  SKIP: weighted_scores.csv already exists for {{lig_name}}/{{res}}')
        continue
    success = run_rosetta_discovery_search(
        '{params.target_pdb}', res, '{params.motifs_file}',
        tp_dir, '{params.discovery_root}',
        '{params.atr}', '{params.rep}', '{params.ddg}',
        extra_args_file='{params.extra_params}' or None,
        work_dir=res_dir
    )
    if not success:
        pdbs = glob.glob(os.path.join(res_dir, '*.pdb'))
        tars = glob.glob(os.path.join(res_dir, '*.tar.gz'))
        if not pdbs and not tars:
            print(f'WARNING: Rosetta produced no output for {{lig_name}}/{{res}}', file=sys.stderr)
        else:
            print(f'Rosetta: placements exist for {{lig_name}}/{{res}} despite non-zero exit')
    # Remove heavy files, keep only CSVs
    for f in glob.glob(os.path.join(res_dir, '*')):
        if not f.endswith('.csv'):
            try:
                if os.path.isdir(f):
                    import shutil; shutil.rmtree(f)
                else:
                    os.remove(f)
            except Exception:
                pass
print('Rosetta done for ' + lig_name)
PYEOF
            JOB_ID=$(bsub \
                -W {params.rosetta_walltime} \
                -q {params.rosetta_queue} \
                -M 8000 \
                -n 1 \
                -R 'span[hosts=1] rusage[mem=8000]' \
                -J "$JOB_NAME" \
                -o "$BATCH_DIR/round1/ros1_$LIG_NAME.out" \
                -e "$BATCH_DIR/round1/ros1_$LIG_NAME.err" \
                {params.python_bin} "$SCRIPT" "$LIG_NAME" "$TP_DIR" "$BATCH_DIR" 2>&1 | grep -oP '<\d+>' | tr -d '<>' || true)
            [ -n "$JOB_ID" ] && JOB_IDS+=("$JOB_ID")
            echo "Submitted $JOB_NAME ($JOB_ID)"
            sleep 0.05
        done

        if [ ${{#JOB_IDS[@]}} -eq 0 ]; then
            echo "No ligands found — nothing to run"
            touch {output.done_flag}
            exit 0
        fi

        # ── Wait for all per-ligand jobs with stuck-job watchdog ─────
        # Two-layer protection:
        #   1. CPU monitor: 0% CPU peak for STUCK_THRESHOLD checks → bkill zombie
        #   2. Hard timeout: MAX_WAIT_HOURS overall → break out (pipeline can resume)
        echo "Waiting for all per-ligand jobs to finish..."
        WAIT_START=$(date +%s)
        STUCK_THRESHOLD=4          # 4 consecutive checks (2h) with 0% CPU → stuck
        MAX_WAIT_HOURS=48          # absolute max wait time for the batch
        declare -A JOB_STUCK_COUNT

        while true; do
            NOW=$(date +%s)
            ELAPSED_HRS=$(( (NOW - WAIT_START) / 3600 ))
            if [ "$ELAPSED_HRS" -ge "$MAX_WAIT_HOURS" ]; then
                echo "WARNING: Maximum wait time ($MAX_WAIT_HOURS h) exceeded — breaking out"
                break
            fi

            STILL_RUNNING=0
            for JID in "${{JOB_IDS[@]}}"; do
                STAT=$(bjobs -o stat -noheader "$JID" 2>/dev/null | tr -d ' ')
                if [ "$STAT" = "RUN" ] || [ "$STAT" = "PEND" ]; then
                    STILL_RUNNING=1

                    # ── Stuck-job detection: 0% CPU peak = hung on NFS or idle ──
                    if [ "$STAT" = "RUN" ]; then
                        CPU_PEAK=$(bjobs -l "$JID" 2>/dev/null | grep "CPU PEAK:" | head -1 | awk -F'[:;]' '{{print $2}}' | tr -d ' ')
                        if [ "$CPU_PEAK" = "0.00" ] || [ -z "$CPU_PEAK" ]; then
                            JOB_STUCK_COUNT[$JID]=$((${{JOB_STUCK_COUNT[$JID]:-0}} + 1))
                            if [ ${{JOB_STUCK_COUNT[$JID]}} -ge "$STUCK_THRESHOLD" ]; then
                                echo "WARNING: Job $JID appears stuck (0% CPU for $STUCK_THRESHOLD checks, ~$((STUCK_THRESHOLD * 5 / 60))h) — killing"
                                bkill "$JID" 2>/dev/null || true
                                unset JOB_STUCK_COUNT[$JID]
                            fi
                        else
                            JOB_STUCK_COUNT[$JID]=0
                        fi
                    fi
                fi
            done

            if [ "$STILL_RUNNING" -eq 0 ]; then
                CSV_COUNT=$(find "$BATCH_DIR"/round1 -name 'weighted_scores.csv' 2>/dev/null | wc -l)
                echo "All per-ligand jobs finished — $CSV_COUNT weighted_scores.csv files produced"
                break
            fi

            # Periodic status report (every ~10 min)
            if [ $(( (NOW - WAIT_START) % 600 )) -lt 30 ]; then
                RUN_COUNT=0; PEND_COUNT=0
                for JID in "${{JOB_IDS[@]}}"; do
                    S=$(bjobs -o stat -noheader "$JID" 2>/dev/null | tr -d ' ')
                    [ "$S" = "RUN" ] && RUN_COUNT=$((RUN_COUNT + 1))
                    [ "$S" = "PEND" ] && PEND_COUNT=$((PEND_COUNT + 1))
                done
                CSV_COUNT=$(find "$BATCH_DIR"/round1 -name 'weighted_scores.csv' 2>/dev/null | wc -l)
                echo "[+${{ELAPSED_HRS}}h] Status: $RUN_COUNT RUN, $PEND_COUNT PEND, $CSV_COUNT CSVs so far"
            fi

            sleep 30
        done

        # ── Per-ligand cleanup: remove logs if ALL CSVs present, touch done ──
        for i in "${{!LIG_NAMES[@]}}"; do
            LIG="${{LIG_NAMES[$i]}}"
            LIG_DIR="$BATCH_DIR/round1/$LIG"
            CSV_COUNT=$(find "$LIG_DIR" -name 'weighted_scores.csv' 2>/dev/null | wc -l)
            if [ "$CSV_COUNT" -ge "$ANCHOR_COUNT" ] && [ "$ANCHOR_COUNT" -gt 0 ]; then
                rm -f "$BATCH_DIR/round1/ros1_$LIG.out" \
                      "$BATCH_DIR/round1/ros1_$LIG.err" \
                      "$BATCH_DIR/round1/ros1_$LIG.py"
                touch "$BATCH_DIR/round1/$LIG.done"
                echo "Ligand $LIG: SUCCESS ($CSV_COUNT/$ANCHOR_COUNT CSVs) — cleaned logs, touched done"
            elif [ "$CSV_COUNT" -gt 0 ]; then
                echo "Ligand $LIG: PARTIAL ($CSV_COUNT/$ANCHOR_COUNT CSVs) — keeping logs for re-run"
            else
                echo "Ligand $LIG: FAILED or no output — keeping logs for debugging"
            fi
        done

        touch {output.done_flag}
        echo 'Rosetta round 1 complete'
        exit 0
        """


rule step2_sentinel:
    """Aggregate sentinel: depends on all rosetta_round1 completing."""
    localrule: True
    input:
        expand(os.path.join(TMP_ROOT, "batch_{batch_id}", "rosetta_round1.done"),
               batch_id=BATCH_IDS),
    output:
        STEP2_DONE,
    resources:
        mem_mb=500,
        cpus=1,
    shell:
        """
        touch {output}
        echo "Step 2 (rosetta_round1) complete"
        """
