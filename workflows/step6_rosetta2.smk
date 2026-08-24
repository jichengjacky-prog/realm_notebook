# ═══════════════════════════════════════════════════════════════════════════════
# STEP 6: Rosetta discovery round 2 — CDPKit conformers, per-ligand LSF jobs
#
# Run:  snakemake -s workflows/step6_rosetta2.smk \
#           --configfile yaml/config_state3.yaml --profile profile/lsf \
#           --resources load=150
# ═══════════════════════════════════════════════════════════════════════════════

import os

configfile: "config.yaml"
include: "shared_config.smk"

# ── Discover batches from step 5 output ──────────────────────────────────
BATCH_IDS = discover_batch_ids("batch_{batch_id}", "conformers_generated.done")

STEP6_DONE = os.path.join(OUTPUT_DIR, ".step6_rosetta2.done")


# ── Rules ─────────────────────────────────────────────────────────────────

rule all:
    input: STEP6_DONE


rule rosetta_discovery_round2:
    """Run Rosetta discovery with CDPKit-generated conformers.
    One bsub job per ligand under round2/."""
    input:
        confs_done       = os.path.join(TMP_ROOT, "batch_{batch_id}", "conformers_generated.done"),
        test_params_done = os.path.join(TMP_ROOT, "batch_{batch_id}", ".round2_ready"),
    output:
        done_flag   = os.path.join(TMP_ROOT, "batch_{batch_id}", "rosetta_round2.done"),
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
        sif_rosetta     = SIF_ROSETTA,
    resources:
        load=10,
        mem_mb=500,
        cpus=1,
        # Controller job: submits the arrays, then polls up to MAX_WAIT_HOURS.
        # Stays on the long queue; the serial Rosetta arrays use both queues.
        queue=LSF_QUEUE_MONITOR,
        walltime="168:00",
    log:
        os.path.join(TMP_ROOT, "batch_{batch_id}", "rosetta_round2.log"),
    shell:
        """
        # NOTE: set +e (not set -e) because some per-ligand Rosetta jobs
        # may exit non-zero (no placements found) which is normal.
        # Explicit error checks are used where needed.
        set +e

        # ── Trap: on ANY exit, reconcile per-ligand .done markers so a
        # ligand whose raw_scores.csv files are all present is always marked
        # done (even on early error exits), then log and re-exit. The batch
        # done flag is NEVER touched here.
        _cleanup_and_exit() {{
            local exit_code=$?
            if [ -n "$BATCH_DIR" ] && [ -n "$ROUND" ] && [ "${{ANCHOR_COUNT:-0}}" -gt 0 ]; then
                for LIG_DIR in "$BATCH_DIR/$ROUND"/*/; do
                    [ -d "$LIG_DIR" ] || continue
                    LIG_NAME=$(basename "$LIG_DIR")
                    CSV_COUNT=$(find "$LIG_DIR" -name 'raw_scores.csv' 2>/dev/null | wc -l)
                    if [ "$CSV_COUNT" -ge "$ANCHOR_COUNT" ]; then
                        touch "$BATCH_DIR/$ROUND/$LIG_NAME.done"
                    fi
                done
            fi
            echo "[$(date)] Shell script exiting (trap) with captured code $exit_code"
            exit "$exit_code"
        }}
        trap '_cleanup_and_exit' EXIT

        BATCH_DIR=$(dirname {output.done_flag})
        mkdir -p $BATCH_DIR
        echo "Starting Rosetta round 2 for batch in $BATCH_DIR"

        # Remove any stale done flag left by a previous (possibly false) run.
        # It is re-created below only after the CSV count has been verified,
        # so a failed re-run can no longer leave the batch marked as done.
        rm -f {output.done_flag}

        ROUND="round2"
        ANCHOR_COUNT=$(echo '{params.anchor_residues}' | tr ',' '\n' | sed '/^[[:space:]]*$/d' | wc -l)

        # ── Collect pending ligands (skip already-complete) ──────────
        PENDING_LIGS=()
        for TP_DIR in "$BATCH_DIR/$ROUND"/*/test_params/; do
            [ -d "$TP_DIR" ] || continue
            LIG_NAME=$(basename $(dirname "$TP_DIR"))
            EXISTING_CSV_COUNT=$(find "$BATCH_DIR/$ROUND/$LIG_NAME" -name 'raw_scores.csv' 2>/dev/null | wc -l)
            if [ "$EXISTING_CSV_COUNT" -ge "$ANCHOR_COUNT" ] && [ "$ANCHOR_COUNT" -gt 0 ]; then
                echo "Ligand $LIG_NAME: already complete ($EXISTING_CSV_COUNT CSVs), skipping submission"
                rm -f "$BATCH_DIR/$ROUND/ros2_$LIG_NAME.out" \
                      "$BATCH_DIR/$ROUND/ros2_$LIG_NAME.err" \
                      "$BATCH_DIR/$ROUND/ros2_$LIG_NAME.py"
                touch "$BATCH_DIR/$ROUND/$LIG_NAME.done"
                continue
            fi
            PENDING_LIGS+=("$LIG_NAME")
            # Remove any stale completion marker from an earlier attempt so
            # the final check only sees chunks that finished THIS run.
            rm -f "$BATCH_DIR/$ROUND/$LIG_NAME.done_pre_exit"
        done

        if [ ${{#PENDING_LIGS[@]}} -eq 0 ]; then
            echo "No pending ligands found — nothing to run"
            touch {output.done_flag}
            exit 0
        fi

        echo "Total pending ligands: ${{#PENDING_LIGS[@]}}"

        # ── Bundle into chunks of 10, submit as LSF job arrays ─────
        CHUNK_SIZE=10
        TOTAL_LIGS=${{#PENDING_LIGS[@]}}
        NUM_CHUNKS=$(( (TOTAL_LIGS + CHUNK_SIZE - 1) / CHUNK_SIZE ))
        ARRAY_JOB_IDS=()

        for ((chunk=0; chunk<NUM_CHUNKS; chunk++)); do
            START=$((chunk * CHUNK_SIZE))
            CHUNK_LIGS=("${{PENDING_LIGS[@]:$START:$CHUNK_SIZE}}")
            N=${{#CHUNK_LIGS[@]}}

            # ── Queue check: wait if >5000 jobs already queued ──
            # Rosetta jobs may land on either queue; count both.
            while true; do
                QUEUE_COUNT=0
                for Q in {params.rosetta_queue}; do
                    QUEUE_COUNT=$(( QUEUE_COUNT + $(bjobs -q "$Q" -noheader 2>/dev/null | wc -l) ))
                done
                if [ "$QUEUE_COUNT" -le 5000 ]; then
                    break
                fi
                echo "[$(date)] Queues '{params.rosetta_queue}' have $QUEUE_COUNT jobs (>5000), waiting 10 min before submitting chunk $((chunk+1))/$NUM_CHUNKS..."
                sleep 600
            done

            # Build space-separated ligand list for script args
            LIG_LIST="${{CHUNK_LIGS[@]}}"

            CHUNK_SCRIPT="$BATCH_DIR/$ROUND/ros2_chunk_${{chunk}}.py"
            rm -f "$CHUNK_SCRIPT"
            cat > "$CHUNK_SCRIPT" << 'PYEOF'
import sys, os, glob, atexit

atexit.register(sys.stdout.flush)
atexit.register(sys.stderr.flush)

# ── Determine which ligand to process from LSF job array index ──────
job_index = int(os.environ.get('LSB_JOBINDEX', '1')) - 1  # LSF is 1-based
batch_dir = sys.argv[1]
round_name = sys.argv[2]
all_ligands = sys.argv[3:]

if job_index >= len(all_ligands):
    print(f'Job index {{job_index+1}} out of range (have {{len(all_ligands)}} ligands), exiting cleanly')
    sys.stdout.flush()
    os._exit(0)

lig_name = all_ligands[job_index]
tp_dir = os.path.join(batch_dir, round_name, lig_name, 'test_params')
print(f'[{{round_name}}] Processing ligand {{job_index+1}}/{{len(all_ligands)}}: {{lig_name}}')

sys.path.insert(0, '{params.discovery_root}/function/discovery')
from utils import run_rosetta_discovery_search

for residue in '{params.anchor_residues}'.split(','):
    res = residue.strip()
    res_dir = os.path.join(batch_dir, round_name, lig_name, res)
    os.makedirs(res_dir, exist_ok=True)
    if os.path.isfile(os.path.join(res_dir, 'raw_scores.csv')):
        print(f'  SKIP: raw_scores.csv already exists for {{lig_name}}/{{res}}')
        continue
    success = run_rosetta_discovery_search(
        '{params.target_pdb}', res, '{params.motifs_file}',
        tp_dir, '{params.discovery_root}',
        '{params.atr}', '{params.rep}', '{params.ddg}',
        extra_args_file='{params.extra_params}' or None,
        work_dir=res_dir,
        rosetta_sif='{params.sif_rosetta}'
    )
    if not success:
        pdbs = glob.glob(os.path.join(res_dir, '*.pdb'))
        tars = glob.glob(os.path.join(res_dir, '*.tar.gz'))
        if not pdbs and not tars:
            print(f'WARNING: Rosetta produced no output for {{lig_name}}/{{res}}', file=sys.stderr)
        else:
            print(f'Rosetta: placements exist for {{lig_name}}/{{res}} despite non-zero exit')
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
sys.stdout.flush()
sys.stderr.flush()
open(os.path.join(batch_dir, round_name, lig_name + '.done_pre_exit'), 'w').close()
os._exit(0)
PYEOF

            # ── Controller decides this chunk's queue: alternate between
            # the queues in {params.rosetta_queue} ("long short") so roughly
            # half the Rosetta work lands on each queue.
            QLIST=({params.rosetta_queue})
            CHUNK_Q="${{QLIST[$(( chunk % ${{#QLIST[@]}} ))]}}"

            ARRAY_JOB_NAME="smk_ros2_c${{chunk}}"
            ARRAY_JOB_ID=$(bsub \
                -W {params.rosetta_walltime} \
                -q "$CHUNK_Q" \
                -M 6000 \
                -n 1 \
                -R 'span[hosts=1] rusage[mem=6000]' \
                -J "${{ARRAY_JOB_NAME}}[1-${{N}}]" \
                -o "$BATCH_DIR/$ROUND/ros2_chunk_${{chunk}}_%I.out" \
                -e "$BATCH_DIR/$ROUND/ros2_chunk_${{chunk}}_%I.err" \
                {params.python_bin} "$CHUNK_SCRIPT" "$BATCH_DIR" "$ROUND" $LIG_LIST \
                2>&1 | grep -oP '<\d+>' | tr -d '<>' || true)
            [ -n "$ARRAY_JOB_ID" ] && ARRAY_JOB_IDS+=("$ARRAY_JOB_ID")
            echo "Submitted chunk $((chunk+1))/$NUM_CHUNKS: ${{ARRAY_JOB_NAME}}[1-${{N}}] ($ARRAY_JOB_ID) on queue $CHUNK_Q — ${{N}} ligands"
            sleep 1
        done

        if [ ${{#ARRAY_JOB_IDS[@]}} -eq 0 ]; then
            echo "ERROR: No array jobs were submitted — NOT marking batch done"
            exit 1
        fi

        # ── Wait for all array jobs ────────────────────────────────
        echo "Waiting for ${{#ARRAY_JOB_IDS[@]}} array jobs to finish..."
        WAIT_START=$(date +%s)
        MAX_WAIT_HOURS=48
        ALL_FINISHED=0

        while true; do
            NOW=$(date +%s)
            ELAPSED_HRS=$(( (NOW - WAIT_START) / 3600 ))
            if [ "$ELAPSED_HRS" -ge "$MAX_WAIT_HOURS" ]; then
                echo "WARNING: Maximum wait time ($MAX_WAIT_HOURS h) exceeded — breaking out"
                break
            fi

            STILL_RUNNING=0
            for ARRAY_JID in "${{ARRAY_JOB_IDS[@]}}"; do
                # NOTE: query by positional job id (bjobs <id>), NOT -J.
                # On this cluster `bjobs -J <numeric-id>` matches job NAMES
                # only and silently returns nothing. Only RUN/PEND/... states
                # keep the loop alive; DONE elements are not running.
                # No stuck-array watchdog: LSF kills the chunk at its 8h
                # walltime, so an idle array cannot hang forever, and the old
                # watchdog's element check compared a multi-line status
                # against "RUN" and bkilled every healthy array.
                ELEM_STATS=$(bjobs -o "stat" -noheader "$ARRAY_JID" 2>/dev/null | tr -d ' ')
                if echo "$ELEM_STATS" | grep -qE 'RUN|PEND|WAIT|USUSP|PSUSP|SSUSP'; then
                    STILL_RUNNING=1
                fi
            done

            if [ "$STILL_RUNNING" -eq 0 ]; then
                ALL_FINISHED=1
                CSV_COUNT=$(find "$BATCH_DIR/$ROUND" -name 'raw_scores.csv' 2>/dev/null | wc -l)
                echo "All array jobs finished — $CSV_COUNT raw_scores.csv files produced"
                break
            fi

            if [ $(( (NOW - WAIT_START) % 600 )) -lt 30 ]; then
                TOT_RUN=0; TOT_PEND=0
                for ARRAY_JID in "${{ARRAY_JOB_IDS[@]}}"; do
                    ES=$(bjobs -o "stat" -noheader "$ARRAY_JID" 2>/dev/null | tr -d ' ')
                    TOT_RUN=$((TOT_RUN + $(echo "$ES" | grep -c "RUN" || echo 0)))
                    TOT_PEND=$((TOT_PEND + $(echo "$ES" | grep -c "PEND" || echo 0)))
                done
                CSV_COUNT=$(find "$BATCH_DIR/$ROUND" -name 'raw_scores.csv' 2>/dev/null | wc -l)
                echo "[+${{ELAPSED_HRS}}h] Status: $TOT_RUN RUN, $TOT_PEND PEND, $CSV_COUNT CSVs so far"
            fi

            sleep 30
        done

        # ── All arrays must have finished before success is possible ─
        if [ "$ALL_FINISHED" -ne 1 ]; then
            echo "ERROR: Wait loop ended without all array jobs finishing — NOT marking batch done"
            exit 1
        fi

        # ── Per-ligand cleanup ────────────────────────────────────
        for LIG in "${{PENDING_LIGS[@]}}"; do
            LIG_DIR="$BATCH_DIR/$ROUND/$LIG"
            CSV_COUNT=$(find "$LIG_DIR" -name 'raw_scores.csv' 2>/dev/null | wc -l)
            if [ "$CSV_COUNT" -ge "$ANCHOR_COUNT" ] && [ "$ANCHOR_COUNT" -gt 0 ]; then
                touch "$BATCH_DIR/$ROUND/$LIG.done"
                echo "Ligand $LIG: SUCCESS ($CSV_COUNT/$ANCHOR_COUNT CSVs)"
            elif [ "$CSV_COUNT" -gt 0 ]; then
                echo "Ligand $LIG: PARTIAL ($CSV_COUNT/$ANCHOR_COUNT CSVs)"
            else
                echo "Ligand $LIG: FAILED or no output"
            fi
        done

        # ── Remove chunk-level artifacts ─────────────────────────
        rm -f "$BATCH_DIR/$ROUND"/ros2_chunk_*.py \
              "$BATCH_DIR/$ROUND"/ros2_chunk_*.out \
              "$BATCH_DIR/$ROUND"/ros2_chunk_*.err \
              "$BATCH_DIR/$ROUND"/ros2_*.out \
              "$BATCH_DIR/$ROUND"/ros2_*.err \
              "$BATCH_DIR/$ROUND"/ros2_*.py 2>/dev/null || true

        # ── Verify completion before declaring the batch done ─────
        # Residue dirs that already have a raw_scores.csv are skipped
        # (never recomputed, never required again). A batch is done when
        # every pending ligand's chunk finished — i.e. a fresh
        # <lig>.done_pre_exit exists — because the chunk script runs Rosetta
        # for EVERY anchor of the ligand; anchors where Rosetta found no
        # placements simply produce no CSV, which is a normal outcome.
        MISSING_LIGS=0
        for LIG in "${{PENDING_LIGS[@]}}"; do
            if [ ! -f "$BATCH_DIR/$ROUND/$LIG.done_pre_exit" ]; then
                echo "ERROR: ligand $LIG chunk did not finish (no fresh .done_pre_exit)"
                MISSING_LIGS=$((MISSING_LIGS + 1))
            fi
        done
        if [ "$MISSING_LIGS" -gt 0 ]; then
            echo "ERROR: $MISSING_LIGS/${{#PENDING_LIGS[@]}} pending ligand chunk(s) did not finish — NOT marking batch done"
            exit 1
        fi
        CSV_COUNT=$(find "$BATCH_DIR/$ROUND" -name 'raw_scores.csv' 2>/dev/null | wc -l)
        echo "All ${{#PENDING_LIGS[@]}} pending ligand chunks finished — $CSV_COUNT raw_scores.csv total (res dirs without placements are normal and skipped)"

        touch {output.done_flag}
        echo 'Rosetta round 2 complete'
        exit 0
        """


rule step6_sentinel:
    """Aggregate sentinel: depends on all rosetta_round2 completing."""
    localrule: True
    input:
        expand(os.path.join(TMP_ROOT, "batch_{batch_id}", "rosetta_round2.done"),
               batch_id=BATCH_IDS),
    output:
        STEP6_DONE,
    resources:
        mem_mb=500,
        cpus=1,
    shell:
        """
        touch {output}
        echo "Step 6 (rosetta_round2) complete"
        """
