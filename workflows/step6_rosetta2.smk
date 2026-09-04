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
        batches_dir     = BATCHES_DIR,
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

        BEST_ANCHOR_MAP=""   # derived below; trap may run before it is set

        # ── Trap: on ANY exit, reconcile per-ligand .done markers so a
        # ligand whose raw_scores.csv is present for its mapped target res
        # (or any res when no map) is always marked done — even on early
        # error exits. The batch done flag is NEVER touched here.
        _cleanup_and_exit() {{
            local exit_code=$?
            if [ -n "$BATCH_DIR" ] && [ -n "$ROUND" ] && [ "${{ANCHOR_COUNT:-0}}" -gt 0 ]; then
                for LIG_DIR in "$BATCH_DIR/$ROUND"/*/; do
                    [ -d "$LIG_DIR" ] || continue
                    LIG_NAME=$(basename "$LIG_DIR")
                    # When a best-anchor map is available, only the mapped
                    # (target) res counts; stale CSVs from older full-anchor
                    # runs must not satisfy completion.
                    LIG_TARGET=$(printf '%s\\n' "$BEST_ANCHOR_MAP" 2>/dev/null | grep -m1 "^$LIG_NAME," | cut -d, -f2)
                    if [ -n "$LIG_TARGET" ]; then
                        CSV_COUNT=$(find "$LIG_DIR/$LIG_TARGET" -name 'raw_scores.csv' 2>/dev/null | wc -l)
                    else
                        CSV_COUNT=$(find "$LIG_DIR" -name 'raw_scores.csv' 2>/dev/null | wc -l)
                    fi
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
        # Round 2 targets ONLY each ligand's round-1 best-score anchor, so a
        # ligand is complete once ONE raw_scores.csv exists (its target res).
        ANCHOR_COUNT=1

        # ── Derive the round-1 best-anchor map from step3 + step4 ──
        # The ligand set comes from step4 (top_ligands_round1.txt, already
        # materialized as round2 test_params by step5).  The best anchor per
        # ligand comes from step3's per-batch scores_round1.csv: its
        # 'ligand' column holds the best round-1 placement name, which
        # embeds the residue ("res<N>_<receptor>_<ligand>_<conf>_<idx>").
        # Map lines are "ligand,res"; missing/unparseable rows → the ligand
        # falls back to the full anchor list (old behavior) with a warning.
        BEST_ANCHOR_MAP=""
        if [ -s "$BATCH_DIR/scores_round1.csv" ]; then
            BEST_ANCHOR_MAP=$({params.python_bin} -c "
import csv, re, sys, os
batch_id = int(sys.argv[2])
batch_file = os.path.join(sys.argv[3], '%04d' % (batch_id // 1000), 'batch_%d.txt' % batch_id)
batch_ligands = []
if os.path.isfile(batch_file):
    with open(batch_file) as fh:
        for line in fh:
            parts = line.strip().split(',')
            if parts:
                batch_ligands.append(parts[0])
out = {{}}
with open(os.path.join(sys.argv[1], 'scores_round1.csv')) as fh:
    for row in csv.DictReader(fh):
        placement = (row.get('ligand') or '').strip()
        m = re.match(r'res(\d+)_', placement)
        if not m:
            continue
        base = placement
        for lig in sorted(batch_ligands, key=len, reverse=True):
            if lig in placement:
                base = lig
                break
        out[base] = m.group(1)
for base, res in out.items():
    print('%s,%s' % (base, res))
" "$BATCH_DIR" '{wildcards.batch_id}' '{params.batches_dir}' 2>/dev/null)
        fi
        if [ -z "$BEST_ANCHOR_MAP" ]; then
            echo "WARNING: no best-anchor map derivable from $BATCH_DIR/scores_round1.csv — round 2 will run ALL anchors for every ligand"
        else
            echo "Best-anchor map derived from scores_round1.csv: $(printf '%s\n' "$BEST_ANCHOR_MAP" | wc -l) ligand(s)"
        fi

        # ── Collect pending ligands (skip already-complete) ──────────
        PENDING_LIGS=()
        for TP_DIR in "$BATCH_DIR/$ROUND"/*/test_params/; do
            [ -d "$TP_DIR" ] || continue
            LIG_NAME=$(basename $(dirname "$TP_DIR"))
            # Only the round-1 best-score anchor is targeted in round 2, so a
            # ligand counts as complete once raw_scores.csv exists for THAT
            # anchor.  No mapping row → old behavior (count across all res
            # dirs, which are also skipped one by one by the chunk script).
            TARGET_RES=$(printf '%s\\n' "$BEST_ANCHOR_MAP" | grep -m1 "^$LIG_NAME," | cut -d, -f2)
            if [ -n "$TARGET_RES" ]; then
                EXISTING_CSV_COUNT=$(find "$BATCH_DIR/$ROUND/$LIG_NAME/$TARGET_RES" -name 'raw_scores.csv' 2>/dev/null | wc -l)
            else
                EXISTING_CSV_COUNT=$(find "$BATCH_DIR/$ROUND/$LIG_NAME" -name 'raw_scores.csv' 2>/dev/null | wc -l)
            fi
            if [ "$EXISTING_CSV_COUNT" -ge "$ANCHOR_COUNT" ] && [ "$ANCHOR_COUNT" -gt 0 ]; then
                if [ -n "$TARGET_RES" ]; then
                    echo "Ligand $LIG_NAME: already complete (target res $TARGET_RES has $EXISTING_CSV_COUNT CSVs), skipping submission"
                else
                    echo "Ligand $LIG_NAME: already complete ($EXISTING_CSV_COUNT CSVs), skipping submission"
                fi
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

            # Per-chunk residue map (ligand,target_residue): round 2 runs ONLY
            # each ligand's round-1 best-score anchor.  "ALL" = fallback (no
            # mapping row) → the chunk script runs the full anchor list.
            RES_MAP_FILE="$BATCH_DIR/$ROUND/ros2_chunk_${{chunk}}_residues.txt"
            : > "$RES_MAP_FILE"
            for LIG in "${{CHUNK_LIGS[@]}}"; do
                TARGET_RES=$(printf '%s\\n' "$BEST_ANCHOR_MAP" | grep -m1 "^$LIG," | cut -d, -f2)
                [ -z "$TARGET_RES" ] && TARGET_RES="ALL"
                echo "$LIG,$TARGET_RES" >> "$RES_MAP_FILE"
            done
            echo "Chunk $((chunk+1))/$NUM_CHUNKS residue map: $(wc -l < "$RES_MAP_FILE") ligands, $(grep -c ',ALL$' "$RES_MAP_FILE" || true) without mapping (will run all anchors)"

            CHUNK_SCRIPT="$BATCH_DIR/$ROUND/ros2_chunk_${{chunk}}.py"
            rm -f "$CHUNK_SCRIPT"
            cat > "$CHUNK_SCRIPT" << 'PYEOF'
import sys, os, glob, atexit, csv, shutil, tarfile

atexit.register(sys.stdout.flush)
atexit.register(sys.stderr.flush)

# ── Determine which ligand to process from LSF job array index ──────
job_index = int(os.environ.get('LSB_JOBINDEX', '1')) - 1  # LSF is 1-based
batch_dir = sys.argv[1]
round_name = sys.argv[2]
res_map_file = sys.argv[-1]
all_ligands = sys.argv[3:-1]

if job_index >= len(all_ligands):
    print(f'Job index {{job_index+1}} out of range (have {{len(all_ligands)}} ligands), exiting cleanly')
    sys.stdout.flush()
    os._exit(0)

lig_name = all_ligands[job_index]
tp_dir = os.path.join(batch_dir, round_name, lig_name, 'test_params')
print(f'[{{round_name}}] Processing ligand {{job_index+1}}/{{len(all_ligands)}}: {{lig_name}}')

sys.path.insert(0, '{params.discovery_root}/function/discovery')
from utils import cache_file_on_node, release_node_cached_files, run_rosetta_discovery_search
atexit.register(release_node_cached_files)

# Stage shared, immutable Rosetta inputs once per compute node so array
# elements on that node avoid repeatedly reading them from NFS.
local_rosetta_sif = cache_file_on_node('{params.sif_rosetta}')
local_motifs_file = cache_file_on_node('{params.motifs_file}')


def keep_best_scoring_pdb(res_dir):
    '''Extract only the highest-weighted placement PDB from Rosetta's archive.'''
    scores_path = os.path.join(res_dir, 'weighted_scores.csv')
    archive_path = os.path.join(res_dir, 'placements.tar.gz')
    if not os.path.isfile(scores_path) or not os.path.isfile(archive_path):
        return None

    best_row = None
    best_score = None
    try:
        with open(scores_path, newline='') as scores_fh:
            for row in csv.DictReader(scores_fh):
                try:
                    score = float(row.get('total', ''))
                except (TypeError, ValueError):
                    continue
                if best_score is None or score > best_score:
                    best_row = row
                    best_score = score

        if best_row is None:
            print('  No scored pose to retain in ' + scores_path)
            return None

        pdb_name = os.path.basename(best_row.get('file', ''))
        if not pdb_name.endswith('.pdb'):
            print('  WARNING: best score row has no valid PDB filename', file=sys.stderr)
            return None

        with tarfile.open(archive_path, 'r:gz') as archive:
            member = next(
                (item for item in archive.getmembers()
                 if item.isfile() and os.path.basename(item.name) == pdb_name),
                None,
            )
            if member is None:
                print('  WARNING: best pose ' + pdb_name + ' is missing from ' + archive_path,
                      file=sys.stderr)
                return None
            source = archive.extractfile(member)
            if source is None:
                print('  WARNING: could not read best pose ' + pdb_name + ' from archive',
                      file=sys.stderr)
                return None
            output_path = os.path.join(res_dir, pdb_name)
            with source, open(output_path, 'wb') as output_fh:
                shutil.copyfileobj(source, output_fh)

        print('  Kept best pose ' + pdb_name + ' (weighted score ' + str(best_score) + ')')
        return output_path
    except (OSError, csv.Error, tarfile.TarError) as exc:
        print('  WARNING: could not retain best-scoring PDB: ' + str(exc), file=sys.stderr)
        return None

# ── Round 2 runs ONLY the round-1 best-score anchor for this ligand ──
# (map derived by the controller from step3's scores_round1.csv: placement
# names embed the residue, "res<N>_...").  "ALL" fallback → run the full
# anchor list, matching the pre-best-anchor behavior.
target_res = 'ALL'
try:
    with open(res_map_file, 'r') as fh:
        for line in fh:
            parts = line.strip().split(',')
            if len(parts) >= 2 and parts[0] == lig_name:
                target_res = parts[1]
                break
except Exception as e:
    print(f'  WARNING: cannot read residue map {{res_map_file}}: {{e}} — running all anchors')

if target_res == 'ALL':
    residues = '{params.anchor_residues}'.split(',')
    print(f'  WARNING: no best-anchor mapping for {{lig_name}} — running ALL anchors')
else:
    residues = [target_res]
    print(f'  Round-1 best anchor for {{lig_name}}: {{target_res}} — round 2 runs only there')

for residue in residues:
    res = residue.strip()
    res_dir = os.path.join(batch_dir, round_name, lig_name, res)
    os.makedirs(res_dir, exist_ok=True)
    if os.path.isfile(os.path.join(res_dir, 'raw_scores.csv')):
        print(f'  SKIP: raw_scores.csv already exists for {{lig_name}}/{{res}}')
        continue
    success = run_rosetta_discovery_search(
        '{params.target_pdb}', res, local_motifs_file,
        tp_dir, '{params.discovery_root}',
        '{params.atr}', '{params.rep}', '{params.ddg}',
        extra_args_file='{params.extra_params}' or None,
        work_dir=res_dir,
        rosetta_sif=local_rosetta_sif
    )
    if not success:
        pdbs = glob.glob(os.path.join(res_dir, '*.pdb'))
        tars = glob.glob(os.path.join(res_dir, '*.tar.gz'))
        if not pdbs and not tars:
            print(f'WARNING: Rosetta produced no output for {{lig_name}}/{{res}}', file=sys.stderr)
        else:
            print(f'Rosetta: placements exist for {{lig_name}}/{{res}} despite non-zero exit')
    best_pdb = keep_best_scoring_pdb(res_dir)
    for f in glob.glob(os.path.join(res_dir, '*')):
        if f.endswith('.csv') or (best_pdb and os.path.abspath(f) == os.path.abspath(best_pdb)):
            continue
        try:
            if os.path.isdir(f):
                shutil.rmtree(f)
            else:
                os.remove(f)
        except Exception:
            pass

print('Rosetta done for ' + lig_name)
release_node_cached_files()
sys.stdout.flush()
sys.stderr.flush()
open(os.path.join(batch_dir, round_name, lig_name + '.done_pre_exit'), 'w').close()
os._exit(0)
PYEOF

            # ── Select the configured queue with the fewest pending jobs
            # owned by this controller's LSF user. Re-query before every
            # submission so new array elements affect the next decision.
            QLIST=({params.rosetta_queue})
            CHUNK_Q=""
            LOWEST_PENDING=""
            for Q in "${{QLIST[@]}}"; do
                Q_JOB_IDS=$(bjobs -p -q "$Q" -o "jobid" -noheader 2>/dev/null)
                Q_QUERY_STATUS=$?
                if [ "$Q_QUERY_STATUS" -eq 0 ]; then
                    Q_PENDING=$(printf '%s\n' "$Q_JOB_IDS" | sed '/^[[:space:]]*$/d' | wc -l)
                    echo "Queue $Q: $Q_PENDING of my jobs pending"
                    if [ -z "$LOWEST_PENDING" ] || [ "$Q_PENDING" -lt "$LOWEST_PENDING" ]; then
                        CHUNK_Q="$Q"
                        LOWEST_PENDING="$Q_PENDING"
                    fi
                else
                    echo "WARNING: could not read my pending jobs for queue $Q"
                fi
            done
            if [ -z "$CHUNK_Q" ]; then
                CHUNK_Q="${{QLIST[0]}}"
                echo "WARNING: no per-user queue data available; falling back to $CHUNK_Q"
            else
                echo "Selected queue $CHUNK_Q for chunk $((chunk+1))/$NUM_CHUNKS ($LOWEST_PENDING pending)"
            fi

            ARRAY_JOB_NAME="smk_ros2_c${{chunk}}"
            ARRAY_JOB_ID=$(bsub \
                -W {params.rosetta_walltime} \
                -q "$CHUNK_Q" \
                -M 1500 \
                -n 1 \
                -R 'span[hosts=1] rusage[mem=1500]' \
                -J "${{ARRAY_JOB_NAME}}[1-${{N}}]" \
                -o "$BATCH_DIR/$ROUND/ros2_chunk_${{chunk}}_%I.out" \
                -e "$BATCH_DIR/$ROUND/ros2_chunk_${{chunk}}_%I.err" \
                {params.python_bin} "$CHUNK_SCRIPT" "$BATCH_DIR" "$ROUND" $LIG_LIST "$RES_MAP_FILE" \
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
        CLEAN_CHECK=0

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
                # `bjobs` without `-a` reports only active array elements. Once
                # every element is DONE/EXIT (or has aged out of mbatchd), it
                # returns no status rows. Therefore empty output means this
                # array has no live elements; it must not keep the controller
                # waiting forever. Two clean polls below guard against a
                # transient empty response, and the final .done_pre_exit check
                # prevents a query glitch from marking incomplete work done.
                ELEM_STATS=$(bjobs -o "stat" -noheader "$ARRAY_JID" 2>/dev/null | tr -d ' ')
                if echo "$ELEM_STATS" | grep -qE 'RUN|PEND|WAIT|USUSP|PSUSP|SSUSP'; then
                    STILL_RUNNING=1
                fi
            done

            if [ "$STILL_RUNNING" -eq 0 ]; then
                # Require two consecutive "all finished" readings (30 s apart) so a
                # transient bjobs glitch cannot end the wait while arrays are PEND.
                CLEAN_CHECK=$((CLEAN_CHECK + 1))
                if [ "$CLEAN_CHECK" -ge 2 ]; then
                    ALL_FINISHED=1
                    CSV_COUNT=$(find "$BATCH_DIR/$ROUND" -name 'raw_scores.csv' 2>/dev/null | wc -l)
                    echo "All array jobs finished — $CSV_COUNT raw_scores.csv files produced"
                    break
                fi
            else
                CLEAN_CHECK=0
            fi

            if [ $(( (NOW - WAIT_START) % 600 )) -lt 30 ]; then
                TOT_RUN=0; TOT_PEND=0
                for ARRAY_JID in "${{ARRAY_JOB_IDS[@]}}"; do
                    ES=$(bjobs -o "stat" -noheader "$ARRAY_JID" 2>/dev/null | tr -d ' ')
                    # NOTE: `|| true` (NOT `|| echo 0`): grep -c already prints
                    # "0" when nothing matches; `|| echo 0` prints a second "0",
                    # and the resulting "0\\n0" makes $(( )) throw an arithmetic
                    # error, which makes bash abort the whole wait loop early.
                    TOT_RUN=$((TOT_RUN + $(echo "$ES" | grep -c "RUN" || true)))
                    TOT_PEND=$((TOT_PEND + $(echo "$ES" | grep -c "PEND" || true)))
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
            # Count only the mapped (target) res when a best-anchor map is
            # available; stale CSVs from older full-anchor runs must not
            # satisfy completion.
            LIG_TARGET=$(printf '%s\\n' "$BEST_ANCHOR_MAP" | grep -m1 "^$LIG," | cut -d, -f2)
            if [ -n "$LIG_TARGET" ]; then
                CSV_COUNT=$(find "$LIG_DIR/$LIG_TARGET" -name 'raw_scores.csv' 2>/dev/null | wc -l)
            else
                CSV_COUNT=$(find "$LIG_DIR" -name 'raw_scores.csv' 2>/dev/null | wc -l)
            fi
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
              "$BATCH_DIR/$ROUND"/ros2_chunk_*_residues.txt \
              "$BATCH_DIR/$ROUND"/ros2_*.out \
              "$BATCH_DIR/$ROUND"/ros2_*.err \
              "$BATCH_DIR/$ROUND"/ros2_*.py 2>/dev/null || true

        # ── Verify completion before declaring the batch done ─────
        # Residue dirs that already have a raw_scores.csv are skipped
        # (never recomputed, never required again). A batch is done when
        # every pending ligand's chunk finished — i.e. a fresh
        # <lig>.done_pre_exit exists — because the chunk script runs Rosetta
        # at the ligand's round-1 best-score anchor only; anchors where
        # Rosetta found no placements simply produce no CSV, which is a
        # normal outcome.
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
