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
    resources:
        load=10,
        mem_mb=2000,
        cpus=1,
        queue=LSF_QUEUE_ROSETTA,
        walltime=LSF_WALLTIME_DEFAULT,
    log:
        os.path.join(TMP_ROOT, "batch_{batch_id}", "rosetta_round2.log"),
    shell:
        """
        # NOTE: set +e (not set -e) because some per-ligand Rosetta jobs
        # may exit non-zero (no placements found) which is normal.
        BATCH_DIR=$(dirname {output.done_flag})
        mkdir -p $BATCH_DIR
        echo "Starting Rosetta round 2 for batch in $BATCH_DIR"

        JOB_IDS=()
        LIG_NAMES=()
        ANCHOR_COUNT=$(echo '{params.anchor_residues}' | tr ',' '\n' | sed '/^[[:space:]]*$/d' | wc -l)
        for TP_DIR in "$BATCH_DIR"/round2/*/test_params/; do
            [ -d "$TP_DIR" ] || continue
            LIG_NAME=$(basename $(dirname "$TP_DIR"))
            # Skip if ligand already has all expected CSVs (idempotent across restarts)
            EXISTING_CSV_COUNT=$(find "$BATCH_DIR/round2/$LIG_NAME" -name 'weighted_scores.csv' 2>/dev/null | wc -l)
            if [ "$EXISTING_CSV_COUNT" -ge "$ANCHOR_COUNT" ] && [ "$ANCHOR_COUNT" -gt 0 ]; then
                echo "Ligand $LIG_NAME: already complete ($EXISTING_CSV_COUNT CSVs), skipping submission"
                rm -f "$BATCH_DIR/round2/ros2_$LIG_NAME.out" \
                      "$BATCH_DIR/round2/ros2_$LIG_NAME.err" \
                      "$BATCH_DIR/round2/ros2_$LIG_NAME.py"
                touch "$BATCH_DIR/round2/$LIG_NAME.done"
                continue
            fi
            LIG_NAMES+=("$LIG_NAME")
            JOB_NAME="smk_ros2_${{LIG_NAME:0:20}}"
            SCRIPT="$BATCH_DIR/round2/ros2_$LIG_NAME.py"
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
    res_dir = os.path.join(batch_dir, 'round2', lig_name, res)
    os.makedirs(res_dir, exist_ok=True)
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
            print(f'WARNING: Round 2 Rosetta produced no output for {{lig_name}}/{{res}}', file=sys.stderr)
        else:
            print(f'R2 Rosetta: placements exist for {{lig_name}}/{{res}} despite non-zero exit')
    for f in glob.glob(os.path.join(res_dir, '*')):
        if not f.endswith('.csv'):
            try:
                if os.path.isdir(f):
                    import shutil; shutil.rmtree(f)
                else:
                    os.remove(f)
            except Exception:
                pass
print('Rosetta round 2 done for ' + lig_name)
PYEOF
            JOB_ID=$(bsub -q {params.rosetta_queue} \
                -W {params.rosetta_walltime} \
                -M 8000 \
                -n 1 \
                -R 'span[hosts=1] rusage[mem=8000]' \
                -J "$JOB_NAME" \
                -o "$BATCH_DIR/round2/ros2_$LIG_NAME.out" \
                -e "$BATCH_DIR/round2/ros2_$LIG_NAME.err" \
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

        echo "Waiting for all per-ligand jobs to finish..."

        while true; do
            STILL_RUNNING=0
            for JID in "${{JOB_IDS[@]}}"; do
                STAT=$(bjobs -o stat -noheader "$JID" 2>/dev/null | tr -d ' ')
                if [ "$STAT" = "RUN" ] || [ "$STAT" = "PEND" ]; then
                    STILL_RUNNING=1
                    break
                fi
            done
            if [ "$STILL_RUNNING" -eq 0 ]; then
                CSV_COUNT=$(find "$BATCH_DIR"/round2 -name 'weighted_scores.csv' 2>/dev/null | wc -l)
                echo "All per-ligand jobs finished — $CSV_COUNT weighted_scores.csv files produced"
                break
            fi
            sleep 30
        done

        # ── Per-ligand cleanup: remove logs if successful, touch done ──
        for i in "${{!LIG_NAMES[@]}}"; do
            LIG="${{LIG_NAMES[$i]}}"
            LIG_DIR="$BATCH_DIR/round2/$LIG"
            CSV_FILE=$(find "$LIG_DIR" -name 'weighted_scores.csv' 2>/dev/null | head -1)
            if [ -n "$CSV_FILE" ]; then
                rm -f "$BATCH_DIR/round2/ros2_$LIG.out" \
                      "$BATCH_DIR/round2/ros2_$LIG.err" \
                      "$BATCH_DIR/round2/ros2_$LIG.py"
                touch "$BATCH_DIR/round2/$LIG.done"
                echo "Ligand $LIG: SUCCESS — cleaned logs, touched done"
            else
                echo "Ligand $LIG: FAILED or no output — keeping logs for debugging"
            fi
        done

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
