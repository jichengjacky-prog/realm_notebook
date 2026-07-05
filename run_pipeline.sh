#!/bin/bash
# ═══════════════════════════════════════════════════════════════════════════════
# Rosetta Discovery Pipeline — Modular Step Controller with Batch Slicing
#
# Each step is an independent Snakemake workflow. For batch-heavy steps
# (1,2,3,5,6,7), batches are processed in slices of BATCH_SLICE_SIZE to keep
# DAG building fast (~2s per slice instead of ~5min for all 30K).
#
# Usage:
#   ./run_pipeline.sh [configfile] [start_step] [end_step]
#   ./run_pipeline.sh yaml/config_state3.yaml          # all steps
#   ./run_pipeline.sh yaml/config_state3.yaml 2        # resume from step 2
#   ./run_pipeline.sh yaml/config_state3.yaml 2 5      # steps 2-5 only
# ═══════════════════════════════════════════════════════════════════════════════

set -euo pipefail

CONFIGFILE="${1:-yaml/config_state3.yaml}"
START_STEP="${2:-1}"
END_STEP="${3:-8}"

PROFILE="profile/lsf"
WORKFLOW_DIR="$(dirname "$0")/workflows"
REALM_DIR="$(dirname "$0")"
BATCH_SLICE_SIZE=500      # batches per Snakemake invocation
MAX_CONCURRENT_SLICES=12  # max slices to run simultaneously

# ── Derive paths from config ─────────────────────────────────────────────
# Use the realm_env Python (which has PyYAML) to parse config.
# The local yaml/ directory shadows the PyYAML package for system python3,
# causing "module 'yaml' has no attribute 'safe_load'" errors.
_REALM_PYTHON="/home/ji.cheng4-umw/miniforge3/envs/realm_env/bin/python3.11"

PYTHON_BIN=$($_REALM_PYTHON -c "
import yaml
with open('$CONFIGFILE') as f:
    cfg = yaml.safe_load(f)
print(cfg.get('python_bin', 'python3.11'))
" 2>/dev/null || echo "python3.11")

SNAKEMAKE=$($_REALM_PYTHON -c "
import yaml
with open('$CONFIGFILE') as f:
    cfg = yaml.safe_load(f)
print(cfg.get('snakemake_bin', 'snakemake'))
" 2>/dev/null || echo "snakemake")

OUTPUT_DIR=$($_REALM_PYTHON -c "
import yaml
with open('$CONFIGFILE') as f:
    cfg = yaml.safe_load(f)
print(cfg.get('output_dir', 'output'))
" 2>/dev/null || echo "output/UNKNOWN_CONFIG_DIR")

# ── Ensure output directory exists ──────────────────────────────────────
mkdir -p "$OUTPUT_DIR"

# ── Discover total batch count ───────────────────────────────────────────
# Use subshells with || true to prevent set -eo pipefail from killing the
# script when the directory doesn't exist yet.
TOTAL_BATCHES=$( (find "$OUTPUT_DIR/batches" -name "batch_*.txt" 2>/dev/null || true) | wc -l)
if [ "$TOTAL_BATCHES" -eq 0 ]; then
    TOTAL_BATCHES=$( (find "$OUTPUT_DIR/tmp" -maxdepth 2 -name "extract_params.done" 2>/dev/null || true) | wc -l)
fi
echo "Total batches discovered: $TOTAL_BATCHES"

# ── Step sentinel files ──────────────────────────────────────────────────
STEP1_DONE="$OUTPUT_DIR/.step1_extract.done"
STEP2_DONE="$OUTPUT_DIR/.step2_rosetta1.done"
STEP3_DONE="$OUTPUT_DIR/.step3_score1.done"
STEP4_DONE="$OUTPUT_DIR/top_ligands_round1.txt"
STEP5_DONE="$OUTPUT_DIR/.step5_conformers.done"
STEP6_DONE="$OUTPUT_DIR/.step6_rosetta2.done"
STEP7_DONE="$OUTPUT_DIR/.step7_score2.done"
STEP8_DONE="$OUTPUT_DIR/cleanup.done"

echo "╔══════════════════════════════════════════════════════════════╗"
echo "║   Rosetta Discovery Pipeline — Modular + Sliced            ║"
echo "╠══════════════════════════════════════════════════════════════╣"
echo "║ Config:      $CONFIGFILE"
echo "║ Output dir:  $OUTPUT_DIR"
echo "║ Steps:       $START_STEP → $END_STEP"
echo "║ Slice size:  $BATCH_SLICE_SIZE batches"
echo "║ Total:       $TOTAL_BATCHES batches"
echo "╚══════════════════════════════════════════════════════════════╝"
echo ""

# ═══════════════════════════════════════════════════════════════════════════════
# Helper: run a sliced step (iterates over batch slices)
# ═══════════════════════════════════════════════════════════════════════════════
run_sliced_step() {
    local step_num="$1"
    local step_name="$2"
    local workflow_file="$3"
    local sentinel="$4"
    local extra_args="${5:-}"

    if [ "$step_num" -lt "$START_STEP" ] || [ "$step_num" -gt "$END_STEP" ]; then
        echo "--- Skipping step $step_num ($step_name) — outside range ---"
        return 0
    fi

    if [ -f "$sentinel" ]; then
        echo "--- Step $step_num ($step_name) already complete ---"
        return 0
    fi

    echo ""
    echo "╔══════════════════════════════════════════════════════════════╗"
    echo "║  STEP $step_num: $step_name (sliced, $BATCH_SLICE_SIZE batches/slice)"
    echo "╚══════════════════════════════════════════════════════════════╝"
    echo ""

    local total_slices=$(( (TOTAL_BATCHES + BATCH_SLICE_SIZE - 1) / BATCH_SLICE_SIZE ))
    local pids=()
    local slice_sentinels=()
    local running=0

    for ((slice=0; slice<total_slices; slice++)); do
        local start=$((slice * BATCH_SLICE_SIZE))
        local end=$((start + BATCH_SLICE_SIZE))
        if [ "$end" -gt "$TOTAL_BATCHES" ]; then
            end=$TOTAL_BATCHES
        fi

        local slice_sentinel="$OUTPUT_DIR/.step${step_num}_${step_name}_slice_${start}_${end}.done"
        slice_sentinels+=("$slice_sentinel")

        # Skip slice if its sentinel already exists
        if [ -f "$slice_sentinel" ]; then
            echo "--- Slice $((slice+1))/$total_slices (batches $start–$((end-1))) [$step_name] — already done ---"
            continue
        fi

        # ── Concurrency limiter: wait if we already have MAX_CONCURRENT_SLICES running ──
        while [ "$running" -ge "$MAX_CONCURRENT_SLICES" ]; do
            # Wait for any one child to finish, then reap it
            wait -n "${pids[@]}" 2>/dev/null || true
            running=$((running - 1))
            # Remove completed PIDs from the tracking array so later
            # "wait for all" doesn't trip over already-reaped children
            local _new_pids=()
            for _p in "${pids[@]}"; do
                if kill -0 "$_p" 2>/dev/null; then
                    _new_pids+=("$_p")
                fi
            done
            pids=("${_new_pids[@]}")
        done

        echo "--- Slice $((slice+1))/$total_slices (batches $start–$((end-1))) [$step_name] launched ---"

        (
            $SNAKEMAKE \
                -s "$WORKFLOW_DIR/$workflow_file" \
                --configfile "$CONFIGFILE" \
                --profile "$PROFILE" \
                --keep-going \
                --rerun-incomplete \
                --config batch_start="$start" batch_end="$end" \
                $extra_args \
                2>&1 | tee "$OUTPUT_DIR/step${step_num}_${step_name}_slice${slice}.log"
            exit_code=${PIPESTATUS[0]}
            if [ "$exit_code" -eq 0 ]; then
                touch "$slice_sentinel"
                echo "--- Slice $((slice+1))/$total_slices DONE [$step_name] ---"
            else
                echo "ERROR: Step $step_num slice $slice failed (exit $exit_code)" >&2
                exit "$exit_code"
            fi
        ) &
        pids+=($!)
        running=$((running + 1))
    done

    # Wait for all remaining slices to finish
    local failed=0
    for pid in "${pids[@]}"; do
        if ! wait "$pid"; then
            failed=1
        fi
    done

    if [ "$failed" -ne 0 ]; then
        echo "ERROR: One or more slices of step $step_num failed"
        exit 1
    fi

    # All slices complete — create aggregate sentinel
    touch "$sentinel"
    echo "--- Step $step_num ($step_name) complete ($total_slices slices) ---"
}

# ═══════════════════════════════════════════════════════════════════════════════
# Helper: run a single (unsliced) step
# ═══════════════════════════════════════════════════════════════════════════════
run_step() {
    local step_num="$1"
    local step_name="$2"
    local workflow_file="$3"
    local sentinel="$4"
    local extra_args="${5:-}"

    if [ "$step_num" -lt "$START_STEP" ] || [ "$step_num" -gt "$END_STEP" ]; then
        echo "--- Skipping step $step_num ($step_name) — outside range ---"
        return 0
    fi

    if [ -f "$sentinel" ]; then
        echo "--- Step $step_num ($step_name) already complete ---"
        return 0
    fi

    echo ""
    echo "╔══════════════════════════════════════════════════════════════╗"
    echo "║  STEP $step_num: $step_name"
    echo "╚══════════════════════════════════════════════════════════════╝"
    echo ""

    $SNAKEMAKE \
        -s "$WORKFLOW_DIR/$workflow_file" \
        --configfile "$CONFIGFILE" \
        --profile "$PROFILE" \
        --keep-going \
        --rerun-incomplete \
        $extra_args \
        2>&1 | tee "$OUTPUT_DIR/step${step_num}_${step_name}.log"

    local exit_code=${PIPESTATUS[0]}
    if [ "$exit_code" -ne 0 ]; then
        echo "ERROR: Step $step_num ($step_name) failed with exit code $exit_code"
        exit "$exit_code"
    fi
    echo "--- Step $step_num ($step_name) completed successfully ---"
}

# ═══════════════════════════════════════════════════════════════════════════════
# Run each step
# ═══════════════════════════════════════════════════════════════════════════════

# Step 1 is NOT sliced — it reads the shapedb list and creates batch files
# internally.  Slicing depends on batch files that don't exist yet.
run_step 1 "parse_and_extract" \
    "step1_parse_and_extract.smk" \
    "$STEP1_DONE" \
    "--resources load=1000"

# Re-discover batch count now that step 1 has created batch files.
TOTAL_BATCHES=$( (find "$OUTPUT_DIR/batches" -name "batch_*.txt" 2>/dev/null || true) | wc -l)
if [ "$TOTAL_BATCHES" -eq 0 ]; then
    TOTAL_BATCHES=$( (find "$OUTPUT_DIR/tmp" -maxdepth 2 -name "extract_params.done" 2>/dev/null || true) | wc -l)
fi
echo "Total batches after step 1: $TOTAL_BATCHES"

if [ "$TOTAL_BATCHES" -eq 0 ]; then
    echo "ERROR: No batch files found in $OUTPUT_DIR/batches/ or $OUTPUT_DIR/tmp/" >&2
    echo "Step 1 (parse_and_extract) may have failed silently — check its log." >&2
    exit 1
fi

run_sliced_step 2 "rosetta_round1" \
    "step2_rosetta1.smk" \
    "$STEP2_DONE" \
    "--resources load=20"

run_sliced_step 3 "score_round1" \
    "step3_score1.smk" \
    "$STEP3_DONE" \
    "--resources load=100"

run_step 4 "filter_top" \
    "step4_filter.smk" \
    "$STEP4_DONE" \
    "--resources mem_mb=8000"

run_sliced_step 5 "generate_conformers" \
    "step5_conformers.smk" \
    "$STEP5_DONE" \
    "--resources load=1000"

run_sliced_step 6 "rosetta_round2" \
    "step6_rosetta2.smk" \
    "$STEP6_DONE" \
    "--resources load=20"

run_sliced_step 7 "score_round2" \
    "step7_score2.smk" \
    "$STEP7_DONE" \
    "--resources load=100"

run_step 8 "aggregate_cleanup" \
    "step8_aggregate.smk" \
    "$STEP8_DONE" \
    "--resources mem_mb=20000"

echo ""
echo "╔══════════════════════════════════════════════════════════════╗"
echo "║   PIPELINE COMPLETE                                        ║"
echo "║   Final results: $OUTPUT_DIR/top_ligands.csv               ║"
echo "╚══════════════════════════════════════════════════════════════╝"
