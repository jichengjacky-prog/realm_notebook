# ═══════════════════════════════════════════════════════════════════════════════
# MASTER CONTROLLER — Orchestrates the 8 pipeline steps as separate Snakemake runs.
#
# Each step is an independent Snakemake workflow in workflows/stepN_*.smk.
# This controller runs them sequentially via shell rules, so each step gets
# its own small DAG and can use step-appropriate --resources limits.
#
# Usage:
#   snakemake -s Snakefile --configfile yaml/config_state3.yaml --profile profile/lsf
#
#   # Or use the bash controller (simpler, recommended):
#   ./run_pipeline.sh yaml/config_state3.yaml
# ═══════════════════════════════════════════════════════════════════════════════

import os

configfile: "config.yaml"
include: "workflows/shared_config.smk"

# ── Profile ─────────────────────────────────────────────────────────────
PROFILE       = config.get("profile", "profile/lsf")
WORKFLOW_DIR  = os.path.join(REALM_LOCATION, "workflows")

# ── Step sentinel files ──────────────────────────────────────────────────
STEP1_DONE = os.path.join(OUTPUT_DIR, ".step1_extract.done")
STEP2_DONE = os.path.join(OUTPUT_DIR, ".step2_rosetta1.done")
STEP3_DONE = os.path.join(OUTPUT_DIR, ".step3_score1.done")
STEP4_DONE = os.path.join(OUTPUT_DIR, "top_ligands_round1.txt")
STEP5_DONE = os.path.join(OUTPUT_DIR, ".step5_conformers.done")
STEP6_DONE = os.path.join(OUTPUT_DIR, ".step6_rosetta2.done")
STEP7_DONE = os.path.join(OUTPUT_DIR, ".step7_score2.done")
STEP8_DONE = os.path.join(OUTPUT_DIR, "cleanup.done")

# ── Derive configfile path ───────────────────────────────────────────────
CONFIGFILE = config.get("_master_configfile", "yaml/config_state3.yaml")


# ═══════════════════════════════════════════════════════════════════════════════
# Master target
# ═══════════════════════════════════════════════════════════════════════════════

rule all:
    input: STEP8_DONE


# ═══════════════════════════════════════════════════════════════════════════════
# Step 1: Parse shapedb & extract library params  (load=1, many parallel jobs)
# ═══════════════════════════════════════════════════════════════════════════════

rule step1_extract:
    localrule: True
    output: STEP1_DONE
    resources:
        mem_mb=1000,
        cpus=1,
    log:
        os.path.join(OUTPUT_DIR, "step1_controller.log"),
    shell:
        """
        set -e
        echo "=== STEP 1: Parse shapedb & extract params ==="
        {SNAKEMAKE_BIN} \
            -s {WORKFLOW_DIR}/step1_parse_and_extract.smk \
            --configfile {CONFIGFILE} \
            --profile {PROFILE} \
            --resources load=1000 \
            --keep-going \
            --rerun-incomplete \
            2>&1 | tee {log}
        echo "Step 1 complete"
        """


# ═══════════════════════════════════════════════════════════════════════════════
# Step 2: Rosetta round 1  (load=10 per job, 1000 total → ~100 concurrent)
# ═══════════════════════════════════════════════════════════════════════════════

rule step2_rosetta1:
    localrule: True
    input: STEP1_DONE
    output: STEP2_DONE
    resources:
        mem_mb=1000,
        cpus=1,
    log:
        os.path.join(OUTPUT_DIR, "step2_controller.log"),
    shell:
        """
        set -e
        echo "=== STEP 2: Rosetta discovery round 1 ==="
        {SNAKEMAKE_BIN} \
            -s {WORKFLOW_DIR}/step2_rosetta1.smk \
            --configfile {CONFIGFILE} \
            --profile {PROFILE} \
            --resources load=1000 \
            --keep-going \
            --rerun-incomplete \
            2>&1 | tee {log}
        echo "Step 2 complete"
        """


# ═══════════════════════════════════════════════════════════════════════════════
# Step 3: Score round 1
# ═══════════════════════════════════════════════════════════════════════════════

rule step3_score1:
    localrule: True
    input: STEP2_DONE
    output: STEP3_DONE
    resources:
        mem_mb=1000,
        cpus=1,
    log:
        os.path.join(OUTPUT_DIR, "step3_controller.log"),
    shell:
        """
        set -e
        echo "=== STEP 3: Score round 1 placements ==="
        {SNAKEMAKE_BIN} \
            -s {WORKFLOW_DIR}/step3_score1.smk \
            --configfile {CONFIGFILE} \
            --profile {PROFILE} \
            --resources load=1000 \
            --keep-going \
            --rerun-incomplete \
            2>&1 | tee {log}
        echo "Step 3 complete"
        """


# ═══════════════════════════════════════════════════════════════════════════════
# Step 4: Filter top-N
# ═══════════════════════════════════════════════════════════════════════════════

rule step4_filter:
    localrule: True
    input: STEP3_DONE
    output: STEP4_DONE
    resources:
        mem_mb=8000,
        cpus=1,
    log:
        os.path.join(OUTPUT_DIR, "step4_controller.log"),
    shell:
        """
        set -e
        echo "=== STEP 4: Filter top-N ligands ==="
        {SNAKEMAKE_BIN} \
            -s {WORKFLOW_DIR}/step4_filter.smk \
            --configfile {CONFIGFILE} \
            --resources mem_mb=8000 \
            2>&1 | tee {log}
        echo "Step 4 complete"
        """


# ═══════════════════════════════════════════════════════════════════════════════
# Step 5: Generate CDPKit conformers
# ═══════════════════════════════════════════════════════════════════════════════

rule step5_conformers:
    localrule: True
    input: STEP4_DONE
    output: STEP5_DONE
    resources:
        mem_mb=1000,
        cpus=1,
    log:
        os.path.join(OUTPUT_DIR, "step5_controller.log"),
    shell:
        """
        set -e
        echo "=== STEP 5: Generate CDPKit conformers ==="
        {SNAKEMAKE_BIN} \
            -s {WORKFLOW_DIR}/step5_conformers.smk \
            --configfile {CONFIGFILE} \
            --profile {PROFILE} \
            --resources load=1000 \
            --keep-going \
            --rerun-incomplete \
            2>&1 | tee {log}
        echo "Step 5 complete"
        """


# ═══════════════════════════════════════════════════════════════════════════════
# Step 6: Rosetta round 2
# ═══════════════════════════════════════════════════════════════════════════════

rule step6_rosetta2:
    localrule: True
    input: STEP5_DONE
    output: STEP6_DONE
    resources:
        mem_mb=1000,
        cpus=1,
    log:
        os.path.join(OUTPUT_DIR, "step6_controller.log"),
    shell:
        """
        set -e
        echo "=== STEP 6: Rosetta discovery round 2 ==="
        {SNAKEMAKE_BIN} \
            -s {WORKFLOW_DIR}/step6_rosetta2.smk \
            --configfile {CONFIGFILE} \
            --profile {PROFILE} \
            --resources load=10000 \
            --keep-going \
            --rerun-incomplete \
            2>&1 | tee {log}
        echo "Step 6 complete"
        """


# ═══════════════════════════════════════════════════════════════════════════════
# Step 7: Score round 2
# ═══════════════════════════════════════════════════════════════════════════════

rule step7_score2:
    localrule: True
    input: STEP6_DONE
    output: STEP7_DONE
    resources:
        mem_mb=1000,
        cpus=1,
    log:
        os.path.join(OUTPUT_DIR, "step7_controller.log"),
    shell:
        """
        set -e
        echo "=== STEP 7: Score round 2 placements ==="
        {SNAKEMAKE_BIN} \
            -s {WORKFLOW_DIR}/step7_score2.smk \
            --configfile {CONFIGFILE} \
            --profile {PROFILE} \
            --resources load=1000 \
            --keep-going \
            --rerun-incomplete \
            2>&1 | tee {log}
        echo "Step 7 complete"
        """


# ═══════════════════════════════════════════════════════════════════════════════
# Step 8: Aggregate & cleanup
# ═══════════════════════════════════════════════════════════════════════════════

rule step8_aggregate:
    localrule: True
    input: STEP7_DONE
    output: STEP8_DONE
    resources:
        mem_mb=20000,
        cpus=1,
    log:
        os.path.join(OUTPUT_DIR, "step8_controller.log"),
    shell:
        """
        set -e
        echo "=== STEP 8: Aggregate results & cleanup ==="
        {SNAKEMAKE_BIN} \
            -s {WORKFLOW_DIR}/step8_aggregate.smk \
            --configfile {CONFIGFILE} \
            --resources mem_mb=20000 \
            2>&1 | tee {log}
        echo "Step 8 complete — pipeline finished!"
        """
