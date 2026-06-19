#!/bin/bash
# Submit the Rosetta discovery Snakemake pipeline to LSF.
#
# Usage:
#   ./run_pipeline.sh [configfile] [batch_slices]            # default: config_state3.yaml  10
#   ./run_pipeline.sh config_state4.yaml 20
#
#   batch_slices: number of LSF jobs to split round-2 across.
#       Lower = fewer master jobs = faster DAG building per job,
#       but each job processes more sub-batches.
#       Recommended: 10-20 for <1000 batch files; 5-10 for >1000.
#
# The master Snakemake process runs inside a 12h LSF job.  It submits
# individual rule jobs (Rosetta, scoring, etc.) as separate LSF jobs via
# the cluster-generic executor defined in profile/lsf/.
#
# Profile settings (profile/lsf/config.yaml):
#   cores: 2000              max concurrent cores
#   jobs: 2000               max parallel LSF jobs

set -euo pipefail

CONFIGFILE="${1:-config_state3.yaml}"
BATCH_SLICES="${2:-10}"
SNAKEMAKE="/home/ji.cheng4-umw/miniforge3/envs/realm_env/bin/snakemake"
PROFILE="profile/lsf"

echo "=== Rosetta Discovery Pipeline ==="
echo "Config:       $CONFIGFILE"
echo "Batch slices: $BATCH_SLICES"
echo "Profile:      $PROFILE"
echo ""



# Submit the master job (2 cores, 16 GB, 30 days walltime)
bsub -n 2 -q long -R "rusage[mem=16000]" -W 720:00 -J "smk_master" \
    "$SNAKEMAKE --keep-going --batch joint_call_rule=1/${BATCH_SLICES} --profile $PROFILE --configfile $CONFIGFILE"

echo "Master job submitted."
echo "Monitor with: bjobs -u $(whoami)"
echo "            cjobs -r"
