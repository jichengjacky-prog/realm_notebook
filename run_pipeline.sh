#!/bin/bash
# Submit the Rosetta discovery Snakemake pipeline to LSF.
#
# Usage:
#   ./run_pipeline.sh [configfile]            # default: config_state3.yaml
#   ./run_pipeline.sh config_state4.yaml
#
# The master Snakemake process runs inside a 12h LSF job.  It submits
# individual rule jobs (Rosetta, scoring, etc.) as separate LSF jobs via
# the cluster-generic executor defined in profile/lsf/.
#
# Profile settings (profile/lsf/config.yaml):
#   cores: 500               max concurrent cores
#   jobs: 500                max parallel LSF jobs
#   max-jobs-per-second: 100 rate limit on job submission
# ./run_pipeline.sh is a simple wrapper that submits the master Snakemake job to LSF.
# usage:
# ./run_pipeline.sh yaml/config_state3.yaml

set -euo pipefail

CONFIGFILE="${1:-config_state3.yaml}"
SNAKEMAKE="/home/ji.cheng4-umw/miniforge3/envs/realm_env/bin/snakemake"
PROFILE="profile/lsf"

echo "=== Rosetta Discovery Pipeline ==="
echo "Config:  $CONFIGFILE"
echo "Profile: $PROFILE"
echo ""

# Submit the master job (1 core, 16 GB, 7 days walltime)
bsub -n 1 -q long -R "rusage[mem=16000]" -W 300:00 -J "smk_master" \
    "$SNAKEMAKE --profile $PROFILE --configfile $CONFIGFILE"

echo "Master job submitted."
echo "Monitor with: bjobs -u $(whoami)"
echo "            cjobs -r"
