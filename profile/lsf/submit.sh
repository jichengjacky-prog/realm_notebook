#!/bin/bash
# Snakemake cluster-generic submit wrapper for LSF
# Snakemake appends the jobscript path as the final argument.
# bsub requires the script via stdin redirection, not as a positional arg.
# Usage: submit.sh [bsub options] <jobscript>

# Extract the jobscript path (last argument) from bsub options (all others)
JOBSCRIPT="${@: -1}"
BSB_ARGS="${@:1:$#-1}"

# Submit via stdin redirection and extract the numeric job ID
bsub -u '' $BSB_ARGS < "$JOBSCRIPT" 2>&1 | sed -n 's/.*Job <\([0-9]*\)>.*/\1/p'
