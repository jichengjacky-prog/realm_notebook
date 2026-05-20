#!/bin/bash
set -euo pipefail

###############################################################################
# LSF SCRATCH WRAPPER (NO MOUNT VERSION)
#
# This script:
#   • Determines a unique scratch directory for each job (array‑safe)
#   • Creates a plain directory under /tmp (no loopback, no mount)
#   • Exports SCRATCH_DIR for the worker script
#   • Runs the worker script inside the scratch directory
#   • Cleans up automatically on exit
###############################################################################

#test print of the job index
echo "LSB_JOBINDEX=$LSB_JOBINDEX"

###############################################################################
# CONFIGURATION
###############################################################################

SCRATCH_SIZE_GB="${SCRATCH_SIZE_GB:-10}"


###############################################################################
# DETERMINE SCRATCH DIRECTORY
###############################################################################

if [ -z "${SCRATCH_DIR:-}" ]; then
    if [ -n "${LSB_JOBINDEX:-}" ]; then
        SCRATCH_DIR="/tmp/${USER}/${LSB_JOBID}_${LSB_JOBINDEX}.scratch"
    else
        SCRATCH_DIR="/tmp/${USER}/${LSB_JOBID}.scratch"
    fi
fi


###############################################################################
# CHECK /tmp SPACE
###############################################################################

TMP_FREE=$(df --output=avail -k /tmp | tail -1)
REQUIRED=$((SCRATCH_SIZE_GB * 1024 * 1024))

if [ "$TMP_FREE" -lt "$REQUIRED" ]; then
    echo "ERROR: Not enough free space in /tmp for ${SCRATCH_SIZE_GB} GB scratch."
    echo "       Required: ${REQUIRED} KB"
    echo "       Available: ${TMP_FREE} KB"
    exit 1
fi


###############################################################################
# CREATE SCRATCH DIRECTORY
###############################################################################

mkdir -p "$SCRATCH_DIR"


###############################################################################
# CLEANUP HANDLER
###############################################################################

cleanup() {
    rm -rf "$SCRATCH_DIR" 2>/dev/null || true
}
trap cleanup EXIT


###############################################################################
# RUN WORKER SCRIPT
###############################################################################
JOBLIST="$1"

if [ -z "$JOBLIST" ]; then
    echo "ERROR: No joblist file passed."
    exit 1
fi

if [ ! -f "$JOBLIST" ]; then
    echo "ERROR: Joblist file not found: $JOBLIST"
    exit 1
fi

echo "LSB_JOBINDEX = $LSB_JOBINDEX"
echo "Using joblist: $JOBLIST"

line="$(sed -n "${LSB_JOBINDEX}p" "$JOBLIST")"

echo "Selected joblist line: $line"

# Convert the line into args safely (simple version):
set -- $line

echo "Passed args:"
printf '  [%s]\n' "$@"

#echo "Passed args:"
#echo "$@"

export SCRATCH_DIR
python /pi/summer.thyme-umw/enamine-REAL-2.6billion/umass_chan_REAL-M_platform/rosetta/chris_hull_optimizations/run_ligand_discovery_search.py "$@"

