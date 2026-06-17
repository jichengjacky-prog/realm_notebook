#!/bin/bash
# Snakemake cluster-generic status script for LSF
# Returns: "running" if job is PEND or RUN, "success" if DONE, "failed" otherwise

JOBID="$1"

if [ -z "$JOBID" ]; then
    echo "failed"
    exit 0
fi

STATUS=$(bjobs -o stat -noheader "$JOBID" 2>/dev/null | tr -d ' ')

case "$STATUS" in
    RUN|PEND|WAIT|PSUSP|USUSP|SSUSP)
        echo "running"
        ;;
    DONE)
        echo "success"
        ;;
    *)
        echo "failed"
        ;;
esac
