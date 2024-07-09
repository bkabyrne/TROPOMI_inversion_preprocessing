#!/bin/bash

# Base name of the job
JOB_BASE_NAME="TROPOMI_preprocess"

# Check if any jobs with the base name are running
if qstat -u $USER | grep -q $JOB_BASE_NAME; then
    echo "A job with the base name $JOB_BASE_NAME is still running. Skipping submission."
else
    # Submit the PBS job
    qsub job.pbs
fi
