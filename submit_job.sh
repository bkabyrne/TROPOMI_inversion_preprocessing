#!/bin/bash

# Set the PATH to include the directory where qstat and qsub are located
export PATH=/PBS/bin:$PATH

# Base name of the job
JOB_BASE_NAME="TROPOMI_preprocess"

# Check if any jobs with the base name are running
if /PBS/bin/qstat -u $USER | grep -q $JOB_BASE_NAME; then
    echo "A job with the base name $JOB_BASE_NAME is still running. Skipping submission."
else
    # Submit the PBS job
    /PBS/bin/qsub /nobackupp19/bbyrne1/TROPOMI_inversion_preprocessing/job.pbs
fi
