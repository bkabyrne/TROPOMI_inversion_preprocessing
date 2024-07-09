# TROPOMI_inversion_preprocessing


## Contents

- submit_job.sh
> Submits job.pbs to the compute node queue
- job.pbs
> PBS script to run download_and_process_TROPOMI
- download_and_process_TROPOMI.py
> Script that loops over days to call TROPOMI data downloading (bash script) and superobs scripts (python function)
- download_TROPOMI.sh
> Downloads TROPOMI data from an s3 bucket
- calc_TROPOMI_mole_fractions_2x25.py
> functions to calculate TROPOMI XCO obs
- create_folder_TROPOMI.py
> function to move TROPOMI data into new directory format (was only used once) 