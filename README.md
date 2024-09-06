# README


**Preprocessing scripts to prepare TROPOMI carbon monoxide for assimilation**

**contact: Brendan Byrne**

**email: brendan.k.byrne@jpl.nasa.gov**

---

Copyright 2024, by the California Institute of Technology. ALL RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any commercial use must be negotiated with the Office of Technology Transfer at the California Institute of Technology.
 
This software may be subject to U.S. export control laws. By accepting this software, the user agrees to comply with all applicable U.S. export laws and regulations. User has the responsibility to obtain export licenses, or other export authority as may be required before exporting such information to foreign countries or providing access to foreign persons.

## Overview

This software automates the downloading of TROPOMI carbon monoxide data, computes 2 x 2.5 degree super-obs, converts data format to one assimilable my CMS-Flux, and writes the data to the format required by CMS-Flux

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