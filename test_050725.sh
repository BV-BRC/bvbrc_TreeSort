#!/bin/bash

# This will be run by the singularity container to:
# 1) Update the path to include the mount point for the bvbrc_TreeSort_dmd directory.
# 2) Change to that directory.
# 3) Run the treesortrunner.run_treesort Python module.

# Update the path
export PATH=$PATH:/project:/project/scripts:/project/external/TreeSort:/home/ac.ddempsey/miniconda3/bin/

# Go to the scripts directory.
cd /project/scripts

# Run the Python module that runs TreeSort.
python -m treesortrunner.run_treesort -j ./test/test_job.json
