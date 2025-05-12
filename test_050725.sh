#!/bin/bash

# This will be called by the singularity container.

# Update the path
export PATH=$PATH:/project:/project/scripts:/project/external/TreeSort:/home/ac.ddempsey/miniconda3/bin/

# TEST
# Update the Python path
export PYTHONPATH=$PYTHONPATH:/project/scripts/treesortrunner

# Go to the scripts directory.
cd /project/scripts

# Run the Python module that runs TreeSort.
python -m treesortrunner.run_treesort -j ./test/swine_H1_test/test_job.json -w ./test/swine_H1_test
