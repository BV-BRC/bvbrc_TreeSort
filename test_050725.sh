#!/bin/bash

# This will be called by the singularity container.

# Update the path
export PATH=$PATH:/project:/project/scripts:/project/external/TreeSort:/home/ac.ddempsey/miniconda3/envs/treesort-env/bin

# TEST
cd /project/external/TreeSort

# Run the Perl script that runs the Python module that runs TreeSort.
perl /project/service-scripts/App-TreeSort.pl https://p3.theseed.org/services/app_service /project/app_specs/TreeSort.json /project/tests/swine_H1/jobdesc.json

# Run the Python module that runs TreeSort.
#run_treesort -j ./tests/swine_H1_test/jobdesc.json -w ./tests/swine_H1_test

# Run TreeSort directly.
#./prepare_dataset.sh ../../tests/swine_H1/swine_H1_HANA.fasta HA ../../tests/swine_H1/test_output
