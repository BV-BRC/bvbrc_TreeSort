#!/bin/bash

# This script initializes the environment for BV-BRC module development.

# Source conda shell integration manually.
if [ -f ~/miniconda3/etc/profile.d/conda.sh ]; then
    source ~/miniconda3/etc/profile.d/conda.sh
else
    echo "Error: Conda shell integration script not found."
    exit 1
fi

# Activate the treesort conda environment.
if ! conda activate treesort-env; then
    echo "Error: Failed to activate Conda environment 'treesort-env'."
    exit 1
fi

# Add all necessary mount points.
mountPoints="/vol,/home,/homes,/disks/patric-common/"
mountPoints+=",/disks/tmp:/tmp"
mountPoints+=",/vol/patric3/production/data-images/patric-data-2022-0119:/opt/patric-common/data"
mountPoints+=",/home/ac.ddempsey/dev/bvbrc_TreeSort_dmd:/project"

# The Singularity container.
container="/vol/patric3/production/containers/ubuntu-063-12.sif"

# Check if the container exists.
if [ ! -f "$container" ]; then
    echo "Error: Singularity container not found at $container."
    exit 1
fi

# Make sure the test script exists.
if [ ! -f "test_050725.sh" ]; then
    echo "Error: Can't find script test_050725.sh."
    exit 1
fi

# Start the Singularity container.
singularity shell -B $mountPoints $container || {
    echo "An error occurred in the Singularity container."
    exit 1
}
