#!/bin/bash
set -ex

source "$HOME/miniconda/etc/profile.d/conda.sh"
conda activate openmc

# Argument List
args=" "

# Check for MPI
if [[ $MPI == 'y' ]]; then
  args="${args} --mpi "
fi
  
# Check for event-based
if [[ $EVENT == 'y' ]]; then
  args="${args} --event "
fi

# Run regression and unit tests
pytest --cov=openmc -v $args tests
