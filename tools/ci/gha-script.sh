#!/bin/bash
set -ex

# Argument List
args=" "

# Check for MPI
if [[ $MPI == 'mpi-y' ]]; then
  args="${args} --mpi "
fi

# Check for event-based
if [[ $EVENT == 'event-y' ]]; then
  args="${args} --event "
fi

# Run regression and unit tests
pytest --cov=openmc -v $args tests
