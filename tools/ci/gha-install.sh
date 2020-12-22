#!/bin/bash
set -ex

# Upgrade pip, pytest, numpy before doing anything else
pip install --upgrade pip
pip install --upgrade pytest
pip install --upgrade numpy

# Install mpi4py for MPI configurations
if [[ $MPI == 'y' ]]; then
    pip install --no-binary=mpi4py mpi4py
fi

# Build and install OpenMC executable
python tools/ci/gha-install.py

# Install Python API in editable mode
pip install -e .[test,vtk]

# For coverage testing of the C++ source files
pip install cpp-coveralls

# For coverage testing of the Python source files
pip install coveralls
