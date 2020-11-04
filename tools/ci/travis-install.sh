#!/bin/bash
set -ex

# Install conda packages
conda install eigen
conda install mpich
conda install -c conda-forge fortran-compiler

if [[ $MPI = 'y' ]]; then
    conda install mpich mpi4py "h5py=*=*mpich*"
else
    conda install h5py
fi

# Python dependencies
conda install pip

# Install NJOY 2016
./tools/ci/travis-install-njoy.sh

# Install DAGMC if needed
if [[ $DAGMC = 'y' ]]; then
    conda install -c conda-forge MOAB
    ./tools/ci/travis-install-dagmc.sh
fi

# Install vectfit for WMP generation if needed
if [[ $VECTFIT = 'y' ]]; then
    ./tools/ci/travis-install-vectfit.sh
fi

# Build and install OpenMC executable
python tools/ci/travis-install.py

# Install Python API in editable mode
pip install -e .[test,vtk]

# For coverage testing of the C++ source files
pip install cpp-coveralls

# For coverage testing of the Python source files
pip install coveralls
