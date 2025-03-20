#!/bin/bash

# Setup python environment
source /home/cbyers/Projects/openmc_rb/mc_PYenv/bin/activate
export PYTHONPATH=/home/cbyers/Projects/dagmc_bld/aegis-deps/MOAB/lib/python3.10/site-packages:$PYTHONPATH

# Setup library paths
export LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu/hdf5/openmpi/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/home/cbyers/Projects/dagmc_bld/aegis-deps/EMBREE/lib:$LD_LIBRARY_PATH

# Setup individual package directories
export HDF5_ROOT=/usr/lib/x86_64-linux-gnu/hdf5/openmpi/
export MOAB_DIR=/home/cbyers/Projects/dagmc_bld/aegis-deps/MOAB
export DAGMC_DIR=/home/cbyers/Projects/dagmc_bld/aegis-deps/DAGMC/bld/

# Setup C++ include paths
export CPLUS_INCLUDE_PATH=/usr/lib/x86_64-linux-gnu/openmpi/include/:$CPLUS_INCLUDE_PATH

# Change for different cross-section files
export OPENMC_CROSS_SECTIONS=/home/cbyers/Projects/crosssecs/endfb-vii.1-hdf5/cross_sections.xml
