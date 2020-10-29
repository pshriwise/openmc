
#!/bin/bash
set -ex

# MOAB Variables
MOAB_BRANCH='Version5.1.0'
MOAB_REPO='https://bitbucket.org/fathomteam/moab/'
MOAB_INSTALL_DIR=$HOME/MOAB/

# DAGMC Variables
DAGMC_BRANCH='develop'
DAGMC_REPO='https://github.com/svalinn/dagmc'
DAGMC_INSTALL_DIR=$HOME/DAGMC/

CURRENT_DIR=$(pwd)

# MOAB Install
conda install -c conda-forge DAGMC

cd $CURRENT_DIR
