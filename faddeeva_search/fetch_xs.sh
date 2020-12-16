#!/bin/bash

set -e

# ENDF/B v.ii data
if [[ ! -d ./endfb71_hdf5 ]] ; then
    wget -O- https://anl.box.com/shared/static/9igk353zpy8fn9ttvtrqgzvw1vtejoz6.xz | tar -Jx
fi

# ENDF/B v.iii data
if [[ ! -d ./endfb80_hdf5 ]] ; then
wget -O- https://anl.box.com/shared/static/uhbxlrx7hvxqw27psymfbhi7bx7s6u6a.xz | tar -Jx
fi
