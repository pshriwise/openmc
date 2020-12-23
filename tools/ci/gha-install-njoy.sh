#!/bin/bash
set -ex

cd $HOME
git clone https://github.com/njoy/NJOY2016
cd NJOY2016
mkdir install
mkdir build && cd build
cmake -Dstatic=on -DCMAKE_INSTALL_PREFIX=$HOME/NJOY2016/install .. && make 2>/dev/null && sudo make install
