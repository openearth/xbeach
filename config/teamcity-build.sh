#!/bin/bash
. /etc/profile
export MODULEPATH=$MODULEPATH:/opt/modules

module load gcc/4.9.1
module load hdf5/1.8.13_gcc_4.9.1
module load netcdf/v4.3.2_v4.4.0_gcc_4.9.1

./autogen.sh

FCFLAGS="-mtune=corei7-avx -funroll-loops --param max-unroll-times=4 -O3 -ffast-math" ./configure  --with-netcdf
make