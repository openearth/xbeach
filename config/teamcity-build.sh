#!/bin/bash
. /etc/profile
export MODULEPATH=$MODULEPATH:/opt/modules
export LD_LIBRARY_PATH=/opt/openmpi/1.8.1_gcc_4.9.1/lib:$LD_LIBRARY_PATH
export PATH=/opt/openmpi/1.8.1_gcc_4.9.1/bin:$PATH

module load gcc/4.9.1
module load hdf5/1.8.13_gcc_4.9.1
module load netcdf/v4.3.2_v4.4.0_gcc_4.9.1

./autogen.sh

mkdir -p /opt/teamcity/work/XBeach_unix/install

FCFLAGS="-mtune=corei7-avx -funroll-loops --param max-unroll-times=4 -O3 -ffast-math" ./configure  --with-netcdf --prefix="/opt/teamcity/work/XBeach_unix/install"

make
make install