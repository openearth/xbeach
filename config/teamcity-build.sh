#!/bin/bash
. /etc/profile
export MODULEPATH=$MODULEPATH:/opt/modules

cd ..
export SVNPATH=$(pwd)
cd trunk/

# start with anaconda since it overwrites gcc !!!
module load anaconda/lnx64_conda

module load gcc/4.9.2
module load hdf5/1.8.14_gcc_4.9.2
module load netcdf/v4.3.2_v4.4.0_gcc_4.9.2
module load openmpi/1.8.3_gcc_4.9.2

make distclean

./autogen.sh

mkdir -p $SVNPATH"/install"

FCFLAGS="-mtune=corei7-avx -funroll-loops --param max-unroll-times=4 -ffree-line-length-none -O3 -ffast-math" ./configure  --with-netcdf --with-mpi --prefix=$SVNPATH"/install"

make
make install

cd $SVNPATH"/install"

tar -cf $SVNPATH"/artifacts.tar" bin lib

cd $SVNPATH"/trunk/config/"

/usr/share/Modules/bin/createmodule.py -p "/opt/xbeach/"$XBEACH_PROJECT_ID"_gcc_4.9.2_1.8.3_HEAD" ./teamcity-env.sh > $SVNPATH"/trunk/config/xbeach-"$XBEACH_PROJECT_ID"_gcc_4.9.2_1.8.3_HEAD"