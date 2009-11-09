#!/bin/sh
export TEAMCITY=yes
./configure --without-mpi
make
make check
make dist
mkdir nompi
mv ./xbeach ./*.tar.gz nompi
make clean
./configure --with-mpi
make
make check
make dist
mkdir mpi
mv ./xbeach ./*.tar.gz nompi
make clean
