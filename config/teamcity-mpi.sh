#!/bin/sh
export TEAMCITY=yes
# we could update the configure script
#autoreconf
echo "##teamcity[progressMessage 'configuring mpi version']"
./configure --with-mpi
echo "##teamcity[progressMessage 'building xbeach mpi version']"
make
echo "##teamcity[testSuiteStarted name='mpi']"
make check
echo "##teamcity[testSuiteFinished name='mpi']"
mkdir mpi
mv ./xbeach mpi
echo "##teamcity[progressMessage 'cleaning']"
make clean
echo "##teamcity[progressMessage 'done']"