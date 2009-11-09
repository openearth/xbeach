#!/bin/sh
export TEAMCITY=yes
# we could update the configure script
#autoreconf
echo "##teamcity[progressMessage 'configuring nompi version']"
./configure --without-mpi
echo "##teamcity[progressMessage 'building source distribution']"
make dist
echo "##teamcity[progressMessage 'building xbeach']"
make
echo "##teamcity[testSuiteStarted name='nompi']"
make check
echo "##teamcity[testSuiteFinished name='nompi']"
mkdir nompi
mv ./xbeach nompi
echo "##teamcity[progressMessage 'cleaning']"
make clean
echo "##teamcity[progressMessage 'configuring mpi version']"
./configure --with-mpi
echo "##teamcity[progressMessage 'building xbeach mpi version']"
make
echo "##teamcity[testSuiteStarted name='nompi']"
make check
echo "##teamcity[testSuiteFinished name='nompi']"
mkdir mpi
mv ./xbeach nompi
echo "##teamcity[progressMessage 'cleaning']"
make clean
echo "##teamcity[progressMessage 'done']"