#!/bin/sh
export TEAMCITY=yes
# we could update the configure script
#autoreconf
./configure --without-mpi
make
echo "##teamcity[testSuiteStarted name='nompi']"
make check
echo "##teamcity[testSuiteFinished name='nompi']"
make dist
mkdir nompi
mv ./xbeach ./*.tar.gz nompi
make clean
./configure --with-mpi
make
echo "##teamcity[testSuiteStarted name='nompi']"
make check
echo "##teamcity[testSuiteFinished name='nompi']"
make dist
mkdir mpi
mv ./xbeach ./*.tar.gz nompi
make clean
