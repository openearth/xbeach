#!/bin/sh
# check if we need to go into the trunk
if [ -d trunk ]
then
    cd trunk
fi
export TEAMCITY=yes
# we could update the configure script
#autoreconf
echo "##teamcity[progressMessage 'configuring nompi version']"
./configure --without-mpi
echo "##teamcity[progressMessage 'cleaning up']"
make clean
echo "##teamcity[progressMessage 'building source distribution']"
make dist
echo "##teamcity[progressMessage 'building xbeach']"
make
echo "##teamcity[testSuiteStarted name='nompi']"
make check
echo "##teamcity[testSuiteFinished name='nompi']"
mkdir nompi
mv ./xbeach nompi
