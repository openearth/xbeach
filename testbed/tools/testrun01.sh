#!/bin/sh
testname=run01
xbeachdir=`pwd`
# goto the testbed tree
cd testbed
# cleanup old run, may not work
rm -rf run/$testname &> /dev/null
# create a new run directory
mkdir -p run/$testname
# copy input
cp -r input/$testname run
# go into the test
cd run/$testname
# run xbeach
$xbeachdir/xbeach > output
# did it run succesful?
status=$?
if [ $status != 0 ] 
then
    # nope
    # don't cleanup, so we can check the output later...
    exit $status
else
    # cleanup
    cd $xbeachdir
    rm -rf testbed/run/$testname
    exit 0
fi
    

