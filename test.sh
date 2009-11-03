#!/bin/sh
testname=`echo $0 | sed -e 's/.*test_\(.*\)/\1/g'`
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
echo running test in `pwd`/run/$testname
cd run/$testname
# use the default parameters
mv params_default.txt params.txt
# cut it short
sed -i -e 's/tstop.*/tstop=20/' params.txt
sed -i -e 's/\(tint.*=\)\(.*\)/\1 1/' params.txt
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
    

