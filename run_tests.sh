#!/bin/bash

# First test that OpenFOAM commands are available
which surfaceFeatureExtract &> /dev/null
retVal=$?
if [ $retVal -ne 0 ]; then
    echo "Could not find OpenFOAM command surfaceFeatureExtract, exiting. OpenFOAM v2312 is probably not sourced correctly?"
    exit 1
fi

# Clean up test folder
test_folder="./run_tests"
echo "Removing $test_folder"
rm -rf $test_folder
echo "Populating $test_folder"
mkdir $test_folder
cp -rv testcase* $test_folder &> /dev/null

cd $test_folder/testcase
folder=`pwd`
echo "Running run_parallel at $folder"
./run_parallel | grep Smoothing\ iteration
echo "Running run_serial at $folder"
./run_serial | grep Smoothing\ iteration
cd ..

echo "Test runs completed"
cd ..
