#!/bin/bash

# First test that OpenFOAM commands are available
which surfaceFeatureExtract &> /dev/null
retVal=$?
if [ $retVal -ne 0 ]; then
    echo "Could not find OpenFOAM command surfaceFeatureExtract, exiting. OpenFOAM is probably not sourced correctly?"
    exit 1
fi

# Stop at errors
set -e

# Clean up test folder and initialize
test_folder="./run_tests"
echo "Removing $test_folder"
rm -rf $test_folder
echo "Populating $test_folder"
mkdir $test_folder
cp -rv testcase* $test_folder &> /dev/null
cd $test_folder

# Function for the tests
function do_test() {
  cd $1
  folder=`pwd`

  echo "Running run_parallel at $folder"
  ./run_parallel
  echo "grepped checkMesh errors (if any):"
  mpirun -np 3 checkMesh -parallel -latestTime

  echo "Running run_serial at $folder"
  ./run_serial
  echo "grepped checkMesh errors (if any):"
  checkMesh -latestTime

  cd ..
}

# Tests
do_test "testcase"
do_test "testcase2"
do_test "testcase3"
do_test "testcase4"
do_test "testcase5"

echo "Test runs completed"
cd ..
