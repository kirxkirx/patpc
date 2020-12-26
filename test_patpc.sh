#!/usr/bin/env bash

# Check if the executable files are present, if not - run make
for EXE_FILE_TO_CHECK in simulator patpc ;do
 if [ ! -x "$EXE_FILE_TO_CHECK" ];then
  make
 fi
done
for EXE_FILE_TO_CHECK in simulator patpc ;do
 if [ ! -x "$EXE_FILE_TO_CHECK" ];then
  echo "ERROR cannot find executable file '$EXE_FILE_TO_CHECK'"
  exit 1
 fi
done

# Simulate photon arrival time data with 400 sec periodicity
echo "Simulating data..."
./simulator > test.dat 
if [ $? -ne 0 ];then
 echo "ERROR generating the test data! Maybe GSL is not installed?"
 exit 1 # assume things are fine
fi

# Run the period search
echo "Searching for a period..." 
TEST_RESULT=`./patpc test.dat 2>&1 | grep 'The peak Hm is at period' | awk 'function abs(x){return ((x < 0.0) ? -x : x)} {if (abs($7-400) > 1 ) print 1 ;else print 0}'`

# Remove files created by the test
for FILE_TO_REMOVE in binned_lightcurve_time.??? Hm.??? phase_folded_and_binned_lightcurve.??? power.??? test.dat ;do
 if [ -f "$FILE_TO_REMOVE" ];then
  rm -f "$FILE_TO_REMOVE"
 fi
done

if [ $TEST_RESULT -eq 0 ];then
 echo "Test passed"
 exit 0
else
 echo "Test failed"
 exit 1
fi
