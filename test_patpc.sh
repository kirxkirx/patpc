#!/usr/bin/env bash

# Check if the executable files are present, if not - run make
for EXE_FILE_TO_CHECK in simulator patpc swift_pointing_evt_splitter ;do
 if [ ! -x "$EXE_FILE_TO_CHECK" ];then
  make
 fi
done
for EXE_FILE_TO_CHECK in simulator patpc swift_pointing_evt_splitter ;do
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
TEST_RESULT=$(./patpc test.dat 2>&1 | grep 'The peak Hm is at period' | awk 'function abs(x){return ((x < 0.0) ? -x : x)} {if (abs($7-400) > 1 ) print 1 ;else print 0}')

# Remove files created by the test
for FILE_TO_REMOVE in binned_lightcurve_time.??? Hm.??? phase_folded_and_binned_lightcurve.??? power.??? test.dat ;do
 if [ -f "$FILE_TO_REMOVE" ];then
  rm -f "$FILE_TO_REMOVE"
 fi
done

if [ $TEST_RESULT -ne 0 ];then
 echo "Period search test failed"
 exit 1
fi
echo "Period search test passed"

# Test swift_pointing_evt_splitter
echo ""
echo "Testing swift_pointing_evt_splitter..."

# Create test data with deliberate gaps > 1000 seconds
# Three chunks: 0-100s, 2000-2100s, 5000-5100s
cat > test_splitter.dat <<EOF
0.0
10.0
20.0
30.0
40.0
50.0
60.0
70.0
80.0
90.0
100.0
2000.0
2010.0
2020.0
2030.0
2040.0
2050.0
2060.0
2070.0
2080.0
2090.0
2100.0
5000.0
5010.0
5020.0
5030.0
5040.0
5050.0
5060.0
5070.0
5080.0
5090.0
5100.0
EOF

# Run the splitter
echo "Running swift_pointing_evt_splitter on test data with gaps..."
./swift_pointing_evt_splitter test_splitter.dat > /dev/null 2>&1

# Check if the expected output files were created
EXPECTED_FILES=3
FOUND_FILES=0
for i in 01 02 03 ;do
 if [ -f "test_splitter_${i}.dat" ];then
  FOUND_FILES=$((FOUND_FILES + 1))
 fi
done

# Clean up test files
rm -f test_splitter.dat test_splitter_01.dat test_splitter_02.dat test_splitter_03.dat

if [ $FOUND_FILES -eq $EXPECTED_FILES ];then
 echo "swift_pointing_evt_splitter test passed (created $FOUND_FILES chunk files as expected)"
 echo ""
 echo "All tests passed"
 exit 0
else
 echo "swift_pointing_evt_splitter test failed (expected $EXPECTED_FILES files, found $FOUND_FILES)"
 exit 1
fi
