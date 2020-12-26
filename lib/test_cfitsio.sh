#!/usr/bin/env bash

# Set the default result
CFITSIO_OK=1

# Write a simple test program using GSL
echo "
#include <fitsio.h>

int main() {
 return 0;
}
" > testprog.c

# check if it compiles
cc -o testprog testprog.c -lm -lcfitsio &>/dev/null
if [ $? -ne 0 ];then
 CFITSIO_OK=0
else
 # check if the executable file is produced
 if [ ! -x testprog ];then
  CFITSIO_OK=0
 else
  # checks if the program can be started and returns 0 exit code
  ./testprog
  if [ $? -ne 0 ];then
   CFITSIO_OK=0
  fi
 fi
fi

# remove temporary files generated by the test
for TEST_FILE_TO_REMOVE in testprog testprog.c ;do
 if [ -f "$TEST_FILE_TO_REMOVE" ];then
  rm -f "$TEST_FILE_TO_REMOVE"
 fi
done

# echo report results
if [ $CFITSIO_OK -eq 1 ];then
 echo "-lcfitsio"
else
 echo "-DPATPC_NOCFITSIO"
fi
