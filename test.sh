#!/bin/bash

if [ "$#" -eq 1 ] || [ "$#" -eq 2 ];
then
  LINE_OF_SIGHT=$1
else
  echo "Missing string LINE_OF_SIGHT! Usage: ./test.sh LINE_OF_SIGHT"
  exit 1
fi

IFS=',' # hyphen (-) is set as delimiter
read -ra ARR <<< "$LINE_OF_SIGHT" # str is read into an array as tokens separated by IFS
IFS=' '
POINTS_NUMBER=${#ARR[@]}

if  [ "$#" -eq 2 ]
then
  if [ "$2" -eq 1 ]
  then
    PROCESSORS=$(python3 -c "from scipy import optimize; import scipy; import numpy; print(int(numpy.ceil(scipy.optimize.fsolve(lambda x: $POINTS_NUMBER - numpy.lib.scimath.log2(x) * x, 5)[0])))")
  elif [ "$2" -eq 2 ]
  then
    PROCESSORS=$(python3 -c "from math import ceil; print(ceil($POINTS_NUMBER/2.0))")
  elif [ "$2" -eq 3 ]
  then
    PROCESSORS=$(python3 -c "from math import ceil; print(ceil($POINTS_NUMBER))")
  else
    PROCESSORS=1
  fi
else
  PROCESSORS=$(python3 -c "from scipy import optimize; import scipy; import numpy; print(int(numpy.ceil(scipy.optimize.fsolve(lambda x: $POINTS_NUMBER - numpy.lib.scimath.log2(x) * x, 5)[0])))")
fi

# compile source code
mpic++ --prefix /usr/local/share/OpenMPI -o vid vid.cpp

# run binary
mpirun --prefix /usr/local/share/OpenMPI --hostfile hostfile -np "$PROCESSORS" vid "$LINE_OF_SIGHT"

# clean directory
rm -f vid