#!/bin/bash

if [ "$#" -eq 1 ];
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

# shellcheck disable=SC2006
PROCESSORS=`python3 -c "from scipy import optimize; import scipy; print(int(scipy.floor(scipy.optimize.fsolve(lambda x: $POINTS_NUMBER - scipy.log2(x) * x, 5)[0])))"`
# shellcheck disable=SC2004
PROCESSORS=$(($POINTS_NUMBER))

# compile source code
mpic++ --prefix /usr/local/share/OpenMPI -o vid vid.cpp

# run binary
mpirun --prefix /usr/local/share/OpenMPI --hostfile hostfile -np "$PROCESSORS" vid "$LINE_OF_SIGHT"

# clean directory
rm -f vid