#!/bin/bash

################################################################
# File:       test.sh
# Author:		  Šimon Stupinský
# University: Brno University of Technology
# Faculty: 	  Faculty of Information Technology
# Course:	    Parallel and Distributed Algorithms
# Date:		    04.04.2020
# Last change:08.04.2020
#
# Subscribe:	The test shell script to run algorithm that solve
#             parallel Line-Of-Sight problem with use max-prescan.
#
################################################################

#####
# @file    test.sh
# @brief   The script firstly check whether the required input was
#          entered on the command line. When was entered also the
#          second argument then is computed the relevant number
#          of processors according to the numbers of the altitudes
#          in the given line-of-sight. Then this script builds the
#          program and subseqenlty runs with the computed numbers
#          of processors. In the end it removes all created files
#          during its execution.
#####

# check whether the input line of sight was entered, otherwise error
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

# comput the numbers of processor according to given second argument or by default
if  [ "$#" -eq 2 ]
then
  # compute the number of processor according to formula n/p >= lg p
  if [ "$2" -eq 1 ]
  then
    PROCESSORS=$(python3 -c "from scipy import optimize; import scipy; import numpy; print(int(numpy.ceil(scipy.optimize.fsolve(lambda x: $POINTS_NUMBER - numpy.lib.scimath.log2(x) * x, 5)[0])))")
  # number of processor is equal to n/2
  elif [ "$2" -eq 2 ]
  then
    PROCESSORS=$(python3 -c "from math import ceil; print(ceil($POINTS_NUMBER/2.0))")
  # number of processor is equal to n
  elif [ "$2" -eq 3 ]
  then
    PROCESSORS=$(python3 -c "from math import ceil; print(ceil($POINTS_NUMBER))")
  # number of processor is equal to 1
  else
    PROCESSORS=1
  fi
# default value is compute according to the rule n/p >= lg p
else
  PROCESSORS=$(python3 -c "from scipy import optimize; import scipy; import numpy; print(int(numpy.ceil(scipy.optimize.fsolve(lambda x: $POINTS_NUMBER - numpy.lib.scimath.log2(x) * x, 5)[0])))")
fi

# compile source code
mpic++ --prefix /usr/local/share/OpenMPI -o vid vid.cpp

# run binary
mpirun --prefix /usr/local/share/OpenMPI -np "$PROCESSORS" vid "$LINE_OF_SIGHT"

# clean directory
rm -f vid