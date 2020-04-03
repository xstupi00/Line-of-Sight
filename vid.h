/**************************************************************
 * File:		vid.h
 * Author:		Šimon Stupinský
 * University: 	Brno University of Technology
 * Faculty: 	Faculty of Information Technology
 * Course:	    Parallel and Distributed Algorithms
 * Date:		03.04.2020
 * Last change:	14.04.2020
 *
 * Subscribe:	TODO
 *
**************************************************************/

/**
 * @file    ots.h
 * @brief   TODO
 */

#include <cmath>
#include <cstring>
#include <vector>

#include <mpi.h>

// The rank of the master processor
#define MASTER 0
// Number of elements send in the MPI_Send and MPI_Recv
#define COUNT 1
// Message tag in the MPI_Send and MPI_Recv
#define TAG 0
#define TAG_1 1
// Flag for measuring the elapsed time in the sorting algorithm
#define MEASURE_TIME

using namespace std;