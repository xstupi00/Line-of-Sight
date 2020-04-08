/**************************************************************
 * File:		vid.h
 * Author:		Šimon Stupinský
 * University: 	Brno University of Technology
 * Faculty: 	Faculty of Information Technology
 * Course:	    Parallel and Distributed Algorithms
 * Date:		03.04.2020
 * Last change:	07.04.2020
 *
 * Subscribe:	The header module of the program implementing Line-of-Sight problem.
 *
**************************************************************/

/**
 * @file    vid.h
 * @brief   The header module contains definitions of the program constants
 *          and macros and also declarations of the functions.
 */

#include <climits>
#include <cmath>
#include <limits>
#include <mpi.h>
#include <vector>

// The rank of the master processor
#define MASTER 0
// Number of elements send in the MPI_Bcast message
#define COUNT 1
// Flag for measuring the elapsed time in the sorting algorithm
//#define MEASURE_TIME
// Size of the integer variable to use as the size of the elements within the shared memory
#define INT_UNIT sizeof(int)
// Size of the float (double) variable to use as the size of the elements within the shared memory
#define FLOAT_UNIT sizeof(float)
// Size of the boolean variable to use as the size of the elements within the shared memory
#define BOOL_UNIT sizeof(bool)
// Macro that finds the nearest power of two according to the given number x
#define NEXT_POWER_2(x) ((x == 1) ? 1 : ((1ULL << sizeof(x) * CHAR_BIT) >> __builtin_clz(x - 1)))
// Define the minimum of the float numbers to use as the neutral element (I) within max-prescan operation
#define FLOAT_MIN -std::numeric_limits<float>::max()

using namespace std;

/**
 * Loads the input line of sight in the specified format and individual numbers
 * convert to the integers and subsequently it stores them to the given vector.
 * The input is in the following format: x_1,x_2,...,x_n, where x_i in N for 1 <= i <= n
 * The result vector will contain the processed numbers: <x1, x2, ..., xn>.
 *
 * @param input_altitudes   input char sequence represents the points in the terrain (altitudes)
 * @param target_altitudes  output vector of the integer points in the terrain (altitudes)
 */
void load_line_of_sight(char *input_altitudes, std::vector<int> *target_altitudes);

/**
 * Computes the number of points which will be processed by the one process.
 * Allocates the shared memory for all process where will be stored the loaded
 * sequence of the altitudes. The master process copies the given vector to
 * the allocated shared memory.
 *
 * @param shared_altitudes  shared window object used for communication (initial address) - handle
 * @param node_altitudes    window node of each process - initial address of the process window (choice)
 * @param window_size       size of the process window in bytes (non-negative integer)
 * @param total_altitudes   number of the altitudes available within all processes
 * @param size              size of a communicator MPI_COMM_WORLD (number of the processes)
 * @param altitudes         vector of the loaded altitudes
 * @param rank              rank of the process in a communicator MPI_COMM_WORLD
 */
void share_points_to_process(
    int **shared_altitudes, MPI_Win *node_altitudes, MPI_Aint *window_size, int total_altitudes, int size,
    std::vector<int> altitudes, int rank
);

/**
 * Computes the relevant index of each process, where start its window with respect
 * to the whole shared allocated window. Queries the size and base pointer for a
 * patch of shared memory windows to obtain all required data for computation.
 * Each process computes own part of the angles from the given altitudes and stores
 * them at the relevant index to the resulting shared window.
 *
 * @param shared_altitudes      shared allocated window containing the loaded altitudes
 * @param shared_angles         shared allocated window to store newly computed angles from altitudes
 * @param max_previous_angles   shared allocated window to store newly computed angles from altitudes (for max-prescan)
 * @param node_altitudes        window object of each process to the shared allocated window of stored altitudes
 * @param node_angles           window object of each process to the shared allocated window of computed angles
 * @param node_prev_angles      window object of each process to the shared allocated window of previous angles
 * @param window_size           size of the process window in bytes (non-negative integer)
 * @param disp_unit             local unit size for displacements, in bytes (non-negative integer)
 * @param total_altitudes       number of the altitudes available within all processes
 * @param size                  size of a communicator MPI_COMM_WORLD (number of the processes)
 * @param rank                  rank of the process in a communicator MPI_COMM_WORLD
 */
void compute_angles(
    int **shared_altitudes, float **shared_angles, float **max_previous_angles, MPI_Win node_altitudes, MPI_Win
    node_angles, MPI_Win node_prev_angles, MPI_Aint window_size, int disp_unit, int total_altitudes, int size, int rank
);

/**
 * Performs the max-prescan operation on the given vector of the computed angles. Firstly,
 * it performs the up-sweep phase of this operation, then the master stores the neutral
 * element to the root of the processing tree and then are performed down-sweep phase.
 *
 * @param shared_angles         shared allocated window containing the computed angles from altitudes
 * @param node_angles           window object of each process to the shared allocated window of previous angles
 * @param window_size           size of the process window in bytes (non-negative integer)
 * @param disp_unit             local unit size for displacements, in bytes (non-negative integer)
 * @param total_angles          number of the altitudes available within all processes
 * @param size                  size of a communicator MPI_COMM_WORLD (number of the processes)
 * @param rank                  rank of the process in a communicator MPI_COMM_WORLD
 */
void max_prescan(
    float **shared_angles, MPI_Win node_angles, MPI_Aint window_size, int disp_unit, int total_angles, int size, int rank
);

/**
 *  Each processor first find the maximum within self n/p section of the angles vector to
 *  generate a processor maximum, then the tree technique is used to max-prescan the
 *  processor sums. The results of the max-prescan of the processor maximums are used
 *  as an offset for each processor to prescan within its n/p section.
 *
 * @param shared_angles         shared allocated window containing the computed angles from altitudes
 * @param node_angles           window object of each process to the shared allocated window of previous angles
 * @param window_size           size of the process window in bytes (non-negative integer)
 * @param disp_unit             local unit size for displacements, in bytes (non-negative integer)
 * @param total_angles          number of the altitudes available within all processes
 * @param size                  size of a communicator MPI_COMM_WORLD (number of the processes)
 * @param rank                  rank of the process in a communicator MPI_COMM_WORLD
 */
void preprocess_subsets(
    float **shared_angles, MPI_Win node_angles, MPI_Aint window_size, int disp_unit, int total_angles, int size, int rank
);

/**
 * Computes the relevant index of each process, where start its window with respect
 * to the whole shared allocated window. Queries the size and base pointer for a
 * patch of shared memory windows to obtain all required data for computation.
 * Each process computes own part of the results from the given data and stores
 * them to the final shared allocated window at the relevant index.
 *
 * @param shared_angles         shared allocated window containing the computed angles from altitudes
 * @param max_previous_angles   shared allocated window containing the computed previous angles with max-prescan
 * @param result                shared allocated window to store the results of the subtraction of angles and previous
 * @param node_angles           window object of each process to the shared allocated window of computed angles
 * @param node_prev_angles      window object of each process to the shared allocated window of previous angles
 * @param node_results          window object of each process to the shared allocated window of final results
 * @param disp_unit             local unit size for displacements, in bytes (non-negative integer)
 * @param window_size           size of the process window in bytes (non-negative integer)
 * @param total_angles          number of the altitudes available within all processes
 * @param size                  size of a communicator MPI_COMM_WORLD (number of the processes)
 * @param rank                  rank of the process in a communicator MPI_COMM_WORLD
 */
void compute_results(
    float **shared_angles, float **max_previous_angles, bool **result, MPI_Win node_angles, MPI_Win node_prev_angles,
    MPI_Win node_results, int disp_unit, int window_size, int total_angles, int size, int rank
);

/**
 * The master process write out to the standard output the final results based on the
 * given shared allocated window of the computed results.
 *
 * @param result            shared allocated window of stored results of the subtraction of angles and previous angles
 * @param node_results      window object of each process to the shared allocated window of final results
 * @param window_size       size of the process window in bytes (non-negative integer)
 * @param disp_unit         local unit size for displacements, in bytes (non-negative integer)
 * @param total_points      number of the altitudes available within all processes
 */
void write_out_result(bool **result, MPI_Win node_results, MPI_Aint window_size, int disp_unit, int total_altitudes);