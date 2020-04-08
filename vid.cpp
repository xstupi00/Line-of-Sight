/**************************************************************
 * File:		vid.cpp
 * Author:		Šimon Stupinský
 * University: 	Brno University of Technology
 * Faculty: 	Faculty of Information Technology
 * Course:	    Parallel and Distributed Algorithms
 * Date:		04.04.2020
 * Last change:	08.04.2020
 *
 * Subscribe:	The main module of the program implementing Line-of-Sight problem.
 *
**************************************************************/

/**
 * @file    vid.cpp
 * @brief   This module contains the implementation of the problem Line-of-Sight
 *          in its parallel version with use the max-prescan operation.
 */


#include "vid.h"


void load_line_of_sight(char *input_altitudes, std::vector<int> *target_altitudes) {
    // Break input line of sight into a series of tokens using the delimiter ',' - individual altitudes
    char *point = strtok(input_altitudes, ",");
    // Store the first altitude to the given vector of altitudes
    (*target_altitudes).push_back(atoi(point));
    // Walk through other altitudes from the input sequence
    while (point) {
        // Obtain the next one altitude or the NULL when the all altitudes was loaded
        point = strtok(NULL, ",");
        // When there is the newly altitude, then store it to the given vector of altitudes
        if (point) {
            (*target_altitudes).push_back(atoi(point));
        }
    }
}

void share_points_to_process(
    int **shared_altitudes, MPI_Win *node_altitudes, MPI_Aint *window_size, int total_altitudes, int size,
    std::vector<int> altitudes, int rank
) {
    // Compute the number of altitudes per each processor - equal number for each processor
    int points_per_processors = floor(total_altitudes / size);
    // Define the final size of the window for each process to each shared memory
    *window_size = (rank < (total_altitudes % size)) ? points_per_processors + 1 : points_per_processors;
    // Create an window for one-sided communication and shared memory access, and allocate memory at each process
    // Allocate shared window to store the loaded angles.
    MPI_Win_allocate_shared(
            (*window_size) * INT_UNIT, INT_UNIT, MPI_INFO_NULL, MPI_COMM_WORLD, &(*shared_altitudes), &(*node_altitudes)
    );
    // Master process copy the angles vector to the allocate shared memory for each process
    if (rank == MASTER) {
        std::copy(altitudes.begin(), altitudes.end(), *shared_altitudes);
    }
}

void compute_angles(
    int **shared_altitudes, float **shared_angles, float **max_previous_angles, MPI_Win node_altitudes, MPI_Win
    node_angles, MPI_Win node_prev_angles, MPI_Aint window_size, int disp_unit, int total_altitudes, int size, int rank
) {
    // Declare auxiliary window size to query at shared memory without the change relevant window size of processes
    MPI_Aint master_window_size;
    // Compute the start index for each process within the shared memory with respect to the common begin
    int start_idx = (rank < total_altitudes % size) ? rank * window_size : rank * window_size + (total_altitudes % size);
    // Query the size and base pointer for a patch of a shared memory window with shared altitudes
    MPI_Win_shared_query(node_altitudes, MASTER, &master_window_size, &disp_unit, &(*shared_altitudes));
    // Query the size and base pointer for a patch of a shared memory window with shared angles
    MPI_Win_shared_query(node_angles, MASTER, &master_window_size, &disp_unit, &(*shared_angles));
    // Query the size and base pointer for a patch of a shared memory window with maximum previous angles (max-prescan)
    MPI_Win_shared_query(node_prev_angles, MASTER, &master_window_size, &disp_unit, &(*max_previous_angles));
    // Each process compute own n/p sub-part of the whole angles and store them to the relevant index
    for (size_t i = start_idx; i < start_idx + window_size; i++) {
        // angle[i] = arctan( (altitude[i] - altitude[0]) / i ); neutral item for i=0
        (*shared_angles)[i] = (*max_previous_angles)[i] =
                (i) ? atan(((*shared_altitudes)[i] - (*shared_altitudes)[0]) / float(i)) : FLOAT_MIN;
    }
    // Blocks until all processes in the communicator have computed own part of the angles
    MPI_Barrier(MPI_COMM_WORLD);
}

void max_prescan(
    float **shared_angles, MPI_Win node_angles, MPI_Aint window_size, int disp_unit, int total_angles, int size, int rank
) {
    // Declare auxiliary window size to query at shared memory without the change relevant window size of processes
    MPI_Aint master_window_size;
    // Up-Sweep phase of the max-prescan operation
    // for d from 0 to (lg n) - 1
    for (size_t d = 0; d < ceil(log2(total_angles)); d++) {
        // in parallel for i from 0 to n-1 by 2^(d+1)
        if (!(rank % (1 << (d)))) {
            // Query the size and base pointer for a patch of a shared memory window with shared angles
            MPI_Win_shared_query(node_angles, MASTER, &master_window_size, &disp_unit, &(*shared_angles));
            // Obtain the value angles[i + 2(d) - 1] if there exists, otherwise replace value by FLOAT_MIN in 1.iter
            float left_node = (rank * 2 + (1 << (d)) - 1 >= total_angles and d == 0) ?
                              FLOAT_MIN : (*shared_angles)[rank * 2 + (1 << (d)) - 1];
            // Obtain the value angles[i + 2(d+1) - 1] if there exists, otherwise replace value by FLOAT_MIN in 1.iter
            float right_node = (rank * 2 + (1 << (d + 1)) - 1 >= total_angles and d == 0) ?
                               FLOAT_MIN : (*shared_angles)[rank * 2 + (1 << (d + 1)) - 1];
            // angles[i + 2(d+1) - 1] = max(angles[i + 2(d) - 1], angles[i + 2(d+1) - 1])
            (*shared_angles)[rank * 2 + (1 << (d + 1)) - 1] = std::max(left_node, right_node);
        }
        // Blocks until relevant processes in the iteration have computed the required results for the next iteration
        MPI_Barrier(MPI_COMM_WORLD);
    }

    // Master process sets the identity to the root to the root of the tree before the down-sweep phase
    if (rank == MASTER) {
        // Query the size and base pointer for a patch of a shared memory window with shared angles
        MPI_Win_shared_query(node_angles, MASTER, &master_window_size, &disp_unit, &(*shared_angles));
        // angles[n - 1] = FLOAT_MIN (I); set the neutral item to the root of the tree
        (*shared_angles)[NEXT_POWER_2(total_angles) - 1] = FLOAT_MIN;
    }

    // Down-Sweep phase of the max-prescan operation
    // for d from (lg n) - 1 downto 0
    for (int d = ceil(log2(total_angles)) - 1; d >= 0; d--) {
        // in parallel for i from 0 to n-1 by 2^(d+1)
        if (!(rank % (1 << (d)))) {
            // Query the size and base pointer for a patch of a shared memory window with shared angles
            MPI_Win_shared_query(node_angles, MASTER, &master_window_size, &disp_unit, &(*shared_angles));
            // Save the temporary - t = angles[i + 2^(d) - 1]
            float tmp = (*shared_angles)[rank * 2 + (1 << (d)) - 1];
            // Set the left child - angles[i + 2^(d) - 1] = angles[i + 2^(d+1) - 1]
            (*shared_angles)[rank * 2 + (1 << (d)) - 1] = (*shared_angles)[rank * 2 + (1 << (d + 1)) - 1];
            // Set the left child - angles[i + 2^(d+1) - 1] = max(t, angles[i + 2^(d+1) - 1])
            (*shared_angles)[rank * 2 + (1 << (d + 1)) - 1] =
                    std::max(tmp, (*shared_angles)[rank * 2 + (1 << (d + 1)) - 1]);
        }
        // Blocks until relevant processes in the iteration have computed the required results for the next iteration
        MPI_Barrier(MPI_COMM_WORLD);
    }
}

void preprocess_subsets(
    float **shared_angles, MPI_Win node_angles, MPI_Aint window_size, int disp_unit, int total_angles, int size, int rank
) {
    // Define the shared pointer to allocate shared window to store the maximums of the processors
    float *sub_max;
    // Create an area of memory for each processors to shared allocated memory (windows within shared array)
    MPI_Win node_sub_max;
    // Declare auxiliary window size to query at shared memory without the change relevant window size of processes
    MPI_Aint master_window_size;

    // Create an window for one-sided communication and shared memory access, and allocate memory at each process
    // Allocate shared window to store the maximums of the processors
    MPI_Win_allocate_shared(FLOAT_UNIT, FLOAT_UNIT, MPI_INFO_NULL, MPI_COMM_WORLD, &sub_max, &node_sub_max);
    // Query the size and base pointer for a patch of a shared memory window with shared angles
    MPI_Win_shared_query(node_angles, MASTER, &master_window_size, &disp_unit, &(*shared_angles));
    // Query the size and base pointer for a patch of a shared memory window with processors maximums
    MPI_Win_shared_query(node_sub_max, MASTER, &master_window_size, &disp_unit, &sub_max);
    // Compute the start index for each process within the shared memory with respect to the common begin
    int start_idx = (rank < total_angles % size) ? rank * window_size : rank * window_size + (total_angles % size);
    // Set the first from angle from the process window as the initial maximum
    // max[i] = angles[(n/p) * i]
    sub_max[rank] = (*shared_angles)[start_idx];
    // Each processor obtain the maximum angle from the its window
    // for j from 1 to n/p
    for (size_t i = start_idx; i < start_idx + window_size; i++) {
        // max[i] = max(max[i], angels[(n/p) * i + j])
        sub_max[rank] = std::max((*shared_angles)[i], sub_max[rank]);
    }

    // Blocks until all processes in the communicator have computed own n/p section maximum
    MPI_Barrier(MPI_COMM_WORLD);
    // Perform the max-prescan operation on the processor maximums - max-prescan(max)
    max_prescan(&sub_max, node_sub_max, window_size, disp_unit, size, size, rank);

    // Query the size and base pointer for a patch of a shared memory window with shared angles
    MPI_Win_shared_query(node_angles, MASTER, &master_window_size, &disp_unit, &(*shared_angles));
    // Query the size and base pointer for a patch of a shared memory window with result of the max-prescan operation
    MPI_Win_shared_query(node_sub_max, MASTER, &master_window_size, &disp_unit, &sub_max);
    // Temporary save of the first angle within process window
    float prev_point = (*shared_angles)[start_idx];
    // Set the maximum as offset for each process base on the results from the max-prescan operation
    (*shared_angles)[start_idx] = sub_max[rank];
    // Each processor prescan self n/p section
    for (size_t i = start_idx + 1; i < start_idx + window_size; i++) {
        // Save the current value of the maximum from the comparison before the rewrites temporary variable
        float max = std::max(prev_point, (*shared_angles)[i - 1]);
        // Save the current value of the angle before the rewrite it with the previous maximum value
        prev_point = (*shared_angles)[i];
        // Set the previous maximum value to the relevant index in the shared memory
        (*shared_angles)[i] = max;
    }
    // Blocks until all processes in the communicator have computed own n/p section within max-prescan operation
    MPI_Barrier(MPI_COMM_WORLD);
}

void compute_results(
    float **shared_angles, float **max_previous_angles, bool **result, MPI_Win node_angles, MPI_Win node_prev_angles,
    MPI_Win node_results, int disp_unit, int window_size, int total_angles, int size, int rank
) {
    // Declare auxiliary window size to query at shared memory without the change relevant window size of processes
    MPI_Aint master_window_size;
    // Compute the start index for each process within the shared memory with respect to the common begin
    int start_idx = (rank < total_angles % size) ? rank * window_size : rank * window_size + (total_angles % size);
    // Query the size and base pointer for a patch of a shared memory window with shared angles
    MPI_Win_shared_query(node_angles, MASTER, &master_window_size, &disp_unit, &(*shared_angles));
    // Query the size and base pointer for a patch of a shared memory window with maximum previous angles (max-prescan)
    MPI_Win_shared_query(node_prev_angles, MASTER, &master_window_size, &disp_unit, &(*max_previous_angles));
    // Query the size and base pointer for a patch of a shared memory window to store the final results
    MPI_Win_shared_query(node_results, MASTER, &master_window_size, &disp_unit, &(*result));
    // in parallel for each process window
    for (size_t i = start_idx; i < start_idx + window_size; i++) {
        // if (angles[i] > max-previous-angles[i]) result[i] = visible else not visible
        (*result)[i] = ((*shared_angles)[i] > (*max_previous_angles)[i]);
    }
    // Blocks until all processes in the communicator have computed own n/p section of final results
    MPI_Barrier(MPI_COMM_WORLD);
}

void write_out_result(
    bool **result, MPI_Win node_results, MPI_Aint window_size, int disp_unit, int total_altitudes
) {
    // Query the size and base pointer for a patch of a shared memory window with stored final results
    MPI_Win_shared_query(node_results, MASTER, &window_size, &disp_unit, &(*result));
    // Write out the place of the observation
    std::cout << "_,";
    // Write out the result for each origina altitude
    for (size_t i = 1; i < total_altitudes - 1; i++) {
        // Write out the result in the specified format - visible (v) or unvisible (u)
        std::cout << ((*result)[i] ? "v," : "u,");
    }
    // Write out the result for the last altitude (without the comma after it)
    std::cout << ((*result)[total_altitudes - 1] ? "v" : "u") << std::endl;
}

int main(int argc, char **argv) {
    // Create the variables to store information data within all processors
    int rank, size, total_altitudes, disp_unit;
    // Create an area of memory for each processors to shared allocated memory (windows within shared array)
    MPI_Win node_altitudes, node_angles, node_prev_angles, node_results;
    // Multiple of the displacement unit that was specified in the window definition for each process
    MPI_Aint window_size;
    // Define the vector to store the loaded altitudes from the input line of sight
    std::vector<int> altitudes;
    // Define the shared pointer to allocate shared window to store loaded altitudes between all processes
    int *shared_altitudes;
    // Define the shared pointers to allocate shared window to store computed angles between all processes
    float *shared_angles, *max_previous_angles;
    // Define the shared pointer to allocate shared window to store final results of the line-of-sight problem
    bool *result;

// Potential definition of the variable to measure the runtime of the algorithm for line-of-sight problem
#ifdef MEASURE_TIME
    double start = 0.0;
#endif

    // Initializes the MPI execution environment
    MPI_Init(&argc, &argv);
    // Returns the size of the group associated with a communicator
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    // Determines the rank of the calling process in the communicator
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // Master process loads the line-of-sight given on the input (first argument on the command line - argv[1])
    if (rank == MASTER) {
        // Load and parse the input line-of-sight to the vector of altitudes
        load_line_of_sight(argv[1], &altitudes);
        // Save the number of the loaded altitudes for subsequent sharing between all processes in this variable
        total_altitudes = altitudes.size();
    }
    // Broadcast a number of altitudes from the master process to all other process of the communicator MPI_COMM_WORLD
    MPI_Bcast(&total_altitudes, COUNT, MPI_INT, MASTER, MPI_COMM_WORLD);
    // Assign the relevant number of altitudes for each process and allocate the relevant shared memory
    share_points_to_process(&shared_altitudes, &node_altitudes, &window_size, total_altitudes, size, altitudes, rank);

    // Create an window for one-sided communication and shared memory access, and allocate memory at each process
    // Allocate shared window to store the computed angles.
    MPI_Win_allocate_shared(
            window_size * FLOAT_UNIT, FLOAT_UNIT, MPI_INFO_NULL, MPI_COMM_WORLD, &shared_angles, &node_angles
    );
    // Create an window for one-sided communication and shared memory access, and allocate memory at each process
    // Allocate shared window to store the maximum previous angles in the next phase of the algorithm
    MPI_Win_allocate_shared(
            window_size * FLOAT_UNIT, FLOAT_UNIT, MPI_INFO_NULL, MPI_COMM_WORLD, &max_previous_angles, &node_prev_angles
    );

// The starting point of measuring the runtime of the line-of-sight algorithm
    if (rank == MASTER) {
#ifdef MEASURE_TIME
        start = MPI_Wtime();
#endif
    }

    // Compute the angles by all processes and store the results in the given allocated shared memories
    compute_angles(
            &shared_altitudes, &shared_angles, &max_previous_angles, node_altitudes, node_angles,
            node_prev_angles, window_size, disp_unit, total_altitudes, size, rank
    );

    // Check whether is the required number of processes to perform only the max-prescan operation itself
    if (size < ceil(total_altitudes / 2.0)) {
        // When the number of processes is not satisfied, then first pre-processing the vector of angles
        preprocess_subsets(&max_previous_angles, node_prev_angles, window_size, disp_unit, total_altitudes, size, rank);
    } else {
        // When the number of processes is satisfied, then perform an only max-prescan operation itself
        max_prescan(&max_previous_angles, node_prev_angles, window_size, disp_unit, total_altitudes, size, rank);
    }

    // Create an window for one-sided communication and shared memory access, and allocate memory at each process
    // Allocate shared window to store the final results of the line-of-sight problem
    MPI_Win_allocate_shared(
            window_size * BOOL_UNIT, BOOL_UNIT, MPI_INFO_NULL, MPI_COMM_WORLD, &result, &node_results
    );
    // Compute the final results of the line-of-sight problem - subtraction of angle and maximum previous angle
    compute_results(
            &shared_angles, &max_previous_angles, &result, node_angles, node_prev_angles,
            node_results, disp_unit, window_size, total_altitudes, size, rank
    );

// The ending point of measuring the runtime of the line-of-sight algorithm
    if (rank == MASTER) {
#ifdef MEASURE_TIME
        double end = MPI_Wtime();
        // Write out the runtime of the algorithm to the standard output in the microseconds
        printf("%lu\n", (size_t)((end - start) * 1000000));
#endif
    }

    // Master process write the out the final results of the line-of sight problem
    if (rank == MASTER) {
        write_out_result(&result, node_results, window_size, disp_unit, total_altitudes);
    }

    // free shared allocated memories of each process
    MPI_Win_free(&node_altitudes);
    MPI_Win_free(&node_angles);
    MPI_Win_free(&node_prev_angles);
    MPI_Win_free(&node_results);

    // Terminates MPI execution environment
    MPI_Finalize();

    // End of the program
    return 0;
}
