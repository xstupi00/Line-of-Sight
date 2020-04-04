#include "vid.h"


void load_line_of_sight(char *input_line_of_sight, std::vector<int> *target_line_of_sight) {
    char *point = strtok(input_line_of_sight, ",");
    (*target_line_of_sight).push_back(atoi(point));
    while (point) {
        point = strtok(NULL, ",");
        if (point) {
            (*target_line_of_sight).push_back(atoi(point));
        }
    }
}

void share_points_to_process(
        int **window_data, MPI_Win *node_window, MPI_Aint *window_size,
        int total_points, int size, std::vector<int> points, int rank
) {
    int points_per_processors = floor(total_points / size);
    int modulo = total_points % size;
    *window_size = (rank < modulo) ? points_per_processors + 1 : points_per_processors;

    MPI_Win_allocate_shared(
            (*window_size) * INT_UNIT, INT_UNIT, MPI_INFO_NULL, MPI_COMM_WORLD, &(*window_data), &(*node_window)
    );
    if (rank == MASTER) {
        std::copy(points.begin(), points.end(), *window_data);
    }
}

void compute_angles(
        int **shared_points, float **shared_angles, float **max_previous_angles, MPI_Win node_points,
        MPI_Win node_angles, MPI_Win node_prev_angles, MPI_Aint points_window_size,
        int disp_unit, int total_points, int size, int rank
) {
    MPI_Aint master_window_size;
    int start_idx = (rank < total_points % size) ?
                    rank * points_window_size : rank * points_window_size + (total_points % size);
    MPI_Win_shared_query(node_points, MASTER, &master_window_size, &disp_unit, &(*shared_points));
    MPI_Win_shared_query(node_angles, MASTER, &master_window_size, &disp_unit, &(*shared_angles));
    MPI_Win_shared_query(node_prev_angles, MASTER, &master_window_size, &disp_unit, &(*max_previous_angles));
    for (size_t i = start_idx; i < start_idx + points_window_size; i++) {
        (*shared_angles)[i] = (*max_previous_angles)[i] =
                (i) ? atan(((*shared_points)[i] - (*shared_points)[0]) / float(i)) : FLOAT_MIN;
    }
}


void max_prescan(
        float **shared_angles, MPI_Win node_angles, MPI_Aint points_window_size,
        int disp_unit, int total_points, int size, int rank
) {
    // up_sweep
    MPI_Aint master_window_size;
    for (size_t d = 0; d < ceil(log2(total_points)); d++) {
        if (!(rank % (1 << (d)))) {
            MPI_Win_shared_query(node_angles, MASTER, &master_window_size, &disp_unit, &(*shared_angles));
            float left_node = (rank * 2 + (1 << (d)) - 1 >= total_points and d == 0) ?
                              FLOAT_MIN : (*shared_angles)[rank * 2 + (1 << (d)) - 1];
            float right_node = (rank * 2 + (1 << (d + 1)) - 1 >= total_points and d == 0) ?
                               FLOAT_MIN : (*shared_angles)[rank * 2 + (1 << (d + 1)) - 1];
            (*shared_angles)[rank * 2 + (1 << (d + 1)) - 1] = std::max(left_node, right_node);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }

    if (rank == MASTER) {
        MPI_Win_shared_query(node_angles, MASTER, &master_window_size, &disp_unit, &(*shared_angles));
        (*shared_angles)[NEXT_POWER_2(total_points) - 1] = FLOAT_MIN;
    }
    MPI_Barrier(MPI_COMM_WORLD);

    // down sweep
    for (int d = ceil(log2(total_points)) - 1; d >= 0; d--) {
        if (!(rank % (1 << (d)))) {
            MPI_Win_shared_query(node_angles, MASTER, &master_window_size, &disp_unit, &(*shared_angles));
            float tmp = (*shared_angles)[rank * 2 + (1 << (d)) - 1];
            (*shared_angles)[rank * 2 + (1 << (d)) - 1] = (*shared_angles)[rank * 2 + (1 << (d + 1)) - 1];
            (*shared_angles)[rank * 2 + (1 << (d + 1)) - 1] =
                    std::max(tmp, (*shared_angles)[rank * 2 + (1 << (d + 1)) - 1]);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
}

void compute_results(
        float **shared_angles, float **max_previous_angles, bool **result,
        MPI_Win node_angles, MPI_Win node_prev_angles, MPI_Win node_results,
        int disp_unit, int window_size, int total_points, int size, int rank
) {
    MPI_Aint master_window_size;
    int start_idx = (rank < total_points % size) ? rank * window_size : rank * window_size + (total_points % size);
    MPI_Win_shared_query(node_angles, MASTER, &master_window_size, &disp_unit, &(*shared_angles));
    MPI_Win_shared_query(node_prev_angles, MASTER, &master_window_size, &disp_unit, &(*max_previous_angles));
    MPI_Win_shared_query(node_results, MASTER, &master_window_size, &disp_unit, &(*result));
    for (size_t i = start_idx; i < start_idx + window_size; i++) {
        (*result)[i] = ((*shared_angles)[i] > (*max_previous_angles)[i]);
    }
}

void write_out_result(
        bool **result, MPI_Win node_results, MPI_Aint window_size, int disp_unit, int total_points
) {
    MPI_Win_shared_query(node_results, MASTER, &window_size, &disp_unit, &(*result));
    std::cout << "_,";
    for (size_t i = 1; i < total_points - 1; i++) {
        std::cout << ((*result)[i] ? "v," : "u,");
    }
    std::cout << ((*result)[total_points - 1] ? "v" : "u") << std::endl;
}

int main(int argc, char **argv) {
    // Create the variables to store info within all processors
    int rank, size, total_points, disp_unit;
    // Structure that represents the status of the received message
    MPI_Status mpi_stat;
    MPI_Win node_points, node_angles, node_prev_angles, node_results;
    MPI_Aint window_size;
    std::vector<int> points;
    int *shared_points;
    float *shared_angles, *max_previous_angles;
    bool *result;

    // Initializes the MPI execution environment
    MPI_Init(&argc, &argv);
    // Returns the size of the group associated with a communicator
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    // Determines the rank of the calling process in the communicator
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank == MASTER) {
        load_line_of_sight(argv[1], &points);
        total_points = points.size();
    }
    MPI_Bcast(&total_points, COUNT, MPI_INT, 0, MPI_COMM_WORLD);
    share_points_to_process(&shared_points, &node_points, &window_size, total_points, size, points, rank);

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Win_allocate_shared(
            window_size * FLOAT_UNIT, FLOAT_UNIT, MPI_INFO_NULL, MPI_COMM_WORLD, &shared_angles, &node_angles
    );
    MPI_Win_allocate_shared(
            window_size * FLOAT_UNIT, FLOAT_UNIT, MPI_INFO_NULL, MPI_COMM_WORLD, &max_previous_angles, &node_prev_angles
    );
    compute_angles(
            &shared_points, &shared_angles, &max_previous_angles, node_points, node_angles,
            node_prev_angles, window_size, disp_unit, total_points, size, rank
    );

    MPI_Barrier(MPI_COMM_WORLD);
    max_prescan(&max_previous_angles, node_prev_angles, window_size, disp_unit, total_points, size, rank);

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Win_allocate_shared(
            window_size * BOOL_UNIT, BOOL_UNIT, MPI_INFO_NULL, MPI_COMM_WORLD, &result, &node_results
    );
    compute_results(
            &shared_angles, &max_previous_angles, &result, node_angles, node_prev_angles,
            node_results, disp_unit, window_size, total_points, size, rank
    );

    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == MASTER) {
        write_out_result(&result, node_results, window_size, disp_unit, total_points);
    }

    // Terminates MPI execution environment
    MPI_Finalize();

    // End of the program
    return 0;
}
