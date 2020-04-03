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


void load_points(char *input_points, std::vector<int> *points, int size, int **send_counts, int **start_idx) {
    load_line_of_sight(input_points, &(*points));
    int points_per_processors = floor(((*points).size() - 1) / size);
    int modulo = ((*points).size() - 1) % size;
    for (size_t i = 0; i < size; i++) {
        (*send_counts)[i] = (i < modulo) ? points_per_processors + 1 : points_per_processors;
        (*start_idx)[i] = (i == 0) ? 0 : (*start_idx)[i - 1] + (*send_counts)[i - 1];
        MPI_Send(&((*send_counts)[i]), COUNT, MPI_INT, i, TAG, MPI_COMM_WORLD);
    }
}

void receive_points(
        int *points_number, int **proc_numbers, std::vector<int> points, int *send_counts,
        int *start_idx, MPI_Status *mpi_stat
) {
    MPI_Recv(&(*points_number), COUNT, MPI_INT, MASTER, TAG, MPI_COMM_WORLD, &(*mpi_stat));
    *proc_numbers = (int *) malloc(sizeof(int) * (*points_number));
    MPI_Scatterv(
            points.data(), send_counts, start_idx, MPI_INT, *proc_numbers, *points_number, MPI_INT, 0, MPI_COMM_WORLD
    );
}

void compute_angles(
        int* points, int count, int master_point, float** proc_angles, int rank, int total_points, int size
) {
    int start_idx = (rank < total_points % size ) ? rank * count + 1: rank * count + (total_points % size) + 1;
    for (size_t i = 0; i < count; i++) {
        (*proc_angles)[i] = atan((points[i] - master_point) / float(start_idx++));
    }
}

int main(int argc, char **argv) {
    // Create the variables to store info within all processors
    int rank, size, points_number, master_point, total_points;
    // Structure that represents the status of the received message
    MPI_Status mpi_stat;
    std::vector<int> points;
    int *proc_numbers;

    // Initializes the MPI execution environment
    MPI_Init(&argc, &argv);
    // Returns the size of the group associated with a communicator
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    // Determines the rank of the calling process in the communicator
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int *send_counts = (int *) malloc(sizeof(int) * size);
    int *start_idx = (int *) malloc(sizeof(int) * size);

    if (rank == MASTER) {
        load_points(argv[1], &points, size, &send_counts, &start_idx);
        master_point = points[0];
        points.erase(points.begin());
        total_points = points.size();
    }
    MPI_Bcast(&master_point, COUNT, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&total_points, COUNT, MPI_INT, 0, MPI_COMM_WORLD);
    receive_points(&points_number, &proc_numbers, points, send_counts, start_idx, &mpi_stat);

    float* proc_angles = (float*) malloc(sizeof(float) * points_number);
    compute_angles(proc_numbers, points_number, master_point, &proc_angles, rank, total_points, size);
    float* angles = (float*) malloc(sizeof(float) * total_points);
    MPI_Gatherv(proc_angles, points_number, MPI_FLOAT, angles, send_counts, start_idx, MPI_FLOAT, MASTER, MPI_COMM_WORLD);

    // Terminates MPI execution environment
    MPI_Finalize();

    // End of the program
    return 0;
}
