#include <mpi.h>

#include "../include/schemes.h"

//LeftAngleSolver

LeftAngleSolver::LeftAngleSolver(int size, int rank, int N, int K, ConvectionDiffusionProblem problem) : 
                          Solver(size, rank, N, K, problem) {
    scheme = "left_angle";
}

void LeftAngleSolver::make_step(int j) {
    std::swap(y_cur, y_prev);

    y_cur[length - 1] = y_prev[length - 1] - c * (y_prev[length - 1] - y_prev[length - 2]);

    if (rank != size - 1 && j < K) {
        MPI_Send(&y_cur[length - 1], 1, MPI_DOUBLE, rank + 1, j + 1, MPI_COMM_WORLD);
    }

    for(int i = 1; i < length - 1; i++) {
        y_cur[i] = y_prev[i] - c * (y_prev[i] - y_prev[i - 1]);
    }

    double y_left;
    if (rank == 0) {
        y_left = problem.y_0(time_coef * j);
    } else if (j == 1){
        y_left = problem.f_0(x[0] - h);
    } else {
        MPI_Status status;
        MPI_Recv(&y_left, 1, MPI_DOUBLE, rank - 1, j, MPI_COMM_WORLD, &status);
    }

    y_cur[0] = y_prev[0] - c * (y_prev[0] - y_left);
}

//RightAngleSolver

