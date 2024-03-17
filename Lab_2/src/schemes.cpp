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
        MPI_Send(&y_cur[length - 1], 1, MPI_DOUBLE, rank + 1, j, MPI_COMM_WORLD);
    }

    for(int i = 1; i < length - 1; i++) {
        y_cur[i] = y_prev[i] - c * (y_prev[i] - y_prev[i - 1]);
    }

    double y_left;
    if (rank == 0) {
        y_left = problem.y_0(time_coef * (j - 1));
    } else if (j == 1){
        y_left = problem.f_0(x[0] - h);
    } else {
        MPI_Status status;
        MPI_Recv(&y_left, 1, MPI_DOUBLE, rank - 1, j - 1, MPI_COMM_WORLD, &status);
    }

    y_cur[0] = y_prev[0] - c * (y_prev[0] - y_left);
}

//RightAngleSolver

RightAngleSolver::RightAngleSolver(int size, int rank, int N, int K, ConvectionDiffusionProblem problem) : 
                          Solver(size, rank, N, K, problem) {
    scheme = "right_angle";
}

void RightAngleSolver::make_step(int j) {
    std::swap(y_cur, y_prev);

    y_cur[0] = y_prev[0] - c * (y_prev[1] - y_prev[0]);

    if (rank != 0 && j < K) {
        MPI_Send(&y_cur[0], 1, MPI_DOUBLE, rank - 1, j, MPI_COMM_WORLD);
    }

    for(int i = 1; i < length - 1; i++) {
        y_cur[i] = y_prev[i] - c * (y_prev[i + 1] - y_prev[i]);
    }

    double y_right;
    if (rank == size - 1) {
        y_right = 0;
    } else if (j == 1){
        y_right = problem.f_0(x[length - 1] + h);
    } else {
        MPI_Status status;
        MPI_Recv(&y_right, 1, MPI_DOUBLE, rank + 1, j - 1, MPI_COMM_WORLD, &status);
    }

    y_cur[length - 1] = y_prev[length - 1] - c * (y_right - y_prev[length - 1]);
}

//ImplicitAngleSolver

ImplicitAngleSolver::ImplicitAngleSolver(int size, int rank, int N, int K, ConvectionDiffusionProblem problem) : 
                          Solver(size, rank, N, K, problem) {
    scheme = "implicit_angle";
}

void ImplicitAngleSolver::make_step(int j) {
    std::swap(y_cur, y_prev);

    double y_left;
    if (rank == 0) {
        y_left = problem.y_0(time_coef * (j - 1));
    } else if (j == 1) {
        y_left = problem.f_0(x[0] - h);
    } else {
        MPI_Status status;
        MPI_Recv(&y_left, 1, MPI_DOUBLE, rank - 1, j, MPI_COMM_WORLD, &status);
    }

    y_cur[0] = (c * y_left + y_prev[0]) / (c + 1);

    for(int i = 1; i < length; i++) {
        y_cur[i] = (c * y_cur[i - 1] + y_prev[i]) / (c + 1);
    }

    if (rank != size - 1) {
        MPI_Send(&y_cur[length - 1], 1, MPI_DOUBLE, rank + 1, j, MPI_COMM_WORLD);
    }
}