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
        MPI_Recv(&y_left, 1, MPI_DOUBLE, rank - 1, j - 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
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
        MPI_Recv(&y_right, 1, MPI_DOUBLE, rank + 1, j - 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    y_cur[length - 1] = y_prev[length - 1] - c * (y_right - y_prev[length - 1]);
}

//ImplicitAngleSolver

ImplicitAngleSolver::ImplicitAngleSolver(int size, int rank, int N, int K, ConvectionDiffusionProblem problem) : 
                          Solver(size, rank, N, K, problem) {
    scheme = "implicit_angle";
    d = 1 / (c + 1);
}

void ImplicitAngleSolver::make_step(int j) {
    std::swap(y_cur, y_prev);

    double y_left;
    if (rank == 0) {
        y_left = problem.y_0(time_coef * (j - 1));
    } else {
        MPI_Recv(&y_left, 1, MPI_DOUBLE, rank - 1, j, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    y_cur[0] = d * (c * y_left + y_prev[0]);

    for(int i = 1; i < length; i++) {
        y_cur[i] = d * (c * y_cur[i - 1] + y_prev[i]);
    }

    if (rank < size - 1) {
        MPI_Send(&y_cur[length - 1], 1, MPI_DOUBLE, rank + 1, j, MPI_COMM_WORLD);
    }
}

//FourPointSolver

FourPointSolver::FourPointSolver(int size, int rank, int N, int K, ConvectionDiffusionProblem problem) : 
                          Solver(size, rank, N, K, problem) {
    scheme = "four_point";
    d = 1 + 2 / c;
}

void FourPointSolver::make_step(int j) {
    std::swap(y_cur, y_prev);
    std::swap(y_cur_left, y_prev_left);

    if (rank == 0) {
        y_cur[0] = problem.y_0(time_coef * j);
        y_cur[1] = y_prev[1] - c * (y_prev[1] - y_prev[0]);
    } else {
        if (j == 1) {
            y_prev_left[0] = problem.f_0(x[0] - 2 * h);
            y_prev_left[1] = problem.f_0(x[0] - h);
        }
        MPI_Recv(&y_cur_left, 2, MPI_DOUBLE, rank - 1, j, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        y_cur[0] = d * y_cur_left[1] - d * y_prev_left[1] + y_prev_left[0];
        y_cur[1] = d * y_cur[0] - d * y_prev[0] + y_prev_left[1];
    }

    for(int i = 2; i < length; i++) {
        y_cur[i] = d * y_cur[i - 1] - d * y_prev[i - 1] + y_prev[i - 2];
    }

    if (rank != size - 1) {
        MPI_Send(&y_cur[length - 2], 2, MPI_DOUBLE, rank + 1, j, MPI_COMM_WORLD);
    }
}

//LaxSolver

LaxSolver::LaxSolver(int size, int rank, int N, int K, ConvectionDiffusionProblem problem) : 
                          Solver(size, rank, N, K, problem) {
    scheme = "lax";
    d_1 = (1 - c) / 2;
    d_2 = (1 + c) / 2;
}

void LaxSolver::make_step(int j) {
    std::swap(y_cur, y_prev);

    double y_left; double y_right;

    if (rank == 0) {
        MPI_Sendrecv(&y_prev[length - 1], 1, MPI_DOUBLE, 1, j - 1, 
                     &y_right, 1, MPI_DOUBLE, 1, K + j - 1,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    } else if (rank < size - 1) {
        MPI_Sendrecv(&y_prev[length - 1], 1, MPI_DOUBLE, rank + 1, j - 1, 
                     &y_left, 1, MPI_DOUBLE, rank - 1, j - 1,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Sendrecv(&y_prev[0], 1, MPI_DOUBLE, rank - 1, K + j - 1, 
                     &y_right, 1, MPI_DOUBLE, rank + 1, K + j - 1,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    } else {
        MPI_Sendrecv(&y_prev[0], 1, MPI_DOUBLE, rank - 1, K + j - 1, 
                     &y_left, 1, MPI_DOUBLE, rank - 1, j - 1,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    if (rank == 0) {
        y_cur[0] = problem.y_0(time_coef * j);
    } else {
        y_cur[0] = d_1 * y_prev[1] + d_2 * y_left;
    }

    for(int i = 1; i < length - 1; i++) {
        y_cur[i] = d_1 * y_prev[i + 1] + d_2 * y_prev[i - 1];
    }

    if (rank == size - 1) {
        y_cur[length - 1] = y_prev[length - 1] - c * (y_prev[length - 1] - y_prev[length - 2]);
    } else {
        y_cur[length - 1] = d_1 * y_right + d_2 * y_prev[length - 2];
    }
}

//LaxWendroffSolver

LaxWendroffSolver::LaxWendroffSolver(int size, int rank, int N, int K, ConvectionDiffusionProblem problem) : 
                          Solver(size, rank, N, K, problem) {
    scheme = "lax_wendroff";
    d_1 = c * (c + 1) / 2;
    d_2 = 1 - c * c;
    d_3 = c * (c - 1) / 2;
}

void LaxWendroffSolver::make_step(int j) {
    std::swap(y_cur, y_prev);
    
    double y_left; double y_right;

    if (rank == 0) {
        MPI_Sendrecv(&y_prev[length - 1], 1, MPI_DOUBLE, 1, j - 1, 
                     &y_right, 1, MPI_DOUBLE, 1, K + j - 1,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    } else if (rank < size - 1) {
        MPI_Sendrecv(&y_prev[length - 1], 1, MPI_DOUBLE, rank + 1, j - 1, 
                     &y_left, 1, MPI_DOUBLE, rank - 1, j - 1,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Sendrecv(&y_prev[0], 1, MPI_DOUBLE, rank - 1, K + j - 1, 
                     &y_right, 1, MPI_DOUBLE, rank + 1, K + j - 1,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    } else {
        MPI_Sendrecv(&y_prev[0], 1, MPI_DOUBLE, rank - 1, K + j - 1, 
                     &y_left, 1, MPI_DOUBLE, rank - 1, j - 1,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    if (rank == 0) {
        y_cur[0] = problem.y_0(time_coef * j);
    } else {
        y_cur[0] = d_1 * y_left + d_2 * y_prev[0] + d_3 * y_prev[1];
    }

    for(int i = 1; i < length - 1; i++) {
        y_cur[i] = d_1 * y_prev[i - 1] + d_2 * y_prev[i] + d_3 * y_prev[i + 1];
    }

    if (rank == size - 1) {
        y_cur[length - 1] = y_prev[length - 1] - c * (y_prev[length - 1] - y_prev[length - 2]);
    } else {
        y_cur[length - 1] = d_1 * y_prev[length - 2] + d_2 * y_prev[length - 1] + d_3 * y_right;
    }
}

//Cross solver

CrossSolver::CrossSolver(int size, int rank, int N, int K, ConvectionDiffusionProblem problem) : 
                          Solver(size, rank, N, K, problem) {
    scheme = "cross";
    y_prev_2 = std::vector(length, 0.0);
}

void CrossSolver::make_step(int j) {
    std::swap(y_prev, y_prev_2);
    std::swap(y_cur, y_prev);

    double y_left; double y_right;

    if (j == 1) {
        if (rank == 0) {
            y_cur[0] = problem.y_0(time_coef * j);
        } else {
            y_cur[0] = y_prev[0] - c * (y_prev[0] - problem.f_0(x[0] - h));
        }
        for (int i = 1; i < length; i++) {
            y_cur[i] = y_prev[i] - c * (y_prev[i] - y_prev[i - 1]);
        }
    } else {

        if (rank == 0) {
            MPI_Sendrecv(&y_prev[length - 1], 1, MPI_DOUBLE, 1, j - 1, 
                         &y_right, 1, MPI_DOUBLE, 1, K + j - 1,
                         MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            
        } else if (rank < size - 1) {
            MPI_Sendrecv(&y_prev[length - 1], 1, MPI_DOUBLE, rank + 1, j - 1, 
                         &y_left, 1, MPI_DOUBLE, rank - 1, j - 1,
                         MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Sendrecv(&y_prev[0], 1, MPI_DOUBLE, rank - 1, K + j - 1, 
                         &y_right, 1, MPI_DOUBLE, rank + 1, K + j - 1,
                         MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        } else {
            MPI_Sendrecv(&y_prev[0], 1, MPI_DOUBLE, rank - 1, K + j - 1, 
                         &y_left, 1, MPI_DOUBLE, rank - 1, j - 1,
                         MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

        if (rank == 0) {
            y_cur[0] = problem.y_0(time_coef * j);
        } else {
            y_cur[0] = y_prev_2[0] - c * (y_prev[1] - y_left);
        }

        for(int i = 1; i < length - 1; i++) {
            y_cur[i] = y_prev_2[i] - c * (y_prev[i + 1] - y_prev[i - 1]);
        }

        if (rank == size - 1) {
            y_cur[length - 1] = y_prev[length - 1] - c * (y_prev[length - 1] - y_prev[length - 2]);
        } else {
            y_cur[length - 1] = y_prev_2[length - 1] - c * (y_right - y_prev[length - 2]);
        }
    }
}