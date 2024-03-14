#include <vector>
#include <math.h>
#include <fstream>
#include <format>
#include <mpi.h>

#include "../include/solver.h"

Solver::Solver(int size, int rank, int N, int K, ConvectionDiffusionProblem problem):
               size{size}, rank{rank}, N{N}, K{K}, problem{problem} {

    left = (rank * (N + 1)) / size;
    right = ((rank + 1) * (N + 1)) / size;
    length = right - left;

    h = problem.L / N;
    tau = problem.T / K;
    c = problem.a * tau / h;
    time_coef = 1 / (K * problem.T);

    y_cur = std::vector(length, 0.0);
    y_prev = std::vector(length, 0.0);

    x = std::vector(length, 0.0);
    for (int i = 0; i < length; i++) {
        x[i] = (problem.L * (i + left)) / N;
    }
}

double Solver::get_error(double t) {
    double max_error = 0;
    double error = 0;

    for (int i = 0; i < length; i++) {
        error = std::abs(y_cur[i] - problem.y_exact(x[i], t));
        if(error > max_error) {
            max_error = error;
        }
    }

    return max_error;
}

double Solver::solve(int frames) {

    double max_error = 0;
    double error = 0;

    for (int i = 0; i < length; i++) {
        y_cur[i] = problem.f_0(x[i]);
    }

    //initializing writing to a file
    std::ofstream output;

    if (frames != 0) {
        output.open("../out/animations/" + scheme + "/x_" + std::to_string(rank) + ".csv");
        for (int i = 0; i < length - 1; i++) {
            output << i << ',';
        }
        output << length - 1 << '\n';
        for (int i = 0; i < length - 1; i++) {
            output << x[i] << ',';
        }
        output << x[length - 1];
        output.close();

        output.open("../out/animations/" + scheme + "/y_" + std::to_string(rank) + ".csv");
        for (int i = 0; i < length - 1; i++) {
                output << i << ',';
            }
            output << length - 1 << '\n';
        for (int i = 0; i < length - 1; i++) {
                output << y_cur[i] << ',';
            }
            output << y_cur[length - 1] << '\n';
    }
    //calculating
    for(int j = 1; j <= K; j++) {
        make_step(j);

        if (frames != 0 && j % (K / frames) == 0) {
            for(int i = 0; i < length - 1; i++) {
                output << y_cur[i] << ',';
            }
            output << y_cur[length - 1] << '\n';
        }

        error = get_error(j / K * problem.T);
        if (error > max_error) {
            max_error = error;
        }
    }

    //closing a file
    if (frames != 0) {
        output.close();
    }
 
    return max_error;
}

void Solver::make_step(int j) {
    //do nothing
}

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