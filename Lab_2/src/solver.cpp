#include <vector>
#include <math.h>
#include <fstream>
#include <format>
#include <string>
#include <mpi.h>

#include "../include/solver.h"

Solver::Solver(int size, int rank, int N, int K, ConvectionDiffusionProblem problem):
               size{size}, rank{rank}, N{N}, K{K}, problem{problem} {

    change_N_K(N, K);
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

double Solver::get_h() {
    return h;
}

double Solver::get_tau() {
    return tau;
}

std::string Solver::get_scheme() {
    return scheme;
}

void Solver::change_N(int N_new) {
    N = N_new;
    left = (rank * (N + 1)) / size;
    right = ((rank + 1) * (N + 1)) / size;
    length = right - left;

    h = problem.L / N;
    c = problem.a * tau / h;

    y_cur = std::vector(length, 0.0);
    y_prev = std::vector(length, 0.0);

    x = std::vector(length, 0.0);
    for (int i = 0; i < length; i++) {
        x[i] = (problem.L * (i + left)) / N;
    }
}

void Solver::change_K(int K_new) {
    K = K_new;
    tau = problem.T / K;
    c = problem.a * tau / h;
    time_coef = 1 / (K * problem.T);
}

void Solver::change_N_K(int N_new, int K_new) {
    change_N(N_new);
    change_K(K_new);
}

double Solver::solve(int frames) {

    double max_error = 0;
    double error = 0;

    for (int i = 0; i < length; i++) {
        y_cur[i] = problem.f_0(x[i]);
    }

    //initializing writing to a file
    std::ofstream output;
    std::ofstream output_exact;

    if (frames != 0) {
        output.open("../out/animations/" + scheme + "/" + problem.name + "/x_" + std::to_string(rank) + ".csv");
        for (int i = 0; i < length - 1; i++) {
            output << i << ',';
        }
        output << length - 1 << '\n';
        for (int i = 0; i < length - 1; i++) {
            output << x[i] << ',';
        }
        output << x[length - 1];
        output.close();

        output.open("../out/animations/" + scheme + "/" + problem.name + "/y_" + std::to_string(rank) + ".csv");
        for (int i = 0; i < length - 1; i++) {
                output << i << ',';
            }
            output << length - 1 << '\n';
        for (int i = 0; i < length - 1; i++) {
                output << y_cur[i] << ',';
            }
            output << y_cur[length - 1] << '\n';

        output_exact.open("../out/animations/" + scheme + "/" + problem.name + "/y_exact_" + std::to_string(rank) + ".csv");
        for (int i = 0; i < length - 1; i++) {
                output_exact << i << ',';
            }
            output_exact << length - 1 << '\n';
        for (int i = 0; i < length - 1; i++) {
                output_exact << problem.y_exact(x[i], 0) << ',';
            }
            output_exact << problem.y_exact(x[length - 1], 0) << '\n';
    }
    //calculating
    for(int j = 1; j <= K; j++) {
        make_step(j);

        if (frames != 0 && j % (K / frames) == 0) {
            for(int i = 0; i < length - 1; i++) {
                output << y_cur[i] << ',';
            }
            output << y_cur[length - 1] << '\n';

            for (int i = 0; i < length - 1; i++) {
                output_exact << problem.y_exact(x[i], time_coef * j) << ',';
            }
            output_exact << problem.y_exact(x[length - 1], time_coef * j) << '\n';
        }

        error = get_error(time_coef * j);
        if (error > max_error) {
            max_error = error;
        }
    }

    //closing a file
    if (frames != 0) {
        output.close();
        output_exact.close();
    }
 
    return max_error;
}

void Solver::make_step(int j) {
    //do nothing
}