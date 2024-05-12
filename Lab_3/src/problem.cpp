#include <eigen3/Eigen/Dense>
#include <fstream>
#include <iostream>

#include "../include/problem.h"
#ifndef progress_bar
#define progress_bar
    #include "../include/progress_bar.h"
#endif

Problem::Problem(double alpha, double L, double T, int n, int k, double sigma,
                 std::function<double(double)> u_x,
                 std::function<double(double)> u_t,
                 std::function<double(double, double)> exact_solution):
                 alpha{alpha}, L{L}, T{T}, sigma{sigma},
                 u_x{u_x}, u_t{u_t}, exact_solution{exact_solution} {
    change_n_k(n, k);
}

void Problem::change_n_k(int n_new, int k_new) {
    n = n_new;
    k = k_new;
    h = L / n;
    tau = T / k;

    numerical_solution = Eigen::MatrixXd(n + 1, k + 1);
}

double Problem::get_error() {
    double error = 0;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < k; j++) {
            error = std::max(error, std::abs(numerical_solution(i, j) - exact_solution(h * i, tau * j)));
        }
    }

    return error;
}

void Problem::save_solution(int frames, int points) {
    if ((n % points != 0) or (k % frames != 0)) {
        std::cout << "Invalid number of frames/points\n";
    }

    std::ofstream output;
    output.open("../out/animation_x.csv");
    for (int i = 0; i < n; i += n / points) {
        output << i << ',';
    }
    output << n << '\n';

    for (int i = 0; i < n; i += n / points) {
        output << (L * i) / n << ',';
    }
    output << L << '\n';
    output.close();

    output.open("../out/animation_y.csv");
    for (int i = 0; i < n; i += n / points) {
        output << i << ',';
    }
    output << n << '\n';

    for (int j = 0; j <= k; j += k / frames) {
        for (int i = 0; i < n; i += n / points) {
            output << numerical_solution(i, j) << ',';
        }
        output << numerical_solution(n, j) << '\n';
    }
    output.close();

    output.open("../out/animation_y_exact.csv");
    for (int i = 0; i < n; i += n / points) {
        output << i << ',';
    }
    output << n << '\n';

    for (int j = 0; j <= k; j += k / frames) {
        for (int i = 0; i < n; i += n / points) {
            output << exact_solution((L * i) / n, (T * j) / k) << ',';
        }
        output << exact_solution(L, (T * j) / k) << '\n';
    }
    output.close();
}