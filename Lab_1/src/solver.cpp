#include <eigen3/Eigen/Dense>

#include "../include/solver.h"


BoundaryValueProblem::BoundaryValueProblem(std::function<double(double)> p,
                                           std::function<double(double)> q, 
                                           std::function<double(double)> f, 
                                           double x_1, double x_2, int N,
                                           double a_1, double b_1, double c_1,
                                           double a_2, double b_2, double c_2):
                                           p{p}, q{q}, f{f},
                                           x_1{x_1}, x_2{x_2}, N{N},
                                           a_1{a_1}, b_1{b_1}, c_1{c_1},
                                           a_2{a_2}, b_2{b_2}, c_2{c_2} {
    h = (x_2 - x_1) / N;
    x = new double[N + 1];
    for(int i = 0; i <= N; i++) {
        x[i] = x_1 + (x_2 - x_1) * i / N;
    }
    y = new double[N + 1];

    Psi = [p, q, f](double x, Eigen::Vector2d u) { return Eigen::Vector2d(u(1), f(x) - p(x) * u(1) - q(x) * u(0)); };
}

BoundaryValueProblem::BoundaryValueProblem(std::function<double(double)> p,
                                           std::function<double(double)> q, 
                                           std::function<double(double)> f, 
                                           double x_1, double x_2, int N,
                                           double a_1, double b_1, double c_1,
                                           double a_2, double b_2, double c_2,
                                           std::function<double(double)> exact_solution):
                                           BoundaryValueProblem(p, q, f, x_1, x_2, N, a_1, b_1, c_1, a_2, b_2, c_2) {
this -> exact_solution = exact_solution;
}

BoundaryValueProblem::~BoundaryValueProblem() {
    delete [] x;
    delete [] y;
}

std::vector<std::vector<double>> BoundaryValueProblem::GetResults(int n) {
    auto tmp = std::vector(n + 1, 0.0);
    auto res = std::vector(2, tmp);

    for(int i = 0; i <= n; i++) {
        res[0][i] = x[(i * N) / n];
        res[1][i] = y[(i * N) / n];
    }

    return res;
}

double BoundaryValueProblem::GetError() {
    double error = 0;
    for(int i = 0; i <= N; i++) {
        double error_cur = std::abs(exact_solution(x[i]) - y[i]);
        if(error_cur > error) {
            error = error_cur;
        }
    }

    return error;
}

void BoundaryValueProblem::CnangeN(int N) {
    this -> N = N;
    h = (x_2 - x_1) / N;
    delete[] x;
    delete[] y;
    x = new double[N + 1];
    for(int i = 0; i <= N; i++) {
        x[i] = x_1 + (x_2 - x_1) * i / N;
    }
    y = new double[N + 1];
}