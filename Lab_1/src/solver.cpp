#include <eigen3/Eigen/Dense>
#include <iostream>

#include "solver.h"


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
        x[i] = (i / N) * h;
    }
    y = new double[N + 1];

    Psi = [p, q, f](double x, Eigen::Vector2d u) { return Eigen::Vector2d(u(1), f(x) - p(x) * u(1) - q(x) * u(0)); };
}

double** BoundaryValueProblem::GetResults(int n) {
    auto res = new double*[2];
    res[0] = new double[n + 1];
    res[1] = new double[n + 1];
    for(int i = 0; i <= n; i++) {
        res[0][i] = x[(i * N) / n];
        res[1][i] = y[(i * N) / n];
    }

    return res;
}

Eigen::Vector2d BoundaryValueProblem::Shoot(Eigen::Vector2d u) {
    y[0] = u(0);
    for(int i = 0; i < N; i++) {
        Eigen::Vector2d k1, k2, k3, k4;
        k1 = Psi(x[i], u);
        k2 = Psi(x[i] + 0.5 * h, u + 0.5 * k1);
        k3 = Psi(x[i] +  0.5 * h, u + 0.5 * k2);
        k4 = Psi(x[i] + h, u + k3);

        u += (h / 6) * ( k1 + 2 * k2 + 2 * k3 + k4);
        y[i + 1] = u(0);
    }

    return u;
}

double FindZero(std::function<double(double)> Func, double l, double r, double f_l, double f_r, double eps) {
    if(f_l >= f_r) {
        std::swap(l, r);
        std::swap(f_l, f_r);
    }
    if(f_l > 0 || f_r < 0) {
        std::cout << "Invalid range!" << '\n';
        return 1;
    }

    while(std::abs(r - l) > eps) {
        double m = (l + r) / 2;
        double f_m = Func(m);
        if(f_m > 0) {
            r = m; f_r = f_m;
        } 
        else {
            l = m; f_l = f_m;
        }
    }
    return (l + r) / 2;
}

Eigen::Vector2d BoundaryValueProblem::FindInitialVals(double eps, int max_iter) {
    std::function<Eigen::Vector2d(double)> Ivals;
    if(b_1 == 0) {
        Ivals = [this](double x) { return Eigen::Vector2d(c_1, x);};
    }
    else {
        Ivals = [this](double x) { return Eigen::Vector2d(x, (c_1 - a_1 * x) / b_1); };
    }

    auto Funk = [this, Ivals](double x){ Eigen::Vector2d fvals = Shoot(Ivals(x));
                                         return a_2 * fvals(0) + b_2 * fvals(1) - c_2; };
    
    double f_l = Funk(0);
    double f_r = Funk(1);
    int i = -1;
    while(std::abs(i) < std::pow(2, max_iter - 1) && (f_l > 0) == (f_r > 0)) {
        i = - (i + (i > 0));
        f_r = Funk(i);
    }
    return Ivals(FindZero(Funk, 0, i, f_l, f_r, eps));
}
