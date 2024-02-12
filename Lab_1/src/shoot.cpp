#include <eigen3/Eigen/Dense>
#include <iostream>
#include <functional>

#include "../include/solver.h"

Eigen::Vector2d BoundaryValueProblem::Shoot(Eigen::Vector2d u) {
    y[0] = u(0);
    for(int i = 0; i < N; i++) {
        Eigen::Vector2d k1, k2, k3, k4;
        k1 = h * Psi(x[i], u);
        k2 = h * Psi(x[i] + 0.5 * h, u + 0.5 * k1);
        k3 = h * Psi(x[i] + 0.5 * h, u + 0.5 * k2);
        k4 = h * Psi(x[i] + h, u + k3);

        u += (k1 + 2 * k2 + 2 * k3 + k4) / 6;
        y[i + 1] = u(0);
    }

    return u;
}

double FindZero(std::function<double(double)> Func, double l, double r, double f_l, double f_r, double eps) {
    if(l >= r) {
        std::swap(l, r);
        std::swap(f_l, f_r);
    }

    while(std::abs(r - l) > eps) {
        double m = (l + r) / 2;
        double f_m = Func(m);
        if((f_l > 0) == (f_m > 0)) {
            l = m; f_l = f_m;
        } 
        else {
            r = m; f_r = f_m;
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
    int i = 1;
    while(std::abs(i) < std::pow(2, max_iter - 1) && (f_l > 0) == (f_r > 0)) {
        i = - (i + (i > 0));
        f_r = Funk(i);
    }
    return Ivals(FindZero(Funk, 0, i, f_l, f_r, eps));
}

Eigen::Vector2d BoundaryValueProblem::FindInitialValsLinear() {
    std::function<Eigen::Vector2d(double)> Ivals;
    if(b_1 == 0) {
        Ivals = [this](double x) { return Eigen::Vector2d(c_1, x);};
    }
    else {
        Ivals = [this](double x) { return Eigen::Vector2d(x, (c_1 - a_1 * x) / b_1); };
    }

    auto Funk = [this, Ivals](double x){ Eigen::Vector2d fvals = Shoot(Ivals(x));
                                         return a_2 * fvals(0) + b_2 * fvals(1) - c_2; };

    double f_0 = Funk(0);
    double f_1 = Funk(1);

    return Ivals(f_0 / (f_0 - f_1));
}