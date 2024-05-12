#include <eigen3/Eigen/Dense>
#include <iostream>

#include "../include/problem.h"

Eigen::Matrix3Xd Problem::generate_system(const Eigen::VectorXd u) {

    Eigen::VectorXd u_alpha(n + 1);
    for (int i = 0; i <= n; i++) {
        u_alpha(i) = std::pow(u(i), alpha);
    }

    Eigen::Matrix3Xd system(3, n + 1);
    double c_1 = (1 - sigma) / (2 * h * h);
    double c_2 = 1 / tau;

    //a(i)
    system(0, 0) = 0;
    for (int i = 1; i <= n; i++) {
        system(0, i) = c_1 * (u_alpha(i) + u_alpha(i - 1));
    }

    //b(i)
    system(1, 0) = -1;
    for (int i = 1 ; i < n; i++) {
        system(1, i) = c_1 * (u_alpha(i + 1) + 2 * u_alpha(i) + u_alpha(i - 1)) + c_2;
    }
    system(1, n) = -1;

    //c(i)
    system(2, 0) = 0;
    for (int i = 1 ; i < n; i++) {
        system(2, i) = c_1 * (u_alpha(i + 1) + u_alpha(i));
    }
    system(2, n) = 0;

    return system;
}

Eigen::VectorXd Problem::generate_d(Eigen::VectorXd u, int j) {
    Eigen::VectorXd u_alpha(n + 1);
    for (int i = 0; i <= n; i++) {
        u_alpha(i) = std::pow(u(i), alpha);
    }

    Eigen::VectorXd d(n + 1);
    double c_1 = 1 / tau;
    double c_2 = sigma / (2 * h * h);
    
    d(0) = u_t((T * j) / k);
    for (int i = 1; i < n; i++) {
        d(i) = -c_1 * numerical_solution(i, j - 1) - 
                c_2 * ((u_alpha(i + 1) + u_alpha(i)) * (numerical_solution(i + 1, j - 1) - numerical_solution(i, j - 1)) -
                       (u_alpha(i) + u_alpha(i - 1)) * (numerical_solution(i, j - 1) - numerical_solution(i - 1, j - 1)));
    }
    d(n) = 0;
    
    return d;
}

Eigen::VectorXd Problem::thomas_solve(Eigen::Matrix3Xd system, Eigen::VectorXd d) {
    Eigen::VectorXd P(n + 1), Q(n + 1);

    P(1) = system(2, 0) / system(1, 0); Q(1) = - d(0) / system(1, 0);
    for(int i = 1; i < n; i++) {
        P(i + 1) = system(2, i) / (system(1, i) - system(0, i) * P(i));
        Q(i + 1) = (system(0, i) * Q(i) - d(i)) / (system(1, i) - system(0, i) * P(i));
    }

    Eigen::VectorXd u(n + 1);
    u(n) = (system(0, n) * Q(n) - d(n)) / (system(1, n) - system(0, n) * P(n));
    for(int i = n - 1; i >= 0; i--) {
        u(i) = P(i + 1) * u(i + 1) + Q(i + 1);
    }

    return u;
}

void Problem::solve() {
    Eigen::VectorXd u;
    Eigen::Matrix3Xd system;
    Eigen::VectorXd d;

    for (int i = 0; i <= n; i++) {
        numerical_solution(i, 0) = u_x((L * i) / n);
    }
    
    for (int j = 1; j <= k; j++) {
        u = numerical_solution.col(j - 1);
        for (int l = 0; l < 2; l++) {
            system = this -> generate_system(u);
            d = this -> generate_d(u, j);
            u = thomas_solve(system, d);
        }
        numerical_solution.col(j) = u;
    }
}