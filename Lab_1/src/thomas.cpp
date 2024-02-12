#include <eigen3/Eigen/Dense>
#include <iostream>

#include "../include/solver.h"

void BoundaryValueProblem::Thomas() {
    std::vector<double> a(N + 1, 0), b(N + 1, 0), c(N + 1, 0), d(N + 1, 0);
    
    for(int i = 1; i < N; i++) {
        a[i] = (1 / h - p(x[i]) / 2) / h;
        c[i] = (1 / h + p(x[i]) / 2) / h;
        b[i] =  a[i] + c[i] - q(x[i]);
        d[i] = f(x[i]);
    }

    int k_1 = 0, k_2 = 0;
    if(b_1 == 0) {
        b[0] = -a_1;
        d[0] = c_1;
    }
    else {
        a[0] = (1 / h - p(x[0]) / 2) / h;
        c[0] = (1 / h + p(x[0]) / 2) / h;
        b[0] =  a[0] + c[0] - q(x[0]);
        d[0] = f(x[0]);

        a.insert(a.begin(), 0);
        b.insert(b.begin(), -b_1 / (2 * h));
        c.insert(c.begin(), -a_1);
        d.insert(d.begin(), c_1);

        double mult = b_1 / (2 * h * c[1]);
        b[0] += mult * a[1];
        c[0] += mult * b[1];
        d[0] -= mult * d[1];

        k_1 = 1;
    }

    if(b_2 == 0) {
        b[N + k_1] = -a_2;
        d[N + k_1] = c_2;
    }
    else {
        a[N + k_1] = (1 / h - p(x[N]) / 2) / h;
        c[N + k_1] = (1 / h + p(x[N]) / 2) / h;
        b[N + k_1] =  a[N + k_1] + c[N + k_1] - q(x[N]);
        d[N + k_1] = f(x[N]);

        a.push_back(a_2);
        b.push_back(-b_2 / (2 * h));
        c.push_back(0);
        d.push_back(c_2);

        double mult = -b_2 / (2 * h * a[N + k_1]);
        a[N + k_1 + 1] += mult * b[N + k_1];
        b[N + k_1 + 1] += mult * c[N + k_1];
        d[N + k_1 + 1] -= mult * d[N + k_1];

        k_2 = 1;
    }

    for(int i = 0; i <= N + k_1 + k_2; i++) {
        std::cout << a[i] << ' ' << b[i] << ' ' << c[i] << ' ' << d[i] << '\n';
    }
    std::cout << '\n';

    int N_new = N + k_1 + k_2;

    std::vector<double> P(N_new + 1, 0), Q(N_new + 1, 0);
    P[1] = c[0] / b[0]; Q[1] = - d[0] / b[0];
    for(int i = 1; i < N_new; i++) {
        P[i + 1] = c[i] / (b[i] - a[i] * P[i]);
        Q[i + 1] = (a[i] * Q[i] - d[i]) / (b[i] - a[i] * P[i]);
    }

    for(int i = 0; i <= N_new; i++) {
        std::cout << P[i] << ' ' << Q[i] << '\n';
    }

    std::vector<double> u(N_new + 1, 0);
    u[N_new] = (a[N_new] * Q[N_new] - d[N_new]) / (b[N_new] - a[N_new] * P[N_new]);
    for(int i = N_new - 1; i >= 0; i--) {
        u[i] = P[i + 1] * u[i + 1] + Q[i + 1];
    }

    for(int i = 0; i <= N; i++) {
        y[i] = u[i + k_1];
    }
}