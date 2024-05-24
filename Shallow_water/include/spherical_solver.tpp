//Функции "дифференцирования" матриц
Eigen::ArrayXXd derivative_e_lambda(const Eigen::ArrayXXd& m, double r, double dlambda) {
    unsigned int n_1 = m.rows();
    unsigned int n_2 = m.cols();

    Eigen::ArrayXXd derivative(n_1, n_2);

    derivative.row(0) = (m.row(1) - m.row(0)) / (r * dlambda);
    derivative.row(1) = (m.row(2) - m.row(0)) / (2 * r * dlambda);

    for (int i = 2; i < n_1 - 2; i++) {
        derivative.row(i) = (8 * (m.row(i + 1) -  m.row(i - 1)) - (m.row(i + 2) - m.row(i - 2))) / (12 * r * dlambda);
    }

    derivative.row(n_1 - 2) = (m.row(n_1 - 1) - m.row(n_1 - 3)) / (2 * r * dlambda);
    derivative.row(n_1 - 1) = (m.row(n_1 - 1) - m.row(n_1 - 2)) / (r * dlambda);

    return derivative;
}

Eigen::ArrayXXd derivative_e_phi(const Eigen::ArrayXXd& m, double r, double dphi) {
    unsigned int n_1 = m.rows();
    unsigned int n_2 = m.cols();

    Eigen::ArrayXXd derivative(n_1, n_2);

    derivative.col(0) = (8 * (m.col(1) -  m.col(n_2 - 1)) - (m.col(2) - m.col(n_2 - 2))) / (12 * r * dphi);
    derivative.col(1) = (8 * (m.col(2) -  m.col(0)) - (m.col(3) - m.col(n_2 - 1))) / (12 * dphi);

    double lambda;
    for (int i = 2; i < n_2 - 2; i++) {
        lambda = - M_PI / 2 + M_PI * (i + 1) / (r + 1);
        derivative.col(i) = (8 * (m.col(i + 1) -  m.col(i - 1)) - (m.col(i - 2) - m.col(i - 2))) / (12 * r * dphi);
    }

    derivative.col(n_2 - 2) = (8 * (m.col(n_2 - 1) -  m.col(n_2 - 3)) - (m.col(0) - m.col(n_2 - 4))) / (12 * dphi);
    derivative.col(n_2 - 1) = (8 * (m.col(0) -  m.col(n_2 - 2)) - (m.col(1) - m.col(n_2 - 3))) / (12 * dphi);

    for (int i = 0; i < r; i++) {
        lambda = - M_PI / 2 + M_PI * (i + 1) / (n_1 + 1);
        derivative.row(i) *= cos(lambda);
    }
    return derivative;
}

// Методы класса SphericalSolver
SphericalSolver::SphericalSolver(SphericalProblem& problem, unsigned int n_1, unsigned int n_2, unsigned int k) : Solver(problem, n_1, n_2, k) {
    double lambda; double phi;
    for (unsigned int j = 0; j < n_2; j++) {
        for (unsigned int i = 0; i < n_1; i++) {
            lambda = - M_PI / 2 + M_PI * (i + 1) / (n_1 + 1);
            phi = - M_PI + 2 * M_PI * j / n_2;
            var.h(i, j) = problem.hInitial(lambda, phi);
            f(i, j) = 2 * problem.omega * std::sin(lambda);
        }
    }

    dlambda = M_PI / (n_1 + 1);
    dphi = 2 * M_PI / n_2;
}

Variables SphericalSolver::equation(const Variables& state) {
    Eigen::ArrayXXd du_dt = -(state.u * derivative_e_lambda(state.u, problem.r, dlambda)
                            + state.v * derivative_e_phi(state.u, problem.r, dphi)
                            + G_EARTH * derivative_e_lambda(state.h, problem.r, dlambda)
                            - f * state.v);
    Eigen::ArrayXXd dv_dt = -(state.u * derivative_e_lambda(state.v, problem.r, dlambda) 
                            + state.v * derivative_e_phi(state.v, problem.r, dphi) 
                            + G_EARTH * derivative_e_phi(state.h, problem.r, dphi)
                            + f * state.u);
    Eigen::ArrayXXd dh_dt = -(derivative_e_lambda(state.u * state.h, problem.r, dlambda) + derivative_e_phi(state.v * state.h, problem.r, dphi));

    du_dt.row(0) = Eigen::VectorXd::Zero(n_2);
    du_dt.row(n_1 - 1) = Eigen::VectorXd::Zero(n_2);

    return Variables(du_dt, dv_dt, dh_dt);
}


