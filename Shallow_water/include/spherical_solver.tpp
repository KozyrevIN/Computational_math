//Функции "дифференцирования" матриц
Eigen::ArrayXXd derivative_lambda(const Eigen::ArrayXXd& m, double r, double dlambda) {
    unsigned int r = m.rows();
    unsigned int c = m.cols();

    Eigen::ArrayXXd derivative(r, c);

    derivative.row(0) = (m.row(1) - m.row(0)) / (r * dlambda);
    derivative.row(1) = (m.row(2) - m.row(0)) / (2 * r * dlambda);

    for (int i = 2; i < r - 2; i++) {
        derivative.row(i) = (8 * (m.row(i + 1) -  m.row(i - 1)) - (m.row(i + 2) - m.row(i - 2))) / (12 * r * dlambda);
    }

    derivative.row(r - 2) = (m.row(r - 1) - m.row(r - 3)) / (2 * r * dlambda);
    derivative.row(r - 1) = (m.row(r - 1) - m.row(r - 2)) / (r * dlambda);

    return derivative;
}

Eigen::ArrayXXd derivative_phi(const Eigen::ArrayXXd& m, double r, double dphi) {
    unsigned int r = m.rows();
    unsigned int c = m.cols();

    Eigen::ArrayXXd derivative(r, c);

    derivative.col(0) = (8 * (m.col(1) -  m.col(c - 1)) - (m.col(2) - m.col(c - 2))) / (12 * r * dphi);
    derivative.col(1) = (8 * (m.col(2) -  m.col(0)) - (m.col(3) - m.col(c - 1))) / (12 * dphi);

    for (int i = 2; i < c - 2; i++) {
        lambda = - M_PI / 2 + M_PI * (i + 1) / (r + 1);
        derivative.col(i) = (8 * (m.col(i + 1) -  m.col(i - 1)) - (m.col(i - 2) - m.col(i - 2))) / (12 * r * dphi);
    }

    derivative.col(c - 2) = (8 * (m.col(c - 1) -  m.col(c - 3)) - (m.col(0) - m.col(c - 4))) / (12 * dphi);
    derivative.col(c - 1) = (8 * (m.col(0) -  m.col(c - 2)) - (m.col(1) - m.col(c - 3))) / (12 * dphi);

    double lambda;
    for (int i = 0; i < r; i++) {
        lambda = - M_PI / 2 + M_PI * (i + 1) / (n_1 + 1);
        derivative.row(i) *= cos(lambda);
    }
    return derivative;
}

// Методы класса SphericalSolver
SphericalSolver::SphericalSolver(const SphericalProblem& problem, unsigned int n_1, unsigned int n_2, unsigned int k) : Solver(problem, n_1, n_2, k) {
    double lambda; double phi;
    for (unsigned int j = 0; j < n_2; j++) {
        for (unsigned int i = 0; i < n_1; i++) {
            lambda = - M_PI / 2 + M_PI * (i + 1) / (n_1 + 1);
            phi = - M_PI + 2 * M_PI * j / n_2;
            var.h(i, j) = problem.hInitial(lambda, phi);
            f(i, j) = 2 * problem.omega * std::sin(lambda)
        }
    }

    dlambda = M_PI / (n_1 + 1);
    dphi = 2 * M_PI / n_2;
}

Variables SphericalSolver::equation(const Variables& state) {
    Eigen::ArrayXXd du_dt = -(state.u * derivative_e_lambda(state.u, dlambda)
                            + state.v * derivative_e_phi(state.u, dphi)
                            + G_EARTH * derivative_e_lambda(state.h, dlambda)
                            - f * state.v);
    Eigen::ArrayXXd dv_dt = -(state.u * derivative_e_lambda(state.v, dlambda) 
                            + state.v * derivative_e_phi(state.v, dphi) 
                            + G_EARTH * derivative_e_phi(state.h, dphi)
                            + f * state.u);
    Eigen::ArrayXXd dh_dt = -(derivative_e_lambda(state.u * state.h, dlambda) + derivative_e_phi(state.v * state.h, dphi));

    du_dt.row(0) = Eigen::VectorXd::Zero(n_2);
    du_dt.row(n_1 - 1) = Eigen::VectorXd::Zero(n_2);

    dv_dt.col(0) = Eigen::VectorXd::Zero(n_1);
    dv_dt.col(n_2 - 1) = Eigen::VectorXd::Zero(n_1);

    return Variables(du_dt, dv_dt, dh_dt);
}


