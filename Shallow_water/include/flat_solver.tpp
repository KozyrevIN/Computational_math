//Функции "дифференцирования" матриц
Eigen::ArrayXXd derivative_x(const Eigen::ArrayXXd& m, double dx) {
    unsigned int r = m.rows();
    unsigned int c = m.cols();

    Eigen::ArrayXXd derivative(r, c);

    derivative.row(0) = (m.row(1) - m.row(0)) / dx;
    derivative.row(1) = (m.row(2) - m.row(0)) / (2 * dx);

    for (int i = 2; i < r - 2; i++) {
        derivative.row(i) = (8 * (m.row(i + 1) -  m.row(i - 1)) - (m.row(i + 2) - m.row(i - 2))) / (12 * dx);
    }

    derivative.row(r - 2) = (m.row(r - 1) - m.row(r - 3)) / (2 * dx);
    derivative.row(r - 1) = (m.row(r - 1) - m.row(r - 2)) / dx;

    return derivative;
}

Eigen::ArrayXXd derivative_y(const Eigen::ArrayXXd& m, double dy) {
    unsigned int r = m.rows();
    unsigned int c = m.cols();

    Eigen::ArrayXXd derivative(r, c);

    derivative.col(0) = (m.col(1) - m.col(0)) / dy;
    derivative.col(1) = (m.col(2) - m.col(0)) / (2 * dy);

    for (int i = 2; i < c - 2; i++) {
        derivative.col(i) = (8 * (m.col(i + 1) -  m.col(i - 1)) - (m.col(i + 2) - m.col(i - 2))) / (12 * dy);
    }

    derivative.col(r - 2) = (m.col(r - 1) - m.col(r - 3)) / (2 * dy);
    derivative.col(r - 1) = (m.col(r - 1) - m.col(r - 2)) / dy;

    return derivative;
}

// Методы класса FlatSolver
FlatSolver::FlatSolver(const FlatProblem& problem, unsigned int n_1, unsigned int n_2, unsigned int k) : Solver(problem, n_1, n_2, k) {
    for (unsigned int j = 0; j < n_2; j++) {
        for (unsigned int i = 0; i < n_1; i++) {
            var.h(i, j) = problem.hInitial((problem.l_x * i) / (n_1 - 1), (problem.l_y * j) / (n_2 - 1));
        }
    }

    dx = problem.l_x / (n_1 - 1);
    dy = problem.l_y / (n_2 - 1);
}

Variables FlatSolver::equation(const Variables& state) {
    Eigen::ArrayXXd du_dt = -(state.u * derivative_x(state.u, dx) + state.v * derivative_y(state.u, dy) + G_EARTH * derivative_x(state.h, dx));
    Eigen::ArrayXXd dv_dt = -(state.u * derivative_x(state.v, dx) + state.v * derivative_y(state.v, dy) + G_EARTH * derivative_y(state.h, dy));
    Eigen::ArrayXXd dh_dt = -(derivative_x(state.u * state.h, dx) + derivative_y(state.v * state.h, dy));

    du_dt.row(0) = Eigen::VectorXd::Zero(n_2);
    du_dt.row(n_1 - 1) = Eigen::VectorXd::Zero(n_2);

    dv_dt.col(0) = Eigen::VectorXd::Zero(n_1);
    dv_dt.col(n_2 - 1) = Eigen::VectorXd::Zero(n_1);

    return Variables(du_dt, dv_dt, dh_dt);
}

