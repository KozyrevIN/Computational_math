//Функции "дифференцирования" матриц

Eigen::ArrayXXd derivative_x(Eigen::ArrayXXd m, double dx) {
    unsigned int r = m.rows();
    unsigned int c = m.cols();

    Eigen::ArrayXXd offset = Eigen::ArrayXXd::Zero(r, c);

    offset.block(0, 0, r - 1, c) = m.block(1, 0, r - 1, c);
    offset.block(1, 0, r - 1, c) -= m.block(0, 0, r - 1, c);

    offset.row(0) *= 2;
    offset.row(r - 1) *= 2;

    return offset / (2 * dx);
}

Eigen::ArrayXXd derivative_y(Eigen::ArrayXXd m, double dy) {
    unsigned int r = m.rows();
    unsigned int c = m.cols();

    Eigen::ArrayXXd offset = Eigen::ArrayXXd::Zero(r, c);

    offset.block(0, 0, r, c - 1) = m.block(0, 1, r, c - 1);
    offset.block(0, 1, r, c - 1) -= m.block(0, 0, r, c - 1);

    offset.col(0) *= 2;
    offset.col(c - 1) *= 2;

    return offset / (2 * dy);
}

// методы абстрактного класса Solver
template <typename Problem>
Solver<Problem>::Solver(Problem& problem, unsigned int n_1, unsigned int n_2, unsigned int k) : problem(problem), n_1(n_1), n_2(n_2), k(k) {
    mesh = CalcMesh(problem, n_1, n_2);

    u = Eigen::ArrayXXd::Zero(n_1, n_2);
    v = Eigen::ArrayXXd::Zero(n_1, n_2);
    h = Eigen::ArrayXXd::Zero(n_1, n_2);
}

template <typename Problem>
void Solver<Problem>::doStep() {
    //do nothing
}

template <typename Problem>
void Solver<Problem>::snapshot(unsigned int frame) {
    //do nothing
}

template <typename Problem>
void Solver<Problem>::solve(unsigned int num_frames) {
    ProgressBar bar((double) k, problem.name);

    if (num_frames != 0) {
        this -> snapshot(0);
    }

    for (unsigned int q = 1; q < k; q++) {
        this -> doStep();
        if (num_frames != 0 and (q % (k / num_frames) == 0)) {
            this -> snapshot(q / (k / num_frames));
        }
        bar.update_and_print_progress(1.0);
    }

    bar.set_100();
}

// Методы класса FlatSolver
FlatSolver::FlatSolver(FlatProblem& problem, unsigned int n_1, unsigned int n_2, unsigned int k) : Solver(problem, n_1, n_2, k) {
    for (unsigned int j = 0; j < n_2; j++) {
        for (unsigned int i = 0; i < n_1; i++) {
            h(i, j) = problem.hInitial((problem.l_x * i) / n_1, (problem.l_y * j) / n_2);
        }
    }

    dx = problem.l_x / (n_1 - 1);
    dy = problem.l_y / (n_2 - 1);
    dt = problem.t / (k - 1);
}

void FlatSolver::doStep() {
    Eigen::ArrayXXd u_new(n_1, n_2);
    Eigen::ArrayXXd v_new(n_1, n_2);
    Eigen::ArrayXXd h_new(n_1, n_2);

    u_new = u - dt * (u * derivative_x(u, dx) + v * derivative_y(u, dy) + G_EARTH * derivative_x(h, dx));
    v_new = v - dt * (u * derivative_x(v, dx) + v * derivative_y(v, dy) + G_EARTH * derivative_y(h, dx));
    h_new = h - dt * (derivative_x(u * h, dx) + derivative_y(v * h, dy));

    std::swap(u, u_new);
    std::swap(v, v_new);
    std::swap(h, h_new);
}

void FlatSolver::snapshot(unsigned int frame) {
    mesh.flatProject(u, v, h);
    mesh.snapshot(frame);
}

// Методы абстрактной фабрики
FlatSolver SolverFactory::getSolver(FlatProblem& problem, unsigned int n_1, unsigned int n_2, unsigned int k) {
    return FlatSolver(problem, n_1, n_2, k);
}

