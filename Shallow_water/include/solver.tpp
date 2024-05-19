// методы абстрактного класса Solver
template <typename Problem>
Solver<Problem>::Solver(Problem& problem, unsigned int n_1, unsigned int n_2, unsigned int k) : problem(problem), n_1(n_1), n_2(n_2), k(k) {
    mesh = CalcMesh(problem, n_1, n_2);

    u = Eigen::MatrixXd(n_1, n_2);
    v = Eigen::MatrixXd(n_1, n_2);
    h = Eigen::MatrixXd(n_1, n_2);
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
    if (num_frames != 0) {
        this -> snapshot(0);
    }
    for (unsigned int q = 1; q < k; q++) {
        this -> doStep();
        if (num_frames != 0 and (q % (k / num_frames) == 0)) {
            this -> snapshot(q / (k / num_frames));
        }
    }
}

// Методы класса FlatSolver
FlatSolver::FlatSolver(FlatProblem& problem, unsigned int n_1, unsigned int n_2, unsigned int k) : Solver(problem, n_1, n_2, k) {
    for (unsigned int j = 0; j < n_2; j++) {
        for (unsigned int i = 0; i < n_1; i++) {
            h(i, j) = problem.hInitial((problem.l_x * i) / n_1, (problem.l_y * j) / n_2);
        }
    }
}

void FlatSolver::doStep() {
    //do nothing
}

void FlatSolver::snapshot(unsigned int frame) {
    mesh.flatProject(u, v, h);
    mesh.snapshot(frame);
}

// Методы абстрактной фабрики
FlatSolver SolverFactory::getSolver(FlatProblem& problem, unsigned int n_1, unsigned int n_2, unsigned int k) {
    return FlatSolver(problem, n_1, n_2, k);
}