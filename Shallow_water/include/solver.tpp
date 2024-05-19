
template <typename Problem>
Solver<Problem>::Solver(Problem& problem, unsigned int n_1, unsigned int n_2, unsigned int k) : problem(problem), n_1(n_1), n_2(n_2), k(k) {
    auto mesh = CalcMesh(problem, n_1, n_2);

    auto u = Eigen::MatrixXd(n_1, n_2);
    auto v = Eigen::MatrixXd(n_1, n_2);
    auto h = Eigen::MatrixXd(n_1, n_2);

    for (unsigned int j = 0; j < n_2; j++) {
        for (unsigned int i = 0; i < n_1; i++) {
            h(i, j) = problem.initial_h((problem.l_x * i) / n_1, (problem.l_y * j) / n_2);
        }
    }
}

template <typename Problem>
void Solver<Problem>::doStep() {
    //do nothing
}

template <typename Problem>
void Solver<Problem>::snapshot(unsigned int frame) {
    mesh.flatProject(u, v, h);
    mesh.snapshot(frame, problem.name);
}

template <typename Problem>
void Solver<Problem>::solve(unsigned int num_frames) {
    for (unsigned int q = 1; q < k; q++) {
        if (q % (k / num_frames) == 0) {
            this -> snapshot(q / (k / num_frames));
        }
        this -> doStep();
    }
}