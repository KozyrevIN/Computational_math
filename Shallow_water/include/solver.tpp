// Методы переменных
Variables::Variables(unsigned int n_1, unsigned int n_2) {
    u = Eigen::ArrayXXd::Zero(n_1, n_2);
    v = Eigen::ArrayXXd::Zero(n_1, n_2);
    h = Eigen::ArrayXXd::Zero(n_1, n_2);
}

Variables::Variables() = default;

Variables::Variables(const Eigen::ArrayXXd& u, const Eigen::ArrayXXd& v, const Eigen::ArrayXXd& h) : u(u), v(v), h(h)
{ }

Variables operator+(const Variables& var_1, const Variables& var_2) {
    return Variables(var_1.u + var_2.u, var_1.v + var_2.v, var_1.h + var_2.h);
}

template <typename T>
Variables operator*(const T& t, const Variables& var) {
    return Variables(t * var.u, t * var.v, t * var.h);
}

template <typename T>
Variables operator/(const Variables& var, const T& div) {
    return Variables(var.u / div, var.v / div, var.h / div);
}

Variables& Variables::operator+=(const Variables& var) {
    this -> u += var.u;
    this -> v += var.v;
    this -> h += var.h;
    return *this;
}

// методы абстрактного класса Solver
template <typename Problem>
Solver<Problem>::Solver(Problem& problem, unsigned int n_1, unsigned int n_2, unsigned int k) : problem(problem), n_1(n_1), n_2(n_2), k(k) {
    dt = problem.t / (k - 1);

    mesh = CalcMesh(problem, n_1, n_2);

    var = Variables(n_1, n_2);
}

template <typename Problem>
void Solver<Problem>::doStep() {
    Variables k1, k2, k3, k4;

    k1 = dt * equation(var);
    k2 = dt * equation(var + 0.5 * k1);
    k3 = dt * equation(var + 0.5 * k2);
    k4 = dt * equation(var + k3);

    var += (k1 + (2 * k2) + (2 * k3) + k4) / 6;
}

template <typename Problem>
Variables Solver<Problem>::equation(Variables state) {
    return Variables();
}

template <typename Problem>
void Solver<Problem>::solve(unsigned int num_frames) {
    ProgressBar bar((double) k, problem.name);

    if (num_frames != 0) {
        mesh.snapshot(u, v, h, 0);
    }

    for (unsigned int q = 1; q < k; q++) {
        doStep();
        if (num_frames != 0 and (q % (k / num_frames) == 0)) {
            mesh.snapshot(u, v, h, q / (k / num_frames));
        }
        bar.update_and_print_progress(1.0);
    }

    bar.set_100();
}

// Методы абстрактной фабрики
FlatSolver SolverFactory::getSolver(const FlatProblem& problem, unsigned int n_1, unsigned int n_2, unsigned int k) {
    return FlatSolver(problem, n_1, n_2, k);
}

SphericalSolver SolverFactory::getSolver(const SphericalProblem& problem, unsigned int n_1, unsigned int n_2, unsigned int k) {
    return SphericalSolver(problem, n_1, n_2, k);
}