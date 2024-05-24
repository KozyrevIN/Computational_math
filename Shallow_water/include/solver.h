#include <eigen3/Eigen/Dense>

#ifndef problems_header
    #define problems_header
    #include "problems.h"
#endif
#ifndef calc_mesh_header
    #define calc_mesh_header
    #include "calc_mesh.h"
#endif
#ifndef progress_bar_header
    #define progress_bar_header
    #include "progress_bar.h"
#endif

// Переменные в уравнении мелкой воды
struct Variables {
    Eigen::ArrayXXd u;
    Eigen::ArrayXXd v;
    Eigen::ArrayXXd h;

    Variables& operator+=(const Variables& var);

    Variables();
    Variables(unsigned int n_1, unsigned int n_2);
    Variables(const Eigen::ArrayXXd& u, const Eigen::ArrayXXd& v, const Eigen::ArrayXXd& h);
};

Variables operator+(const Variables& var_1, const Variables& var_2);

template <typename T>
Variables operator*(const T& m, const Variables& var);

template <typename T>
Variables operator/(const Variables& var, const T& div);

// Абстрактный класс солвер
template <typename Problem>
class Solver 
{
friend class SolverFactory;

protected:
    Problem problem;

    unsigned int n_1;
    unsigned int n_2;
    unsigned int k;

    double dt;

    Variables var;

    CalcMesh mesh;

    virtual Variables equation(const Variables& state);
    void doStep();

    Solver(const Problem& problem, unsigned int n_1, unsigned int n_2, unsigned int k);

public:
    void solve(unsigned int num_frames);
};

// Солвер для плоской задачи
Eigen::ArrayXXd derivative_x(const Eigen::ArrayXXd& m, double dx);

Eigen::ArrayXXd derivative_y(const Eigen::ArrayXXd& m, double dy);

class FlatSolver : public Solver<FlatProblem> 
{
friend class SolverFactory;

private:
    double dx;
    double dy;

    Variables equation(const Variables& state);

protected:
    FlatSolver(const FlatProblem& problem, unsigned int n_1, unsigned int n_2, unsigned int k);
};

//Солвер для сферической задачи
Eigen::ArrayXXd derivative_lambda(const Eigen::ArrayXXd& m, double r, double dlambda);

Eigen::ArrayXXd derivative_phi(const Eigen::ArrayXXd& m, double r, double dphi);

class SphericalSolver : public Solver<SphericalProblem> 
{
friend class SolverFactory;

private:
    double dlambda;
    double dphi;

    Eigen::ArrayXXd f;

    Variables equation(const Variables& state);

protected:
    SphericalSolver(const SphericalProblem& problem, unsigned int n_1, unsigned int n_2, unsigned int k);
};

// Фабрика солверов
class SolverFactory
{
public:
    SolverFactory() = default;
    FlatSolver getSolver(const FlatProblem& problem, unsigned int n_1, unsigned int n_2, unsigned int k);
    SphericalSolver getSolver(const SphericalProblem& problem, unsigned int n_1, unsigned int n_2, unsigned int k);
};

#include "solver.tpp"
#include "flat_solver.tpp"
#include "spherical_solver.tpp"