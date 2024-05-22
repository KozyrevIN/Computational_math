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
    Variables(Eigen::ArrayXXd u, Eigen::ArrayXXd v, Eigen::ArrayXXd h);
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

    virtual Variables equation(Variables state);
    void doStep();
    virtual void snapshot(unsigned int frame);

    Solver(Problem& problem, unsigned int n_1, unsigned int n_2, unsigned int k);

public:
    void solve(unsigned int num_frames);
};

// Солвер для плоской задачи
class FlatSolver : public Solver<FlatProblem> 
{
friend class SolverFactory;

private:
    double dx;
    double dy;
    double dt;

    Variables equation(Variables state);
    void snapshot(unsigned int frame);

protected:
    FlatSolver(FlatProblem& problem, unsigned int n_1, unsigned int n_2, unsigned int k);
};

// Фабрика солверов
class SolverFactory
{
public:
    SolverFactory() = default;
    FlatSolver getSolver(FlatProblem& problem, unsigned int n_1, unsigned int n_2, unsigned int k);
};

#include "solver.tpp"
#include "flat_solver.tpp"