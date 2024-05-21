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

template <typename Problem>
class Solver 
{
friend class SolverFactory;

protected:
    Problem problem;

    unsigned int n_1;
    unsigned int n_2;
    unsigned int k;
    CalcMesh mesh;

    Eigen::ArrayXXd u;
    Eigen::ArrayXXd v;
    Eigen::ArrayXXd h;

    virtual void doStep();
    virtual void snapshot(unsigned int frame);

    Solver(Problem& problem, unsigned int n_1, unsigned int n_2, unsigned int k);

public:
    void solve(unsigned int num_frames);
};

class FlatSolver : public Solver<FlatProblem> 
{
friend class SolverFactory;

private:
    double dt;
    double dx;
    double dy;

    void doStep();
    void snapshot(unsigned int frame);

protected:
    FlatSolver(FlatProblem& problem, unsigned int n_1, unsigned int n_2, unsigned int k);
};

class SolverFactory
{
public:
    SolverFactory() = default;
    FlatSolver getSolver(FlatProblem& problem, unsigned int n_1, unsigned int n_2, unsigned int k);
};

#include "solver.tpp"