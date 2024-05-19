#include <eigen3/Eigen/Dense>

#ifndef problems_header
    #define problems_header
    #include "problems.h"
#endif
#ifndef calc_mesh_header
    #define calc_mesh_header
    #include "calc_mesh.h"
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

    Eigen::MatrixXd u;
    Eigen::MatrixXd v;
    Eigen::MatrixXd h;

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