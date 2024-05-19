#include <eigen3/Eigen/Dense>

#ifndef calc_mesh
    #define calc_mesh
    #include "calc_mesh.h"
#endif

template <typename Problem>

class Solver 
{
private:
    Problem problem;

    unsigned int n_1;
    unsigned int n_2;
    unsigned int k;
    CalcMesh mesh;

    Eigen::VectorXd u;
    Eigen::VectorXd v;
    Eigen::VectorXd h;

    void doStep();
    void snapshot(unsigned int frame);

public:
    Solver(Problem& problem, unsigned int n_1, unsigned int n_2, unsigned int k);

    void solve(unsigned int num_frames);
};

#include "solver.tpp"