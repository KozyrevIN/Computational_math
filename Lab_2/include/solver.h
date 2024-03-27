#include <vector>
#include <string>

#ifndef convection_diffusion_problem
#define convection_diffusion_problem
    #include "convection_diffusion_problem.h"
#endif

class Solver
{
    protected:
    std::string scheme; //solver scheme name

    int rank, size; //rank and size of MPI_Comm

    int N, K; //cells in computational mesh (n - space, k - time)
    int left, right, length; //vertexes passed to this process
    double h, tau, c, time_coef; //coordinate and time steps and current parameter

    std::vector<double> y_cur; //2 last calculation steps
    std::vector<double> y_prev; 
    std::vector<double> x;

    ConvectionDiffusionProblem problem; //problem to solve

    public:
    //constructor
    Solver(int, int, int, int, ConvectionDiffusionProblem);

    //common functions
    double get_error(double);
    double get_h();
    double get_tau();
    std::string get_scheme();
    virtual void change_N(int);
    void change_K(int);
    void change_N_K(int, int);
    double solve(int);
    //specific for scheme functions
    virtual void make_step(int);
};