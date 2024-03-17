#include <iostream>
#include <math.h>
#include <mpi.h>

#ifndef convection_diffusion_problem
#define convection_diffusion_problem
    #include "../include/convection_diffusion_problem.h"
#endif
#ifndef solver
#define solver
    #include "../include/solver.h"
#endif
#include "../include/schemes.h"


int main(int argc, char **argv)
{
    //initializing MPI
    MPI_Init(&argc, &argv);

    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    double a = 0.5; //parameter in equation
    double L = 1; double T = 1; //space and time sizes of study area
    
    //Setting up the rectangle problem
    std::string name_rect = "rectangle";
    auto f_0_rect = []( double x ) { if (std::abs(x - 0.25) < 0.1) return 1.0;
                                else return 0.0; }; //initial values
    auto y_0_rect = []( double t ) { return 0; }; //boundary conditions

    auto rectangle = ConvectionDiffusionProblem(a, L, T, f_0_rect, y_0_rect, name_rect);

    //setting up hat problem
    std::string name_hat = "hat";
    auto f_0_hat = []( double x ) { if (std::abs(x - 0.25) < 0.1) return std::exp(- 1 / (1 - 100 * std::pow(x - 0.25, 2)));
                                else return 0.0; }; //initial values
    auto y_0_hat = []( double t ) { return 0; }; //boundary conditions

    auto hat = ConvectionDiffusionProblem(a, L, T, f_0_hat, y_0_hat, name_hat);

    //LeftAngleSolver------------------------------------------------------------------------------------------------------------------

    //setting up solvers
    int N = 800; int K = 800;
    auto rect_solver = LeftAngleSolver(size, rank, N, K, rectangle);
    auto hat_solver = LeftAngleSolver(size, rank, N, K, hat);

    //testing solver
    rect_solver.solve(180);
    hat_solver.solve(180);

    MPI_Finalize();

    return 0;
}
