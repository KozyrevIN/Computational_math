#include <iostream>
#include <mpi.h>

#ifndef convection_diffusion_problem
#define convection_diffusion_problem
    #include "../include/convection_diffusion_problem.h"
#endif
#include "../include/solver.h"


int main(int argc, char **argv)
{
    //initializing MPI
    MPI_Init(&argc, &argv);

    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    //Setting up the problem
    double a = 0.5; //parameter in equation
    double L = 1; double T = 1; //space and time sizes of study area
    auto f_0 = []( double x ) { if (std::abs(x - 0.5) < 0.1) return 1.0;
                                else return 0.0; }; //initial values
    auto y_0 = []( double t ) { return 0; }; //boundary conditions

    auto rectangle = ConvectionDiffusionProblem(a, L, T, f_0, y_0);

    //setting up solver
    int N = 800; int K = 800;
    auto rect_solver = LeftAngleSolver(size, rank, N, K, rectangle);

    //testing solver
    double error;
    error = rect_solver.solve(180);
    std::cout << rank << ' ' << error << '\n';

    MPI_Finalize();

    return 0;
}
