#include <iostream>
#include <fstream>
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

template <typename T>

void plot_errors(T error_solver, int rank, int min_deg, int max_deg, bool stable) {
    double error;
    std::ofstream error_file;

    error_file.open("../out/errors/" + error_solver.get_scheme() + "/error_h_" + std::to_string(rank) + ".csv");
    error_file << "h,tau,error\n";

    error_solver.change_K(std::pow(2, max_deg));
    for(int n = std::pow(2, min_deg); n <= std::pow(2, max_deg - 2); n = n * 2) {
        error_solver.change_N(n);
        error = error_solver.solve(0);
        if (!isinf(error)) {
            error_file << error_solver.get_h() << ',' << error_solver.get_tau() << ',' << error << '\n';
        }
    }
    error_file.close();
    
    if (stable) {
        error_file.open("../out/errors/" + error_solver.get_scheme() + "/error_tau_" + std::to_string(rank) + ".csv");
        error_file << "h,tau,error\n";

        error_solver.change_N(std::pow(2, max_deg));
        for(int n = std::pow(2, min_deg); n <= std::pow(2, max_deg - 2); n = n * 2) {
            error_solver.change_K(n);
            error = error_solver.solve(0);
            if (!isinf(error)) {
                error_file << error_solver.get_h() << ',' << error_solver.get_tau() << ',' << error << '\n';
            }
        }
        error_file.close();
    }

    error_file.open("../out/errors/" + error_solver.get_scheme() + "/error_h_tau_" + std::to_string(rank) + ".csv");
    error_file << "h,tau,error\n";

    for(int n = std::pow(2, min_deg); n <= std::pow(2, max_deg); n = n * 2) {
        error_solver.change_N_K(n, n);
        error = error_solver.solve(0);
        if (!isinf(error)) {
            error_file << error_solver.get_h() << ',' << error_solver.get_tau() << ',' << error << '\n';
        }
    }
    error_file.close();
}

int main(int argc, char **argv)
{
    //initializing MPI
    MPI_Init(&argc, &argv);

    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    double a = 0.5; //parameter in equation
    double L = 4; double T = 4; //space and time sizes of study area
    
    //Setting up the rectangle problem
    std::string name_rect = "rectangle";
    auto f_0_rect = []( double x ) { if (std::abs(x - 1) < 0.5) return 1.0;
                                     else return 0.0; }; //initial values
    auto y_0_rect = []( double t ) { return 0; }; //boundary conditions

    auto rectangle = ConvectionDiffusionProblem(a, L, T, f_0_rect, y_0_rect, name_rect);

    //setting up hat problem
    std::string name_hat = "hat";
    
    auto f_0_hat = []( double x ) { if (std::abs(x - 1.25) < 1) return std::exp(- 1 / (1 - std::pow(x - 1.25, 2)));
                                else return 0.0; }; //initial values
    auto y_0_hat = []( double t ) { return 0; }; //boundary conditions

    auto hat = ConvectionDiffusionProblem(a, L, T, f_0_hat, y_0_hat, name_hat);

    //initial paraneters for solvers
    int N = 500; int K = 500;
    int deg_min = 5; int deg_max = 18;

    /*{ //LeftAngleSolver------------------------------------------------------------------------------------------------------------------

        //setting up solvers
        auto rect_solver = LeftAngleSolver(size, rank, N, K, rectangle);
        auto hat_solver = LeftAngleSolver(size, rank, N, K, hat);

        //testing solver
        rect_solver.solve(180);
        hat_solver.solve(180);

        //getting error
        plot_errors(hat_solver, rank, deg_min, deg_max, false);

        if (rank == 0) {
            std::cout << "finished left angle\n";
        }
    }

    { //RightAngleSolver------------------------------------------------------------------------------------------------------------------

        //setting up solvers
        auto rect_solver = RightAngleSolver(size, rank, N, K, rectangle);
        auto hat_solver = RightAngleSolver(size, rank, N, K, hat);

        //testing solver
        rect_solver.solve(180);
        hat_solver.solve(180);

        if (rank == 0) {
            std::cout << "finished right angle\n";
        }
    }

    { //ImplicitAngleSolver------------------------------------------------------------------------------------------------------------------

        //setting up solvers
        auto rect_solver = ImplicitAngleSolver(size, rank, N, K, rectangle);
        auto hat_solver = ImplicitAngleSolver(size, rank, N, K, hat);

        //testing solver
        rect_solver.solve(180);
        hat_solver.solve(180);

        //getting error
        plot_errors(hat_solver, rank, deg_min, deg_max, true);

        if (rank == 0) {
            std::cout << "finished implicit angle\n";
        }
    }

    { //FourPointSolver------------------------------------------------------------------------------------------------------------------

        //setting up solvers
        auto rect_solver = FourPointSolver(size, rank, N, K, rectangle);
        auto hat_solver = FourPointSolver(size, rank, N, K, hat);

        //testing solver
        rect_solver.solve(180);
        hat_solver.solve(180);

        if (rank == 0) {
            std::cout << "finished four point solver\n";
        }
    }

    { //LaxSolver------------------------------------------------------------------------------------------------------------------

        //setting up solvers
        auto rect_solver = LaxSolver(size, rank, N, K, rectangle);
        auto hat_solver = LaxSolver(size, rank, N, K, hat);

        //testing solver
        rect_solver.solve(180);
        hat_solver.solve(180);

        //getting errors
        plot_errors(hat_solver, rank, deg_min, deg_max, false);

        if (rank == 0) {
            std::cout << "finished lax\n";
        }
    }

    { //LaxWendroffSolver------------------------------------------------------------------------------------------------------------------

        //setting up solvers
        auto rect_solver = LaxWendroffSolver(size, rank, N, K, rectangle);
        auto hat_solver = LaxWendroffSolver(size, rank, N, K, hat);

        //testing solver
        rect_solver.solve(180);
        hat_solver.solve(180);

        //getting errors
        plot_errors(hat_solver, rank, deg_min, deg_max, false);

        if (rank == 0) {
            std::cout << "finished lax-wendroff\n";
        }
    }*/

    { //CrossSolver------------------------------------------------------------------------------------------------------------------

        //setting up solvers
        auto rect_solver = CrossSolver(size, rank, N, K, rectangle);
        auto hat_solver = CrossSolver(size, rank, N, K, hat);

        //testing solver
        rect_solver.solve(180);
        hat_solver.solve(180);

        //getting errors
        plot_errors(hat_solver, rank, deg_min, deg_max, false);

        if (rank == 0) {
            std::cout << "finished cross\n";
        }
    }

    MPI_Finalize();

    return 0;
}
