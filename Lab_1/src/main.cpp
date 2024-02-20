#include <iostream>
#include <fstream>
#include <math.h>
#include <eigen3/Eigen/Dense>

#include "../include/solver.h"

//Boundary value problem parameters y'' + py' + qy = f
double p(double x) { return 0; };
double q(double x) { return 4; };
double f(double x) { return -4 * x; };
double x_1 = 0, x_2 = 4; int N = 1000;
double a_1 = 1, b_1 = -1, c_1 = 2;
double a_2 = 1, b_2 = 0, c_2 = std::cos(8) - 4;

double exact_solution(double x) {return std::cos(2 * x) - x;}

int main()
{
    BoundaryValueProblem problem(p, q, f, x_1, x_2, N, a_1, b_1, c_1, a_2, b_2, c_2, exact_solution);

    //testing shooting method
    
    //nonlinear initial parameters serch

    auto ivals = problem.FindInitialVals(1e-10, 10000);
    problem.Shoot(ivals);
    auto res = problem.GetResults(1000);
    std::ofstream shoot;
    shoot.open ("../out/shoot.csv");
    shoot << 'x' << ',' << 'y' << '\n';

    for(int i = 0; i <= 1000; i++) {
        shoot << res[0][i] << ',' << res[1][i] << '\n';
    }

    shoot.close();


    //linear initial patameters search

    ivals = problem.FindInitialValsLinear();
    problem.Shoot(ivals);
    res = problem.GetResults(1000);
    std::ofstream shoot_linear;
    shoot_linear.open ("../out/shoot_linear.csv");
    shoot_linear << 'x' << ',' << 'y' << '\n';

    for(int i = 0; i <= 1000; i++) {
        shoot_linear << res[0][i] << ',' << res[1][i] << '\n';
    }

    shoot_linear.close();

    //getting error for shooting method depending on step size
    std::ofstream shoot_error;
    shoot_error.open ("../out/shoot_error.csv");
    shoot_error << 'h' << ',' << "error" << '\n';
    std::ofstream shoot_linear_error;
    shoot_linear_error.open ("../out/shoot_linear_error.csv");
    shoot_linear_error << 'h' <<',' << "error" << '\n';

    for(int i = 10; i <= 1e5; i += i / 10) {
        problem.CnangeN(i);
        double h = (x_2 - x_1) / i;

        ivals = problem.FindInitialVals(1e-14, 1000);
        problem.Shoot(ivals);
        double error = problem.GetError();
        shoot_error << h << ',' << error << '\n';

        ivals = problem.FindInitialValsLinear();
        problem.Shoot(ivals);
        error = problem.GetError();
        shoot_linear_error << h << ',' << error << '\n';
    }

    shoot_error.close();
    shoot_linear_error.close();
    
    //testing thomas method
    problem.Thomas();
    res = problem.GetResults(1000);
    std::ofstream thomas;
    thomas.open ("../out/thomas.csv");
    thomas << 'x' << ',' << 'y' << '\n';

    for(int i = 0; i <= 1000; i++) {
        thomas << res[0][i] << ',' << res[1][i] << '\n';
    }

    thomas.close();

    //getting error for thomas method depending on step size
    std::ofstream thomas_error;
    thomas_error.open ("../out/thomas_error.csv");
    thomas_error << 'h' << ',' << "error" << '\n';

    for(int i = 10; i <= 1e7; i += i / 10) {
        problem.CnangeN(i);
        double h = (x_2 - x_1) / i;

        problem.Thomas();
        double error = problem.GetError();
        thomas_error << h << ',' << error << '\n';
    }

    thomas_error.close();

    return 0;
}