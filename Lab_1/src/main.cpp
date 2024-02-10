#include <iostream>
#include <fstream>
#include <math.h>
#include <eigen3/Eigen/Dense>

#include "solver.h"

using namespace std;

//Boundary value problem parameters
double p(double x) {return -1;};
double q(double x) {return 0;};
double f(double x) {return 1 + x;};
double x_0 = 0, x_1 = 1; int N = 10000;
double a_1 = 1, b_1 = 0, c_1 = 1;
double a_2 = 1, b_2 = 0, c_2 = 1;

int main()
{
    BoundaryValueProblem problem(p, q, f, x_0, x_1, N, a_1, b_1, c_1, a_2, b_2, c_2);
    problem.FindInitialVals(1e-9, 1000);

    auto res = problem.GetResults(1000);
    ofstream out_shoot;
    out_shoot.open ("../out_shoot.csv");
    out_shoot << 'x' << ',' << 'y' << '\n';
    for(int i = 0; i <= 1000; i++) {
        out_shoot << res[0][i] << ',' << res[1][i] << '\n';
    }
    out_shoot.close();

    return 0;
}