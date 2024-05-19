#include <iostream>

#ifndef problems
    #define problems
    #include "../include/problems.h"
#endif
#include "../include/solver.h"
#ifndef calc_mesh
    #define calc_mesh
    #include "../include/calc_mesh.h"
#endif

int main()
{
    auto problem = FlatGaussianProblem();

    int n_x = 10; int n_y = 10; int k = 10;
    auto mesh = CalcMesh(problem, 10, 10);
    mesh.snapshot(0);

    return 0;
}