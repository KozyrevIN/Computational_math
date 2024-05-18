#include <iostream>

#ifndef problems
    #define problems
    #include "../include/problems.h"
#endif
#include "../include/solver.h"

int main()
{
    auto problem = FlatGaussianProblem();

    int n_x = 1000; int n_y = 1000; int k = 1000;
    auto solver = Solver(problem, n_x, n_y, k);

    unsigned int num_frames = 100;
    solver.solve(num_frames);

    return 0;
}