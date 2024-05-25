#include <iostream>

#ifndef problems_header
    #define problems_header
    #include "../include/problems.h"
#endif
#ifndef solver_header
    #define solver_header
    #include "../include/solver.h"
#endif
#ifndef calc_mesh_header
    #define calc_mesh_header
    #include "../include/calc_mesh.h"
#endif

int main()
{
    // Задаем параметры задачи и разрешение для солвера
    auto problem = SphericalProblem(spherical_gaussian, "spherical_gaussian_60fps");
    unsigned int n_x = 500; 
    unsigned int n_y = 500; 
    unsigned int k = 5000;
    unsigned int frames = 500;

    // Дальше код работает автоматически
    auto factory = SolverFactory();
    auto solver = factory.getSolver(problem, n_x, n_y, k);
    solver.solve(frames);

    return 0;
}