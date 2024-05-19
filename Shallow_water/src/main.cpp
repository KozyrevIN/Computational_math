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
    auto problem = FlatProblem(flat_gaussian, "flat_gaussian");
    unsigned int n_x = 20; 
    unsigned int n_y = 10; 
    unsigned int k = 10;
    unsigned int frames = 10;

    // Дальше код работает автоматически
    auto factory = SolverFactory();
    auto solver = factory.getSolver(problem, n_x, n_y, k);
    solver.solve(frames);

    return 0;
}