#include <iostream>
#include <chrono>

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
    auto t1 = std::chrono::high_resolution_clock::now();

    // Задаем параметры задачи и разрешение для солвера
    auto problem = SphericalProblem(spherical_gaussian, "spherical_gaussian");
    unsigned int n_x = 500; 
    unsigned int n_y = 500; 
    unsigned int k = 45000;
    unsigned int frames = 900;

    // Дальше код работает автоматически
    auto factory = SolverFactory();
    auto solver = factory.getSolver(problem, n_x, n_y, k);
    solver.solve(frames);

    auto t2 = std::chrono::high_resolution_clock::now();
    std::cout << "runtime " << std::chrono::duration_cast<std::chrono::seconds>(t2 - t1).count() << " seconds\n";

    return 0;
}