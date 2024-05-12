#include <functional>
#include <iostream>

#include "../include/problem.h"
#include "../include/progress_bar.h"

int main(int argc, char **argv)
{
    /*
    Setting up the problem
    */
    double alpha = 2; 
    double L = 10; double T = 10;
    double D = 0.9;
    
    auto u_x = []( double x ) { return 0; }; //initial values
    auto u_t = [alpha, D]( double t ) { return std::pow(alpha * D, 1 / alpha) * std::pow(D * t, 1 / alpha); }; //boundary condition
    auto exact_solution = [alpha, D](double x, double t) { if (x < D * t) return std::pow(alpha * D, 1 / alpha) * std::pow(D * t - x, 1 / alpha);
                                                           else return 0.0; };
    
    int n = 100; int k = 10000;
    double sigma = 0.5;

    auto problem = Problem(alpha, L, T, n, k, sigma, u_x, u_t, exact_solution);
    
    /*
    Testing
    */
    problem.solve();
    std::cout << problem.get_error() << '\n';
    problem.save_solution(100, 100);

    return 0;
}