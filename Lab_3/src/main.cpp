#include <functional>
#include <iostream>
#include <fstream>

#include "../include/problem.h"
#ifndef progress_bar
#define progress_bar
    #include "../include/progress_bar.h"
#endif

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
    
    /*
    Testing
    */
    int n_max = 256; int k_max = n_max * n_max;
    double sigma = 0.5;

    //initializing progress bar
    double total_points = (n_max + 1) * (k_max + 1);
    for (int n = 8; n <= n_max; n = n * 2) {
        total_points += std::pow(n + 1, 3);
    }
    ProgressBar bar(64, total_points);

    //plotting 1 solution
    auto problem = Problem(alpha, L, T, n_max, k_max, sigma, u_x, u_t, exact_solution);
    problem.solve(&bar);
    problem.save_solution(128, 128);
    
    //calculating error
    
    for (int n = 8; n <= n_max; n = n * 2) {
        int k = n * n;
        problem.change_n_k(n, k);
        problem.solve(&bar);
        double error = problem.get_error();
    }
    
    bar.set_100();
    return 0;
}