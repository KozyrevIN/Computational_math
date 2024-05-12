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
    int n_max = 512; int k_max = n_max * n_max;

    //initializing progress bar
    double total_points = 0;
    for (int n = 8; n <= n_max; n = n * 2) {
        total_points += 2 * std::pow(n + 1, 3);
    }
    total_points += 2 * std::pow(64 + 1, 3);
    total_points += 65 * 4097;
    ProgressBar bar(64, total_points);

    //initializing problem
    auto problem = Problem(alpha, L, T, 8, 64, 0, u_x, u_t, exact_solution);
    
    //sigma = 0
    problem.change_n_k(64, 512);
    problem.solve(&bar);
    problem.save_solution(64, 64);

    //calculating errors for sigma = 0
    std::ofstream output;
    output.open("../out/sigma_0.0/errors.csv");
    output << "h,error\n";
    for (int n = 8; n < n_max; n = n * 2) {
        int k = n * n;
        problem.change_n_k(n, k);
        problem.solve(&bar);
        double error = problem.get_error();
        output << L / n << ',' << error << '\n';
    }
    problem.change_n_k(n_max, k_max);
    problem.solve(&bar);
    double error = problem.get_error();
    output << L / n_max << ',' << error;
    output.close();

    //sigma = 0.5
    problem.change_sigma(0.5);
    problem.change_n_k(64, 512);
    problem.solve(&bar);
    problem.save_solution(64, 64);

    //calculating errors for sigma = 0.5
    output.open("../out/sigma_0.5/errors.csv");
    output << "h,error\n";
    for (int n = 8; n < n_max; n = n * 2) {
        int k = n * n;
        problem.change_n_k(n, k);
        problem.solve(&bar);
        double error = problem.get_error();
        output << L / n << ',' << error << '\n';
    }
    problem.change_n_k(n_max, k_max);
    problem.solve(&bar);
    error = problem.get_error();
    output << L / n_max << ',' << error;
    output.close();

    //sigma = 1
    problem.change_sigma(1);
    problem.change_n_k(64, 4096);
    problem.solve(&bar);
    problem.save_solution(64, 64);

    bar.set_100();
    return 0;
}