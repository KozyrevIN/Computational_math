#include <functional>

#include "../include/convection_diffusion_problem.h"

ConvectionDiffusionProblem::ConvectionDiffusionProblem(double a, double L, double T,
                                                         std::function<double(double)> f_0, std::function<double(double)> y_0, std::string name):
                                                         a{a}, L{L}, T{T}, f_0{f_0}, y_0{y_0}, name{name} { }

double ConvectionDiffusionProblem::y_exact(double x, double t) {
    if (x > a * t) {
        return(f_0(x - a * t));
    }
    else {
        return(y_0(t - x / a));
    }
}