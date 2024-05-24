#ifndef problems
    #define problems
    #include "../include/problems.h"
#endif

//Плоская задача
double flat_gaussian(double x, double y, double l_x, double l_y) {
    double x_0 = l_x / 2; double y_0 = l_y / 2;
    return std::exp( -(std::pow(x - x_0, 2) + std::pow(y - y_0, 2)) / (2 * (l_x * l_x + l_y * l_y) / 100));
}

FlatProblem::FlatProblem(std::function<double(double, double, double, double)> hFunction, std::string name) : hFunction(hFunction), name(name) 
{ }

double FlatProblem::hInitial(double x, double y) {
    return hDefault + hFunction(x, y, l_x, l_y);
}

//Сферическая задача
double spherical_gaussian(double lambda, double phi) {
    double delta_lambda = lambda - M_PI / 4;
    double phi_0 = 0;
    double a = std::pow(std::sin(phi - phi_0 / 2.0), 2) + std::cos(phi) * std::cos(phi_0) * std::pow(std::sin(delta_lambda / 2.0), 2);
    double c = 2 * std::atan2(std::pow(a, 0.5), std::pow(1 - a, 0.5));
    return std::exp( -std::pow(c, 2));
}

SphericalProblem::SphericalProblem(std::function<double(double, double)> hInitial, std::string name) : hFunction(hFunction), name(name) 
{ }

double SphericalProblem::hInitial(double lambda, double phi) {
    return hDefault + hFunction(lambda, phi);
}


