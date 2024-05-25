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
    double lambda_0 = M_PI / 4;
    double phi_0 = 0;
    double angle = std::acos(std::sin(lambda) * std::sin(lambda_0) + std::cos(lambda) * std::cos(lambda_0) * std::cos(phi - phi_0));
    return std::exp( - 2 * std::pow(angle, 2));
}

SphericalProblem::SphericalProblem(std::function<double(double, double)> hFunction, std::string name) : hFunction(hFunction), name(name) 
{ }

double SphericalProblem::hInitial(double lambda, double phi) {
    return hDefault + hFunction(lambda, phi);
}


