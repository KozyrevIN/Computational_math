#include <functional>
#include <cmath>
#include <string>

const double G_EARTH = 9.81;

// Плоская задача
double flat_gaussian(double x, double y, double l_x, double l_y);

struct FlatProblem {
    const double l_x = 1;
    const double l_y = 1;
    const double t = 0.1;

    const std::function<double(double x, double y, double l_x, double l_y)> hFunction;
    double hInitial(double x, double y);

    const std::string name;

    FlatProblem(std::function<double(double, double, double, double)> hFunction, std::string name);
};

//Сферическая задача
double spherical_gaussian(double lambda, double phi);

struct SphericalProblem {
    const double r = 1;
    const double omega;
    const double t = 0.1;

    const std::function<double(double, double)>  hInitial;

    const std::string name;

    SphericalProblem(std::function<double(double, double)> hInitial, std::string name);
};