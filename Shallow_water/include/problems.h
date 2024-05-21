#include <functional>
#include <cmath>
#include <string>

const double G_EARTH = 9.81;

double flat_gaussian(double x, double y);

struct FlatProblem {
    const double l_x = 1e5;
    const double l_y = 1e5;

    const double t = 1e5;

    std::function<double(double, double)> hInitial;

    std::string name;

    FlatProblem(std::function<double(double, double)> hInitial, std::string name) : hInitial(hInitial), name(name)
    { }
};