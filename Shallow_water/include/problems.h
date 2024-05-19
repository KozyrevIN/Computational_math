#include <functional>
#include <cmath>
#include <string>

double flat_gaussian(double x, double y);

struct FlatProblem {
    const double l_x = 1e6;
    const double l_y = 1e6;

    const double T = 1e6;

    std::function<double(double, double)> hInitial;

    std::string name;

    FlatProblem(std::function<double(double, double)> hInitial, std::string name) : hInitial(hInitial), name(name)
    { }
};