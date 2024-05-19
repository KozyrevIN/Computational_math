#include <functional>
#include <cmath>
#include <string>

struct FlatProblem {
    const double l_x = 1e6;
    const double l_y = 1e6;

    const double T = 1e6;

    virtual double hInital(double x, double y) { return 0.0; };

    virtual std::string name() { return ""; };
};

struct FlatGaussianProblem : FlatProblem {
    double hInitial(double x, double y) {
        double x_0 = l_x / 2; double y_0 = l_y / 2;
        return 100 + std::exp( -(std::pow(x - x_0, 2) + std::pow(y - y_0, 2)) / (2 * (l_x * l_x + l_y * l_y) / 100));
    }
    std::string name() {
        return "flat_gaussian";
    }
};