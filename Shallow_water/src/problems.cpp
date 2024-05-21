#ifndef problems
    #define problems
    #include "../include/problems.h"
#endif

double flat_gaussian(double x, double y) {
    double l_x = 1e5; double l_y = 1e5;
    double x_0 = l_x / 2; double y_0 = l_y / 2;
    return std::exp( -(std::pow(x - x_0, 2) + std::pow(y - y_0, 2)) / (2 * (l_x * l_x + l_y * l_y) / 100));
}
