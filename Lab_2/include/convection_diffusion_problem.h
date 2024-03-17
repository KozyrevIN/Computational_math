#include <functional>
#include <string>

struct ConvectionDiffusionProblem //equation of form u_t + au_x = 0 with initial values f_0 and boundary conditions y_0
{
    double a; //parameter in equation
    double L, T; //space and time sizes of study area
    
    std::function<double(double)> f_0; //initial values
    std::function<double(double)> y_0; //boundary conditions

    std::string name;

    //constructor
    ConvectionDiffusionProblem(double, double, double, std::function<double(double)>, std::function<double(double)>, std::string name);

    //method that gives exact solution y(x, t)
    double y_exact(double, double);
};
