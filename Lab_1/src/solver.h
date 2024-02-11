#include <functional>
#include <eigen3/Eigen/Dense>

class BoundaryValueProblem //y'' + p(x)y' + q(x)y = f(x) with boundary conditions
{
private:
    std::function<double(double)> p;
    std::function<double(double)> q;
    std::function<double(double)> f;

    std::function<Eigen::Vector2d(double, Eigen::Vector2d)> Psi; //vector representation of equation: u' = psi(x, u) where u = (y, y')^T

    double x_1, x_2; //bounds
    int N; //number of vertixies in partition
    double h; //gridstep

    double a_1, b_1, c_1, a_2, b_2, c_2; //boundary conditions in form a_1*y(x_1) + b_1*y'(x_1) = c_1, a_2*y(x_2) + b_2*y'(x_2) = c_2

    double* x, * y; //arrays to store function calculated on last step

public:
    BoundaryValueProblem(std::function<double(double)>,
                         std::function<double(double)>, 
                         std::function<double(double)>, 
                         double, double, int,
                         double, double, double,
                         double, double, double);

    ~BoundaryValueProblem();

    //common methods
    double** GetResults(int n);

    //methods to deal with shooting method
    Eigen::Vector2d Shoot(Eigen::Vector2d);
    Eigen::Vector2d FindInitialVals(double, int);
    Eigen::Vector2cd FindALinear();

    //methods to deal with ...
};