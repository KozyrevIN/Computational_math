#include <functional>
#include <eigen3/Eigen/Dense>

class Problem {
private:
    double alpha, L, T;
    int n, k;
    double h, tau;
    double sigma;

    std::function<double(double)> u_x;
    std::function<double(double)> u_t;
    std::function<double(double, double)> exact_solution;

    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> numerical_solution;

    Eigen::Matrix3Xd generate_system(Eigen::VectorXd);
    Eigen::VectorXd generate_d(Eigen::VectorXd, int);
    Eigen::VectorXd thomas_solve(Eigen::Matrix3Xd, Eigen::VectorXd);

public:
    Problem(double, double, double, int, int, double,
            std::function<double(double)>,
            std::function<double(double)>,
            std::function<double(double, double)>);

    void change_n_k(int, int);

    double get_error();
    void save_solution(int, int);

    void solve();
};