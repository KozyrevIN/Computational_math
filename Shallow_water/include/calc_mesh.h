#include <iostream>
#include <cmath>
#include <vector>
#include <eigen3/Eigen/Dense>

#include <iostream>

#ifndef problems
    #define problems
    #include "problems.h"
#endif

// Класс расчётной сетки
class CalcMesh
{
private:
    // 3D-сетка из расчётных точек
    unsigned int n_1;
    unsigned int n_2;

    Eigen::MatrixXd x;
    Eigen::MatrixXd y;
    Eigen::MatrixXd z;

    Eigen::MatrixXd vx;
    Eigen::MatrixXd vy;
    Eigen::MatrixXd vz;

    Eigen::MatrixXd h;

public:
    // Конструктор для плоской задачи
    CalcMesh(FlatProblem problem, unsigned n_x, unsigned n_y);

    // Метод отвечает за проектирование плоского решения на меш
    void flatProject(Eigen::MatrixXd h, Eigen::MatrixXd u, Eigen::MatrixXd v);

    // Метод отвечает за запись текущего состояния сетки в снапшот в формате VTK
    void snapshot(unsigned int snap_number, std::string problem_name);
};