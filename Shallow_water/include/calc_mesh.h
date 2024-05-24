#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include <eigen3/Eigen/Dense>
#include <iostream>

#include <vtkDoubleArray.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkQuad.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkUnstructuredGrid.h>
#include <vtkSmartPointer.h>

#ifndef problems_header
    #define problems_header
    #include "problems.h"
#endif

// Класс расчётной сетки
class CalcMesh
{
private:
    // Размеры сетки
    unsigned int n_1;
    unsigned int n_2;

    // Название задачи, для которой строится сетка
    std::string problemName;

    // Координаты точек сетки
    Eigen::ArrayXX<Eigen::Vector3d> points;

    // Локальные орты
    Eigen::ArrayXX<Eigen::Matrix3d> localUnitVectors;

    // Scale factor для отображения поля h на меше
    double scaleFactor;

    // Дефолтная высота воды
    double hDefault;

public:
    // Пустой конструктор
    CalcMesh();

    // Конструктор для плоской задачи
    CalcMesh(FlatProblem& problem, unsigned int n_x, unsigned int n_y);

    // Конструктор для сферической задачи
    CalcMesh(SphericalProblem& problem, unsigned int n_lambda, unsigned int n_phi);

    // Метод отвечает за проектирование текущего состояния сетки в снапшот в формате VTK
    void snapshot(const Eigen::ArrayXXd& u, const Eigen::ArrayXXd& v, const Eigen::ArrayXXd& h, unsigned int snap_number);
};