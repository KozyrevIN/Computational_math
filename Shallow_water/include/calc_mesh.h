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
    // 3D-сетка из расчётных точек
    unsigned int n_1;
    unsigned int n_2;

    std::string problemName;

    vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid;

public:
    // Дефолтный конструктор
    CalcMesh() = default;
    // Конструктор для плоской задачи
    CalcMesh(FlatProblem problem, unsigned n_x, unsigned n_y);

    // Метод отвечает за проектирование плоского решения на меш
    void flatProject(Eigen::ArrayXXd& h, Eigen::ArrayXXd& u, Eigen::ArrayXXd& v);

    // Метод отвечает за запись текущего состояния сетки в снапшот в формате VTK
    void snapshot(unsigned int snap_number);
};