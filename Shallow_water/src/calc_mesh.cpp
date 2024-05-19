#include <iostream>
#include <cmath>
#include <eigen3/Eigen/Dense>

#include <vtkDoubleArray.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkQuad.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkUnstructuredGrid.h>
#include <vtkSmartPointer.h>

#ifndef calc_mesh
    #define calc_mesh
    #include "../include/calc_mesh.h"
#endif

CalcMesh::CalcMesh(FlatProblem problem, unsigned n_x, unsigned n_y) : n_1(n_x), n_2(n_y) {
    x = Eigen::MatrixXd(n_1, n_2);
    y = Eigen::MatrixXd(n_1, n_2);
    z = Eigen::MatrixXd(n_1, n_2);
    for (unsigned int j = 0; j < n_2; j++) {
        for (unsigned int i = 0; i < n_1; i++) {
            x(i, j) = (problem.l_x * i) / (n_1 - 1);
            y(i, j) = (problem.l_y * j) / (n_2 - 1);
            z(i, j) = 0;
        }
    }

    vx = Eigen::MatrixXd(n_1, n_2);
    vy = Eigen::MatrixXd(n_1, n_2);
    vz = Eigen::MatrixXd(n_1, n_2);

    h = Eigen::MatrixXd(n_1, n_2);

    problemName = problem.name;
}

void CalcMesh::flatProject(Eigen::MatrixXd& u, Eigen::MatrixXd& v, Eigen::MatrixXd& h)
{
    this -> h = h;
    vx = u;
    vy = v;

}

void CalcMesh::snapshot(unsigned int snap_number) {
    // Сетка в терминах VTK
    vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
    // Точки сетки в терминах VTK
    vtkSmartPointer<vtkPoints> dumpPoints = vtkSmartPointer<vtkPoints>::New();

    // Скалярное поле на точках сетки
    auto h_field = vtkSmartPointer<vtkDoubleArray>::New();
    h_field -> SetName("h");

    // Векторное поле на точках сетки
    auto vel = vtkSmartPointer<vtkDoubleArray>::New();
    vel->SetName("velocity");
    vel->SetNumberOfComponents(3);

    // Обходим все точки нашей расчётной сетки
    for(unsigned int i = 0; i < n_1; i++) {
        for(unsigned int j = 0; j < n_2; j++) {
            // Вставляем новую точку в сетку VTK-снапшота
            dumpPoints->InsertNextPoint(x(i, j), y(i, j), z(i, j));

            // Добавляем значение векторного поля в этой точке
            double _vel[3] = {vx(i, j), vy(i, j), vz(i, j)};
            vel -> InsertNextTuple(_vel);

            // И значение скалярного поля тоже
            h_field -> InsertNextValue(h(i, j));
        }
    }

    // Грузим точки в сетку
    unstructuredGrid->SetPoints(dumpPoints);

    // Присоединяем векторное и скалярное поля к точкам
    unstructuredGrid->GetPointData()->AddArray(vel);
    unstructuredGrid->GetPointData()->AddArray(h_field);

    // А теперь пишем, как наши точки объединены в квадраты
    for(unsigned int i = 0; i < n_1 - 1; i++) {
        for (unsigned int j = 0; j < n_2 - 1; j++) {
            auto quad = vtkSmartPointer<vtkQuad>::New();
            quad->GetPointIds()->SetId(0, n_2 * (i + 1) + j );
            quad->GetPointIds()->SetId(1, n_2 * (i + 1) + j + 1 );
            quad->GetPointIds()->SetId(2, n_2 * i + j + 1 );
            quad->GetPointIds()->SetId(3, n_2 * i + j );
            unstructuredGrid->InsertNextCell(quad->GetCellType(), quad->GetPointIds());
        }
    }

    // Создаём снапшот в файле с заданным именем
    std::string fileName = "../out/" + problemName + "_animation/" + std::to_string(snap_number) + ".vts";
    vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
    writer->SetFileName(fileName.c_str());
    writer->SetInputData(unstructuredGrid);
    writer->Write();
}
