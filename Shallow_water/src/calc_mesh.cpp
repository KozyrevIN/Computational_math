#include <iostream>
#include <cmath>
#include <eigen3/Eigen/Dense>

#include <vtkDoubleArray.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkXMLStructuredGridWriter.h>
#include <vtkStructuredGrid.h>
#include <vtkSmartPointer.h>

#ifndef calc_mesh
    #define calc_mesh
    #include "../include/calc_mesh.h"
#endif

CalcMesh::CalcMesh(FlatProblem& problem, unsigned n_x, unsigned n_y) : n_1(n_x), n_2(n_y) {
    auto x = Eigen::MatrixXd(n_1, n_2);
    auto y = Eigen::MatrixXd(n_1, n_2);
    auto z = Eigen::MatrixXd(n_1, n_2);
    for (unsigned int j = 0; j < n_2; j++) {
        for (unsigned int i = 0; i < n_1; i++) {
            x(i, j) = (problem -> l_x * i) / (n_1 - 1);
            y(i, j) = (problem -> l_y * j) / (n_2 - 1);
            z(i, j) = 0;
        }
    }

    auto vx = Eigen::MatrixXd(n_1, n_2);
    auto vy = Eigen::MatrixXd(n_1, n_2);
    auto vz = Eigen::MatrixXd(n_1, n_2);

    auto h = Eigen::MatrixXd(n_1, n_2);
}

CalcMesh::flatProject(Eigen::MatrixXd h, Eigen::MatrixXd u, Eigen::MatrixXd v) h(h), vx(u), vy(v)
{ }

void CalcMesh::snapshot(unsigned int snap_number, std::string problem_name) {
    // Сетка в терминах VTK
    vtkSmartPointer<vtkStructuredGrid> structuredGrid = vtkSmartPointer<vtkStructuredGrid>::New();
    // Точки сетки в терминах VTK
    vtkSmartPointer<vtkPoints> dumpPoints = vtkSmartPointer<vtkPoints>::New();

    // Скалярное поле на точках сетки
    auto h = vtkSmartPointer<vtkDoubleArray>::New();
    h->SetName("h");

    // Векторное поле на точках сетки
    auto vel = vtkSmartPointer<vtkDoubleArray>::New();
    vel->SetName("velocity");
    vel->SetNumberOfComponents(3);

    // Обходим все точки нашей расчётной сетки
    for(unsigned i = 0; i < n_1; i++) {
        for(unsigned j = 0; j < n_1; j++) {
            // Вставляем новую точку в сетку VTK-снапшота
            dumpPoints->InsertNextPoint(x(i, j), y(i, j), z(i, j));

            // Добавляем значение векторного поля в этой точке
            double _vel[3] = {vx(i, j), vy(i, j), vz(i, j)};
            vel->InsertNextTuple(_vel);

            // И значение скалярного поля тоже
            h->InsertNextValue(h(i, j));
        }
    }

    // Задаём размеры VTK-сетки (в точках, по трём осям)
    structuredGrid->SetDimensions(n_1, n_2, 1);
    // Грузим точки в сетку
    structuredGrid->SetPoints(dumpPoints);

    // Присоединяем векторное и скалярное поля к точкам
    structuredGrid->GetPointData()->AddArray(vel);
    structuredGrid->GetPointData()->AddArray(smth);

    // Создаём снапшот в файле с заданным именем
    string fileName = "../out/" + problem_name + "_animation-" + std::to_string(snap_number) + ".vts";
    vtkSmartPointer<vtkXMLStructuredGridWriter> writer = vtkSmartPointer<vtkXMLStructuredGridWriter>::New();
    writer->SetFileName(fileName.c_str());
    writer->SetInputData(structuredGrid);
    writer->Write();
}
