#ifndef calc_mesh
    #define calc_mesh
    #include "../include/calc_mesh.h"
#endif

CalcMesh::CalcMesh() 
{ };

CalcMesh::CalcMesh(FlatProblem& problem, unsigned int n_x, unsigned int n_y) : n_1(n_x), n_2(n_y)
{
    problemName = problem.name;
    hDefault = problem.hDefault;
    scaleFactor = 1;

    points = Eigen::ArrayXX<Eigen::Vector3d>(n_1, n_2);
    for (int j = 0; j < n_2; j++) {
        for (int i = 0; i < n_1; i++) {
            points(i, j) << (problem.l_x * i) / (n_1 - 1),
                            (problem.l_y * j) / (n_2 - 1),
                            0;
        }
    }

    localUnitVectors = Eigen::ArrayXX<Eigen::Matrix3d>(n_1, n_2);
    for (int j = 0; j < n_2; j++) {
        for (int i = 0; i < n_1; i++) {
            localUnitVectors(i, j) << 1, 0, 0,
                                      0, 1, 0,
                                      0, 0, 1;
        }
    }
}

CalcMesh::CalcMesh(SphericalProblem& problem, unsigned int n_lambda, unsigned int n_phi) : n_1(n_lambda), n_2(n_phi)
{
    problemName = problem.name;
    hDefault = problem.hDefault;


    double lambda; double phi;
    double h_max = 0; double h;
    for (unsigned int j = 0; j < n_2; j++) {
        for (unsigned int i = 0; i < n_1; i++) {
            lambda = - M_PI / 2 + M_PI * (i + 1) / (n_1 + 1);
            phi = - M_PI + 2 * M_PI * j / n_2;
            h = std::abs(problem.hInitial(lambda, phi) - problem.hDefault);
            if (h > h_max) {
                h_max = h;
            }
        }
    }
    scaleFactor = 0 * 0.1 * problem.r / h_max;

    points = Eigen::ArrayXX<Eigen::Vector3d>(n_1, n_2);
    for (int j = 0; j < n_2; j++) {
        for (int i = 0; i < n_1; i++) {
            lambda = - M_PI / 2 + M_PI * (i + 1) / (n_1 + 1);
            phi = - M_PI + 2 * M_PI * j / n_2;
            points(i, j) <<   problem.r * std::cos(lambda) * std::cos(phi),
                            - problem.r * std::cos(lambda) * std::sin(phi),
                              problem.r * std::sin(lambda);
        }
    }

    localUnitVectors = Eigen::ArrayXX<Eigen::Matrix3d>(n_1, n_2);
    for (int j = 0; j < n_2; j++) {
        for (int i = 0; i < n_1; i++) {
            lambda = -M_PI / 2 + M_PI * (i + 1) / (n_1 + 1);
            phi = -M_PI + 2 * M_PI * j / n_2;
            localUnitVectors(i, j) << -std::sin(lambda) * std::cos(phi), std::sin(phi),  std::cos(lambda) * std::cos(phi),
                                      -std::sin(lambda) * std::cos(phi), std::cos(phi), -std::cos(lambda) * std::sin(phi),
                                       std::cos(lambda)                , 0            ,  std::sin(lambda);
        }
    }
}

void CalcMesh::snapshot(const Eigen::ArrayXXd& u, const Eigen::ArrayXXd& v, const Eigen::ArrayXXd& h, unsigned int snap_number) {
    // Сетка в терминах VTK
    vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();

    // Точки сетки в терминах VTK
    vtkSmartPointer<vtkPoints> dumpPoints = vtkSmartPointer<vtkPoints>::New();

    // Обходим все точки нашей расчётной сетки
    for(unsigned int i = 0; i < n_1; i++) {
        for(unsigned int j = 0; j < n_2; j++) {
            dumpPoints->InsertNextPoint(points(i, j)(0) + scaleFactor * (h(i, j) - hDefault) * localUnitVectors(i, j)(0, 2),
                                        points(i, j)(1) + scaleFactor * (h(i, j) - hDefault) * localUnitVectors(i, j)(1, 2),
                                        points(i, j)(2) + scaleFactor * (h(i, j) - hDefault) * localUnitVectors(i, j)(2, 2));
        }
    }

    // Грузим точки в сетку
    unstructuredGrid->SetPoints(dumpPoints);

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

    for(unsigned int i = 0; i < n_1 - 1; i++) {
        auto quad = vtkSmartPointer<vtkQuad>::New();
        quad->GetPointIds()->SetId(0, n_2 * (i + 1) + n_2 - 1);
        quad->GetPointIds()->SetId(1, n_2 * (i + 1));
        quad->GetPointIds()->SetId(2, n_2 * i);
        quad->GetPointIds()->SetId(3, n_2 * i + n_2 - 1);
        unstructuredGrid->InsertNextCell(quad->GetCellType(), quad->GetPointIds());
    }

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
            // Добавляем значение векторного поля в этой точке
            double _vel[3] = {u(i, j) * localUnitVectors(i, j)(0, 0) + v(i, j) * localUnitVectors(i, j)(0, 1),
                              u(i, j) * localUnitVectors(i, j)(1, 0) + v(i, j) * localUnitVectors(i, j)(1, 1),
                              u(i, j) * localUnitVectors(i, j)(2, 0) + v(i, j) * localUnitVectors(i, j)(2, 1)};
            vel -> InsertNextTuple(_vel);

            // И значение скалярного поля тоже
            h_field -> InsertNextValue(h(i, j));
        }
    }

    // Присоединяем векторное и скалярное поля к точкам
    unstructuredGrid->GetPointData()->AddArray(vel);
    unstructuredGrid->GetPointData()->AddArray(h_field);

    // Создаём снапшот в файле с заданным именем
    std::string fileName = "../out/" + problemName + "_animation/frame-" + std::to_string(snap_number) + ".vts";
    vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
    writer->SetFileName(fileName.c_str());
    writer->SetInputData(unstructuredGrid);
    writer->Write();
}
