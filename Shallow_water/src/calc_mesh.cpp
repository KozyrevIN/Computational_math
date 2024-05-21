#ifndef calc_mesh
    #define calc_mesh
    #include "../include/calc_mesh.h"
#endif

CalcMesh::CalcMesh(FlatProblem problem, unsigned n_x, unsigned n_y) : n_1(n_x), n_2(n_y) {
    // Сетка в терминах VTK
    unstructuredGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
    // Точки сетки в терминах VTK
    vtkSmartPointer<vtkPoints> dumpPoints = vtkSmartPointer<vtkPoints>::New();

    // Обходим все точки нашей расчётной сетки
    for(unsigned int i = 0; i < n_1; i++) {
        for(unsigned int j = 0; j < n_2; j++) {
            dumpPoints->InsertNextPoint((problem.l_x * i) / (n_1 - 1), (problem.l_y * j) / (n_2 - 1), 0);
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

    problemName = problem.name;
}

void CalcMesh::flatProject(Eigen::ArrayXXd& u, Eigen::ArrayXXd& v, Eigen::ArrayXXd& h)
{
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
            double _vel[3] = {u(i, j), v(i, j), 0};
            vel -> InsertNextTuple(_vel);

            // И значение скалярного поля тоже
            h_field -> InsertNextValue(h(i, j));
        }
    }

    // Присоединяем векторное и скалярное поля к точкам
    unstructuredGrid->GetPointData()->AddArray(vel);
    unstructuredGrid->GetPointData()->AddArray(h_field);
}

void CalcMesh::snapshot(unsigned int snap_number) {
    // Создаём снапшот в файле с заданным именем
    std::string fileName = "../out/" + problemName + "_animation/frame-" + std::to_string(snap_number) + ".vts";
    vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
    writer->SetFileName(fileName.c_str());
    writer->SetInputData(unstructuredGrid);
    writer->Write();
}
