cmake_minimum_required(VERSION 3.27)

project(Shallow_water)

set(SOURCES 
    src/main.cpp
    src/problems.cpp
    src/calc_mesh.cpp
)

add_executable(main ${SOURCES})

find_package(VTK
COMPONENTS
    CommonCore
    CommonDataModel
    IOXML
)

target_link_libraries(main
  PRIVATE
   ${VTK_LIBRARIES}
)

vtk_module_autoinit(
  TARGETS main
  MODULES ${VTK_LIBRARIES}
)