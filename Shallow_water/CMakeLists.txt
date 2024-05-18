cmake_minimum_required(VERSION 3.27)

project(Shallow_water)

set(SOURCES 
    src/main.cpp
    src/solver.cpp
    src/calc_mesh.cpp
    src/progress_bar.cpp
)

add_executable(main ${SOURCES})

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

target_link_libraries(main ${VTK_LIBRARIES})