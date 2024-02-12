cmake_minimum_required(VERSION 3.27.4)

project(Lab_1)

add_executable(${PROJECT_NAME} src/main.cpp src/solver.cpp src/shoot.cpp src/thomas.cpp)

target_link_libraries(${PROJECT_NAME} ${EIGEN_LIBRARY})