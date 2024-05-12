cmake_minimum_required(VERSION 3.27)

project(Lab_3)

set(SOURCES 
    src/main.cpp
    src/problem.cpp
    src/problem_solve.cpp
    src/progress_bar.cpp
)

add_executable(main ${SOURCES})