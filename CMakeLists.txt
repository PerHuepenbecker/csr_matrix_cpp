cmake_minimum_required(VERSION 3.30)
project(csr_matrix_cpp)

set(CMAKE_CXX_STANDARD 23)

add_library(csr_matrix_cpp STATIC library.cpp
        matrix/Matrix.h
        )

add_executable(example example/main.cpp)
target_link_libraries(example csr_matrix_cpp)



