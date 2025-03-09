cmake_minimum_required(VERSION 3.30)
project(csr_matrix_cpp)

set(CMAKE_CXX_STANDARD 23)

# Automatically detect Homebrew on macOS
if(APPLE)
    message(STATUS "Configuring for macOS with Homebrew")

    list(APPEND CMAKE_PREFIX_PATH /opt/homebrew /usr/local /usr/local/opt)

    find_package(PkgConfig REQUIRED)

    pkg_check_modules(GTEST REQUIRED gtest)

    include_directories(${GTEST_INCLUDE_DIRS})
    link_directories(${GTEST_LIBRARY_DIRS})
endif()


add_library(csr_matrix_cpp STATIC csr_matrix_library.cpp
        matrix/Matrix.h
        csr/CSR_matrix.h
        test/test_matrix_constructors.cpp
)

add_executable(test_matrix_constructors test/test_matrix_constructors.cpp)
target_link_libraries(test_matrix_constructors csr_matrix_cpp ${GTEST_LIBRARIES} pthread)

add_executable(example example/main.cpp)
target_link_libraries(example csr_matrix_cpp)



