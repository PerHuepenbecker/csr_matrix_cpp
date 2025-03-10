//
// Created by Per Hüpenbecker on 08.03.25.
//

#include "../csr_matrix_library.h"

// Some example code to demonstrate the usage of the library

int main(){

    // Option 1: Construct a matrix from a vector of values and the number of columns

    std::vector<float> values = {0,2,0,4,0,0,0,0,9};
    Matrix<float> matrix(values, 3);

    // Convert the matrix to a CSR matrix

    CSRMatrix<float> csr_matrix(matrix);

    csr_matrix.display();

    std::cout << std::endl;

    csr_matrix.display_full();

    // Multiply the CSR matrix by a scalar

    csr_matrix * 2;

    std::cout << std::endl;

    csr_matrix.display_full();

    std::vector<float> values2 = {1,0,0,1};
    std::vector<float> values3 = {1,1,1,1};

    Matrix matrix_01 = Matrix(values2, 2);
    Matrix matrix_02 = Matrix(values3, 2);

    std::cout << std::endl;

    CSRMatrix csr_01 = CSRMatrix(matrix_01);

    csr_01.display_full();

    auto result = csr_01 * matrix_02;

    std::cout << std::endl;

    result.display();

    std::vector<double> values4 = {4,1,2,1,3,0,2,0,5};
    Matrix matrix_03 = Matrix(values4, 3);

    std::cout << std::endl;

    matrix_03.display();

    std::cout << "Cholesky Decomposition" << std::endl;

    // Perform Cholesky decomposition on the matrix to egt lower triangular matrix
    auto cholesky = matrix_03.cholesky_decomposition();
    // Get the transpose of the lower triangular matrix to get the upper triangular matrix
    auto cholesky_t = cholesky.T();

    cholesky.display();

    std::cout << std::endl;
    std::cout << "Test" << std::endl;

    // Multiply the lower triangular matrix with its transpose to get the original matrix
    (cholesky*cholesky_t).display();

    // Use the Cholesky decomposition to solve a linear system of equations
    std::cout << "Cholesky Solver" << std::endl;
    std::vector<double> b = {1,2,3};

    auto result_cholesky = matrix_03.solve_cholesky(b);

    for(auto& i: result_cholesky){
        std::cout << i << " ";
    }

    std::cout << std::endl;

    // Calculate the inverse of the matrix using the Cholesky decomposition
    auto inverse = matrix_03.inverse_cholesky();

    inverse.display();

    std::cout << std::endl;
    // Multiply the original matrix with its inverse to get the identity matrix
    std::cout << "Test" << std::endl;

    (matrix_03 * inverse).display();

    return 0;
}