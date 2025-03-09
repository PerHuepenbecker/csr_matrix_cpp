## Basic C++ library for matrix operations  

This tiny library contains a simple implementation of basic matrix operations, including the Compressed Sparse Row  (CSR) approach for efficient operations on sparse matrices.
I wrote this library as a tiny refresher on linear algebra and specifically on matrices and as a preparation for my bachelors thesis. 

### Features

- templated code using concept constraints
- matrix addition, subtraction, multiplication implemented via operator overloading
- transpose, inverse calculation
- CSR matrix representation
- basic CSR matrix operations (addition, subtraction, multiplication)
- basic matrix decomposition (Cholesky)
- basic testing using googletest

### Wishlist

- more matrix decompositions (LU, QR)
- matrix solvers (Gauss-Seidel, SOR)
- performance improvements

### Usage

Include the header file `csr_matrix_library.h` in your project and use the `Matrix` class to create matrices. The `CSRMatrix` class is used to create matrices in the CSR format. 

### Examples

#### Basic matrix operations
```cpp
#include "csr_matrix_library.h"

std::vector<double> data = {1, 2, 3, 4, 5, 6, 7, 8, 9};

// Deep copy of the data vector evaded here but also possible
// Constructs a 3x3 matrix from the data vector in row-major order
Matrix<double> matrix(std::move(data), 3);

// Display the matrix
matrix_t.display();

// Transpose the matrix
Matrix<double> matrix_t = matrix.T()
        
// Display the transposed matrix
matrix_t.display();

// Perform scalar multiplication in place or return a new matrix
matrix_t *= 2.0; // In place
Matrix<double> matrix_scaled = matrix_t * 2.0; // New matrix

// Use a different constructor to create a 1x3 matrix
Matrix<double> matrix_2({{2}{2}{2}});

// Perform matrix multiplication
Matrix<double> matrix_result = matrix * matrix_2;
```

#### CSR matrix operations
```cpp
#include "csr_matrix_library.h"

std::vector<double> data = {1, 0, 0, 4, 0, 6, 0, 8, 0};
Matrix<double> matrix(std::move(data), 3);

// Convert the matrix to CSR format

CSRMatrix<double> csr_matrix(matrix);

// Display the CSR matrix non-zero values
csr_matrix.display();

// Display the CSR matrix full representation
csr_matrix.display_full();

// Perform matrix multiplication of a CSR matrix with a dense matrix
Matrix<double> matrix_2({{2},{2},{2}});

// Perform matrix multiplication only iteration over the non-zero values of the CSR matrix
Matrix<double> matrix_result = csr_matrix * matrix_2;
```
