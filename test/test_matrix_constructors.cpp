//
// Created by Per Hüpenbecker on 09.03.25.
//


#include <gtest/gtest.h>
#include "../csr_matrix_library.h"

// Test case for matrix construction
TEST(MatrixTest, Construction) {
Matrix<double> mat({{1.0, 2.0}, {3.0, 4.0}});
EXPECT_EQ(mat.get_rows_count(), 2);
EXPECT_EQ(mat.get_cols_count(), 2);
EXPECT_EQ(mat.get_data()[0], 1.0);
EXPECT_EQ(mat.get_data()[3], 4.0);

std::vector<double> values = {1, 2, 3, 4, 5, 6, 7, 8, 9};
Matrix<double> matrix(values, 3);
EXPECT_EQ(matrix.get_rows_count(), 3);
EXPECT_EQ(matrix.get_cols_count(), 3);
EXPECT_EQ(matrix.get_data()[0], 1.0);
EXPECT_EQ(matrix.get_data()[8], 9.0);
}

// Test case for addition
TEST(MatrixTest, Addition) {
Matrix<double> mat1({{1.0, 2.0}, {3.0, 4.0}});
Matrix<double> mat2({{5.0, 6.0}, {7.0, 8.0}});
Matrix<double> result = mat1 + mat2;

EXPECT_EQ(result.get_data(), std::vector<double>({6.0, 8.0, 10.0, 12.0}));
}

// Test case for subtraction
TEST(MatrixTest, Subtraction) {
Matrix<double> mat1({{5.0, 6.0}, {7.0, 8.0}});
Matrix<double> mat2({{1.0, 2.0}, {3.0, 4.0}});
Matrix<double> result = mat1 - mat2;

EXPECT_EQ(result.get_data(), std::vector<double>({4.0, 4.0, 4.0, 4.0}));
}

// Test case for scalar multiplication
TEST(MatrixTest, ScalarMultiplication) {
Matrix<double> mat({{1.0, 2.0}, {3.0, 4.0}});
mat *= 2.0;

EXPECT_EQ(mat.get_data(), std::vector<double>({2.0, 4.0, 6.0, 8.0}));
}

// Test case for matrix multiplication
TEST(MatrixTest, MatrixMultiplication) {
Matrix<double> mat1({{1.0, 2.0}, {3.0, 4.0}});
Matrix<double> mat2({{2.0, 0.0}, {1.0, 2.0}});
Matrix<double> result = mat1 * mat2;

EXPECT_EQ(result.get_data(), std::vector<double>({4.0, 4.0, 10.0, 8.0}));
}

// Test case for transposition
TEST(MatrixTest, Transpose) {
Matrix<double> mat({{1.0, 2.0}, {3.0, 4.0}});
Matrix<double> result = mat.T();

EXPECT_EQ(result.get_data(), std::vector<double>({1.0, 3.0, 2.0, 4.0}));
}

// Test case for Cholesky Decomposition
TEST(MatrixTest, CholeskyDecomposition) {
Matrix<double> mat({{4.0, 12.0, -16.0}, {12.0, 37.0, -43.0}, {-16.0, -43.0, 98.0}});
Matrix<double> result = mat.cholesky_decomposition();

std::vector<double> expected = {
        2.0, 0.0, 0.0,
        6.0, 1.0, 0.0,
        -8.0, 5.0, 3.0
};

EXPECT_EQ(result.get_data(), expected);
}

// Test case for inverse using Cholesky Decomposition
TEST(MatrixTest, InverseCholesky) {
Matrix<double> mat({{4.0, 12.0, -16.0}, {12.0, 37.0, -43.0}, {-16.0, -43.0, 98.0}});
Matrix<double> result = mat.inverse_cholesky();

std::vector<double> expected = {
        49.361, -13.556, 2.111,
        -13.556, 3.778, -0.556,
        2.111, -0.556, 0.111
};

EXPECT_NEAR(result.get_data()[0], expected[0], 1e-3);
EXPECT_NEAR(result.get_data()[4], expected[4], 1e-3);
EXPECT_NEAR(result.get_data()[8], expected[8], 1e-3);
}

// Test case for regularization
TEST(MatrixTest, Regularization) {
Matrix<double> mat({{0.0, 0.0}, {0.0, 0.0}});
mat.regularize();
Matrix<double> mat_02({{1.0, 1.0}, {1.0, 1.0}});
mat_02.regularize();

EXPECT_EQ(mat.get_data()[0], 0.0);
EXPECT_EQ(mat.get_data()[3], 0.0);
EXPECT_GT(mat_02.get_data()[0], 1.0);
EXPECT_GT(mat_02.get_data()[3], 1.0);

}

// Test case for edge cases
TEST(MatrixTest, EdgeCases) {
// Empty matrix
EXPECT_THROW(Matrix<double> mat({}, 1), std::invalid_argument);

// Incompatible dimensions for addition
Matrix<double> mat1({{1.0, 2.0}});
Matrix<double> mat2({{1.0, 2.0}, {3.0, 4.0}});
EXPECT_THROW(mat1 + mat2, std::invalid_argument);

// Non-square matrix for Cholesky
Matrix<double> non_square({{1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}});
EXPECT_THROW(non_square.cholesky_decomposition(), std::invalid_argument);

// Non-positive-definite matrix for Cholesky
Matrix<double> non_pos_def({{1.0, 2.0}, {2.0, 1.0}});
EXPECT_THROW(non_pos_def.cholesky_decomposition(), std::invalid_argument);
}

// Main function for running all tests
int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
