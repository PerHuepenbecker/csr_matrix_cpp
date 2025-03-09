//
// Created by Per Hüpenbecker on 07.03.25.
//
#ifndef CSR_MATRIX_CPP_MATRIX_H
#define CSR_MATRIX_CPP_MATRIX_H

#include <vector>
#include <memory>
#include <iostream>
#include <algorithm>

template <typename U>
class Matrix {
private:
    std::vector<U> data;
    size_t rows_count;
    size_t cols_count;

public:

    // Generic Constructor accepting any container that provides begin and end iterators
    template<typename Iterable>
    Matrix(const Iterable& values, int cols_count) :
            data(values.begin(), values.end()),
            rows_count(data.size() / cols_count),
            cols_count(cols_count)
    {
        if (data.size() % cols_count != 0) {
            throw std::invalid_argument("The number of elements must be a multiple of the row length.");
        }
    }

    // Constructor accepting a vector of vectors of type U
    Matrix(const std::vector<std::vector<U>>& values) {
        cols_count = values.at(0).size();
        for (const auto& row: values){
            if (row.size() != cols_count) {
                data.clear();
                throw std::invalid_argument("All rows must have the same length.");
            }
            data.insert(data.end(), row.begin(), row.end());
        }
        rows_count = values.size();
    }

    // Constructor accepting a right hand value to avoid unnecessary copies

    Matrix(std::vector<U>&& values, int cols_count)
        :data(std::move(values)), rows_count(values.size()/cols_count), cols_count(cols_count) {

        if (data.size() % cols_count != 0) {
            throw std::invalid_argument("The number of elements must be a multiple of the row length.");
        }
    }

    // Constructor accepting a vector of type U and the number of columns
    Matrix(const std::vector<U>& values, int cols_count) :
            data(values),
            rows_count(data.size() / cols_count),
            cols_count(cols_count)
    {
        if (values.size() % cols_count != 0) {
            throw std::invalid_argument("The number of elements must be a multiple of the row length.");
        }
    }

    // Constructor accepting a smart pointer to an array of type U, the number of expected values
    // and the number of columns. Currently not as efficient as possible since it performs a deep copy
    // of the values. Future updates will address this issue and provide a more efficient solution using
    // shared pointers and different data holders.

    Matrix(const std::shared_ptr<U>& values, size_t values_count, int num_cols) :
            data(values.get(), values.get() + values_count),
            rows_count(values_count / num_cols),
            cols_count(num_cols)
    {
        if (values_count % num_cols != 0) {
            throw std::invalid_argument("The number of elements must be a multiple of the row length.");
        }
    }

    // Overloaded operator for in place addition of two matrices. The left hand matrice will hold the combined value.

    void operator+=(const Matrix& rhs) {
        if(rows_count != rhs.rows_count || cols_count != rhs.cols_count) {
            throw std::invalid_argument("Matrices must have the same dimensions.");
        }
        std::transform(data.begin(), data.end(), rhs.data.begin(), data.begin(), std::plus<U>());
    }

    // Overloaded operator for in place subtraction of two matrices. The left hand matrice will hold the combined value.

    void operator-=(const Matrix& rhs){
        if(rows_count != rhs.rows_count || cols_count != rhs.cols_count){
            throw std::invalid_argument("Matrices must have the same dimensions.");
        }

        std::transform(data.begin(), data.end(), rhs.data.begin(), data.begin(), std::minus<U>());
    }

    // Overloaded operator for scalar multiplication of a matrix. The matrix will be modified in place.

    void operator*=(const U scalar){
        std::transform(data.begin(), data.end(), data.begin(), [scalar](U value) {return value*scalar;});
    }

    void operator*=(const Matrix& rhs){
        if(cols_count != rhs.rows_count){
            throw std::invalid_argument("Incompatible matrix dimensions for multiplication");
        }

        // Temporary buffer to store the result of the matrix multiplication - necessary
        // to avoid modification of the original data buffer during the multiplication which
        // would corrupt the result.

        std::vector<U> tmp_result(rows_count * rhs.cols_count, 0);

        for (size_t i = 0; i < rows_count; i++){
            for(size_t k = 0; k < cols_count; k++){
                U value = data[i * cols_count + k];

                for(size_t j = 0; j < rhs.cols_count; j++){
                    tmp_result[i*rhs.cols_count + j] += value * rhs.data[k * rhs.cols_count + j];
                }
            }
        }

        // move the result back to the original data buffer to avoid unnecessary copies
        data = std::move(tmp_result);
    }

    // Overloaded operator for addition of two matrices. The matrices must have the same dimensions.
    // The result is a new Matrix object.

    friend Matrix operator+(const Matrix& lhs, const Matrix& rhs) {
        if(lhs.rows_count != rhs.rows_count || lhs.cols_count != rhs.cols_count){
            throw std::invalid_argument("Matrices must have the same dimensions.");
        }
        std::vector<U> result(lhs.data.size());
        std::transform(lhs.data.begin(), lhs.data.end(), rhs.data.begin(), result.begin(), std::plus<U>());

        return Matrix(result, lhs.cols_count);
    }

    friend Matrix operator-(const Matrix&lhs, const Matrix&rhs) {
        if(lhs.rows_count != rhs.rows_count || lhs.cols_count != rhs.cols_count){
            throw std::invalid_argument("Matrices must have the same dimensions.");
        }
        std::vector<U> result(lhs.data.size());
        std::transform(lhs.data.begin(), lhs.data.end(), rhs.data.begin(), result.begin(), std::minus<U>());

        return Matrix(result, lhs.cols_count);
    }

    // Overloaded operator for matrix multiplication of two matrices. Uses the ikj order matrix multiplication algorithm
    // to achieve better cache performance. The result is a new Matrix object.

    friend Matrix operator*(const Matrix&lhs, const Matrix&rhs) {
        if (lhs.cols_count != rhs.rows_count) {
            throw std::invalid_argument(
                    "The number of columns in the first matrix must be equal to the number of rows in the second matrix.");
        }

        std::vector<U> result(lhs.rows_count * rhs.cols_count);
        for (size_t i = 0; i < lhs.rows_count; i++) {

            for (size_t k = 0; k < lhs.cols_count; k++) {
                U lhs_value = lhs.data[i * lhs.cols_count + k];
                for (size_t j = 0; j < rhs.cols_count; j++) {
                    result[i * rhs.cols_count + j] += lhs_value * rhs.data[k * rhs.cols_count + j];
                }
            }
        }
        return Matrix(std::move(result), rhs.cols_count);
    }

    // Scalar multiplication of a matrix. The result is a new Matrix object.

    friend Matrix operator*(Matrix& lhs, const U scalar) {
        std::vector<U> result(lhs.data.size());
        std::transform(lhs.data.begin(), lhs.data.end(),result.begin(), [scalar](U value) {return value*scalar;});
    }

    // Method to transpose the matrix. The method uses a block transpose algorithm to
    // achieve better cache performance. The block size is set to 64 and the method returns a new Matrix object.

    Matrix T(){
        const size_t blockSize = 64;
        std::vector<U> result(data.size());

        for(size_t i=0; i < rows_count; i += blockSize){
            for(size_t j = 0; j < cols_count; j += blockSize){
                for(size_t ii = i; ii < std::min(i + blockSize, rows_count); ++ii){
                    for(size_t jj = j; jj < std::min(j + blockSize, cols_count); ++jj) {
                        result[jj * rows_count +  ii] = data[ii * cols_count + jj];
                    }
                }
            }
        }

        for(size_t i = 0; i < rows_count; i++){
            for(size_t j = 0; j < cols_count; j++){
                result[j * rows_count + i] = data[i * cols_count + j];
            }
        }
        return Matrix(result, rows_count);
    }

    std::vector<U> get_data() const {
        return data;
    }

    size_t get_rows_count() const {
        return rows_count;
    }
    size_t get_cols_count() const {
        return cols_count;
    }

    void display() {
        for(size_t i = 0; i < rows_count; i++){
            for(size_t j = 0; j < cols_count; j++){
                std::cout << data[i*cols_count+j] << " ";
            }
            std::cout << std::endl;
        }
    }
};

#endif //CSR_MATRIX_CPP_MATRIX_H