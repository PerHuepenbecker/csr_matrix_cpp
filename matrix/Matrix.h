//
// Created by Per Hüpenbecker on 07.03.25.
//

#ifndef CSR_MATRIX_CPP_MATRIX_H
#define CSR_MATRIX_CPP_MATRIX_H

#include <vector>
#include <memory>
#include <iostream>
#include <algorithm>

template <typename T>
class Matrix {
private:
    std::vector<T> data;
    size_t rows_count;
    size_t cols_count;

public:

    // Generic Constructor accepting any container that provides begin and end iterators

    template<typename Iterable>
Matrix(const Iterable& values, int row_length)  {
        data = std::vector<T>(values.begin(), values.end());
        if (data.size() % row_length != 0) {
            throw std::invalid_argument("The number of elements must be a multiple of the row length.");
        }
        rows_count = data.size() / row_length;
        cols_count = row_length;
    }

    // Constructor accepting a vector of vectors of type T

    Matrix(const std::vector<std::vector<T>>& values) {
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

    // Constructor accepting a vector of type T and the number of columns

    Matrix(const std::vector<T>& values, int row_length) {
        if (values.size() % row_length != 0) {
            throw std::invalid_argument("The number of elements must be a multiple of the row length.");
        }
        data = values;
        rows_count = data.size() / row_length;
        cols_count = row_length;
    }

    // Constructor accepting a smart pointer to an array of type T, the number of expected values
    // and the number of columns. Currently not as efficient as possible since it performs a deep copy
    // of the values. Future updates will address this issue and provide a more efficient solution using
    // shared pointers and different data holders.

    Matrix(const std::shared_ptr<T>& values, size_t values_count, int num_cols) {
        if (values_count % num_cols != 0) {
            throw std::invalid_argument("The number of elements must be a multiple of the row length.");
        }
        data = std::vector<T>(values.get(), values.get() + values_count);
        rows_count = values_count / num_cols;
        cols_count = num_cols;
    }

    friend Matrix operator+(const Matrix& lhs, const Matrix& rhs){
        if(lhs.rows_count != rhs.rows_count || lhs.cols_count != rhs.cols_count){
            throw std::invalid_argument("Matrices must have the same dimensions.");
        }
        std::vector<T> result(lhs.data.size());
        std::transform(lhs.data.begin(), lhs.data.end(), rhs.data.begin(), result.begin(), std::plus<T>());

        return Matrix(result, lhs.cols_count);
    }

    friend Matrix operator-(const Matrix&lhs, const Matrix&rhs){
        if(lhs.rows_count != rhs.rows_count || lhs.cols_count != rhs.cols_count){
            throw std::invalid_argument("Matrices must have the same dimensions.");
        }
        std::vector<T> result(lhs.data.size());
        std::transform(lhs.data.begin(), lhs.data.end(), rhs.data.begin(), result.begin(), std::minus<T>());

        return Matrix(result, lhs.cols_count);
    }

    friend Matrix operator*(const Matrix&lhs, const Matrix&rhs){
        if(lhs.cols_count != rhs.rows_count){
            throw std::invalid_argument("The number of columns in the first matrix must be equal to the number of rows in the second matrix.");
        }
        std::vector<T> result(lhs.rows_count * rhs.cols_count);
        for(size_t i = 0; i < lhs.rows_count; i++){
            for(size_t j = 0; j < rhs.cols_count; j++){
                T sum = 0;
                for(size_t k = 0; k < lhs.cols_count; k++){
                    sum += lhs.data[i * lhs.cols_count + k] * rhs.data[k * rhs.cols_count + j];
                }
                result[i * rhs.cols_count + j] = sum;
            }
        }
        return Matrix(result, rhs.cols_count);
    }

    // method for displaying the matrix

    void display(){
        for(size_t i = 0; i < rows_count; i++){
            for(size_t j = 0; j < cols_count; j++){
                std::cout << data[i*cols_count+j] << " ";
            }
            std::cout << std::endl;
        }
    }

};


#endif //CSR_MATRIX_CPP_MATRIX_H
