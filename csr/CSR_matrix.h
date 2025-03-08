//
// Created by Per Hüpenbecker on 08.03.25.
//
#ifndef CSR_MATRIX_CPP_CSR_MATRIX_H
#define CSR_MATRIX_CPP_CSR_MATRIX_H

#include "../matrix/Matrix.h"

template <typename U>
class CSRMatrix {
    private:
    std::vector<U> data;
    std::vector<size_t> col_indices;
    std::vector<size_t> row_ptrs;
    size_t rows_count;
    size_t cols_count;

    public:
    // Constructor accepting a Matrix object
    CSRMatrix(const Matrix<U>& matrix) {
        rows_count = matrix.get_rows_count();
        cols_count = matrix.get_cols_count();
        size_t row_ptr = 0;
        row_ptrs.push_back(row_ptr);
        auto matrix_data = matrix.get_data();

        std::cout << "Inner Checkpoint" << std::endl;

        for(size_t i = 0; i < rows_count; i++){
            for(size_t j = 0; j < cols_count; j++){
                std::cout << "Inner Checkpoint 2" << std::endl;
                if(matrix_data[i* cols_count + j] != 0){
                    data.push_back(matrix_data[i*cols_count + j]);
                    col_indices.push_back(j);
                    row_ptr++;
                }
            }
            row_ptrs.push_back(row_ptr);
        }
    }

    // display function of CSR Matrix
    void display() {
        for (const auto &i: data) {
            std::cout << i << " ";
        }
    }

    // display method for the full matrix including all zero values

    void display_full(){
        for(size_t i = 0; i < rows_count; i++){

            auto start_index = row_ptrs[i];
            auto end_index = row_ptrs[i+1];

            size_t k = start_index;

            for (size_t j = 0; j < cols_count; j++) {
                if(j == col_indices[k]){
                    std::cout << data[k] << " ";
                    if(k < end_index) k++;

                } else {
                    std::cout << "0 ";
                }
            }
            std::cout<< std::endl;
        }
    }

    void operator*(U scalar){
        for (auto& element: data){
            element *= scalar;
        }
    }

};


#endif //CSR_MATRIX_CPP_CSR_MATRIX_H
