#pragma once

#include "hoArmadillo.h"
#include "hoNDArray.h"
#include "vector_td.h"

namespace Gadgetron {
    namespace NFFT_internal {

        template<class REAL> struct NFFT_Matrix {

            NFFT_Matrix(size_t cols, size_t rows) : n_cols(cols),n_rows(rows) {
                weights = std::vector<std::vector<REAL>>(n_cols);
                indices = std::vector<std::vector<size_t>>(n_cols);
            }
            NFFT_Matrix(){}
            std::vector<std::vector<REAL>> weights;
            std::vector<std::vector<size_t>> indices;
            size_t n_cols, n_rows;
        };


        template<class REAL> NFFT_Matrix<REAL> transpose(const NFFT_Matrix<REAL>& matrix);

        template<class REAL, unsigned int D>
        NFFT_Matrix<REAL>
        make_NFFT_matrix(const hoNDArray<vector_td<REAL, D>> trajectories, const vector_td<size_t, D> &image_dims,
                         REAL W, const vector_td<REAL, D> &beta);
    }
}


