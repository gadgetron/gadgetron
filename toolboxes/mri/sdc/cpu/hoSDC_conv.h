/**
 * \file hoSDC_conv.h
 * \brief Convolution for sampling density compensation (host).
 */

#pragma once

#include "hoSDC_kernel.h"

#include "hoArmadillo.h"
#include "hoNDArray.h"
#include "vector_td.h"


namespace Gadgetron
{
    namespace SDC_internal
    {
        template<class REAL>
        struct hoConvMatrix
        {
            hoConvMatrix()
            { }

            hoConvMatrix(size_t cols, size_t rows) : n_cols(cols),n_rows(rows)
            {
                weights = std::vector<std::vector<REAL>>(n_cols);
                indices = std::vector<std::vector<size_t>>(n_cols);
            }

            std::vector<std::vector<REAL>> weights;
            std::vector<std::vector<size_t>> indices;
            size_t n_cols, n_rows;
        };

        template<class REAL>
        hoConvMatrix<REAL> transpose(const hoConvMatrix<REAL>& matrix);

        template<class REAL, unsigned int D>
        hoConvMatrix<REAL> make_conv_matrix(const hoNDArray<vector_td<REAL, D>> trajectories, const hoSDC_kernel<REAL, D>& kernel);

    }   // namespace SDC_internal
}   // namespace Gadgetron