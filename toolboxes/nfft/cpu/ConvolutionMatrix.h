#pragma once

#include "hoArmadillo.h"
#include "hoNDArray.h"
#include "vector_td.h"

#include "ConvolutionKernel.h"

namespace Gadgetron
{
    namespace ConvInternal
    {
        template<class REAL>
        struct ConvolutionMatrix
        {
            ConvolutionMatrix()
            {
                
            }

            ConvolutionMatrix(size_t cols, size_t rows)
              : n_cols(cols),n_rows(rows)
            {
                weights = std::vector<std::vector<REAL>>(n_cols);
                indices = std::vector<std::vector<size_t>>(n_cols);
            }

            std::vector<std::vector<REAL>> weights;
            std::vector<std::vector<size_t>> indices;
            size_t n_cols, n_rows;
        };


        template<class REAL> ConvolutionMatrix<REAL> transpose(
            const ConvolutionMatrix<REAL>& matrix);

        template<class REAL, unsigned int D, template<class, unsigned int> class K>
        ConvolutionMatrix<REAL> make_conv_matrix(
            const hoNDArray<vector_td<REAL, D>> trajectory,
            const vector_td<size_t, D> &matrix_size,
            const ConvolutionKernel<REAL, D, K>& kernel);
    }
}

