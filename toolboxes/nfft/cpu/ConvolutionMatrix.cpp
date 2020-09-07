//
// Created by dchansen on 10/12/18.
//
#include "ConvolutionMatrix.h"

#include <GadgetronTimer.h>
#include <numeric>
#include "vector_td_utilities.h"
#include <boost/range/irange.hpp>

namespace
{
    using namespace Gadgetron;

    template<int N>
    struct iteration_counter { };

    template<class REAL, unsigned int D, template<class, unsigned int> class K>
    void iterate_body(
        const vector_td<REAL, D> &point,
        const vector_td<size_t, D> &matrix_size,
        std::vector<size_t> &indices,
        std::vector<REAL> &weights,
        vector_td<REAL, D> &image_point,
        size_t index,
        const ConvolutionKernel<REAL, D, K>& kernel,
        iteration_counter<-1>)
    {
        auto delta = abs(image_point - point);
        indices.push_back(index);
        weights.push_back(kernel.get(delta));
    }

    template<class REAL, unsigned int D, template<class, unsigned int> class K, int N>
    void iterate_body(
        const vector_td<REAL, D> &point,
        const vector_td<size_t, D> &matrix_size,
        std::vector<size_t> &indices,
        std::vector<REAL> &weights,
        vector_td<REAL, D> &image_point,
        size_t index,
        const ConvolutionKernel<REAL, D, K>& kernel,
        iteration_counter<N>)
    {
        size_t frame_offset = std::accumulate(&matrix_size[0], &matrix_size[N], 1, std::multiplies<size_t>());

        for (int i = std::ceil(point[N] - kernel.get_radius());
             i <= std::floor(point[N] + kernel.get_radius());
             i++)
        {
            auto wrapped_i = (i + matrix_size[N]) % matrix_size[N];
            size_t index2 = index + frame_offset * wrapped_i;
            image_point[N] = i;
            iterate_body(point, matrix_size, indices, weights, image_point,
                         index2, kernel, iteration_counter<N - 1>());
        }
    }

    template<class REAL, unsigned int D, template<class, unsigned int> class K>
    std::tuple<std::vector<size_t>, std::vector<REAL>> get_indices(
        const vector_td<REAL, D> &point,
        const vector_td<size_t, D> &matrix_size,
        const ConvolutionKernel<REAL, D, K>& kernel)
    {
        std::vector<size_t> indices;
        indices.reserve(size_t(std::pow(std::ceil(kernel.get_width()), D)));

        std::vector<REAL> weights;
        weights.reserve(size_t(std::pow(std::ceil(kernel.get_width()), D)));

        vector_td<REAL, D> image_point;
        size_t index = 0;
        iterate_body(point, matrix_size, indices, weights, image_point, index,
                     kernel, iteration_counter<D - 1>());

        return std::make_tuple(std::move(indices), std::move(weights));
    }
}


template<class REAL, unsigned int D, template<class, unsigned int> class K>
Gadgetron::ConvInternal::ConvolutionMatrix<REAL>
Gadgetron::ConvInternal::make_conv_matrix(
    const Gadgetron::hoNDArray<Gadgetron::vector_td<REAL, D>> trajectory,
    const Gadgetron::vector_td<size_t, D> &matrix_size,
    const ConvolutionKernel<REAL, D, K>& kernel)
{
    ConvolutionMatrix<REAL> matrix(trajectory.get_number_of_elements(),
                                   prod(matrix_size));
    
    #pragma omp parallel for 
    for (int i = 0; i < (int)trajectory.get_number_of_elements(); i++)
    {
        std::tie(matrix.indices[i], matrix.weights[i]) = get_indices(
            trajectory[i], matrix_size, kernel);
    }

    return matrix;
}


template<class REAL>
Gadgetron::ConvInternal::ConvolutionMatrix<REAL>
Gadgetron::ConvInternal::transpose(const Gadgetron::ConvInternal::ConvolutionMatrix<REAL> &matrix) {

    auto work_index = boost::irange(size_t(0),matrix.n_cols);

    ConvolutionMatrix<REAL> transposed(matrix.n_rows, matrix.n_cols);

    std::vector<int> counts(matrix.n_rows, 0);

    for (size_t i : work_index)
    {
        auto &rows = matrix.indices[i];
        for (auto &row : rows)
            counts[row]++;
    }

    for (size_t i = 0; i < counts.size(); i++)
    {
        transposed.indices[i].reserve(counts[i]);
        transposed.weights[i].reserve(counts[i]);
    }

    for (size_t i : work_index )
    {
        auto &rows = matrix.indices[i];
        auto &weights = matrix.weights[i];
        for (size_t n = 0; n < rows.size(); n++)
        {
            transposed.indices[rows[n]].push_back(i);
            transposed.weights[rows[n]].push_back(weights[n]);
        }
    }

    return transposed;
}


template Gadgetron::ConvInternal::ConvolutionMatrix<float>
Gadgetron::ConvInternal::make_conv_matrix<float, 1, Gadgetron::KaiserKernel>(
    const Gadgetron::hoNDArray<Gadgetron::vector_td<float, 1>> trajectory,
    const Gadgetron::vector_td<size_t, 1> &matrix_size,
    const ConvolutionKernel<float, 1, Gadgetron::KaiserKernel>& kernel);

template Gadgetron::ConvInternal::ConvolutionMatrix<float>
Gadgetron::ConvInternal::make_conv_matrix<float, 2, Gadgetron::KaiserKernel>(
    const Gadgetron::hoNDArray<Gadgetron::vector_td<float, 2>> trajectory,
    const Gadgetron::vector_td<size_t, 2> &matrix_size,
    const ConvolutionKernel<float, 2, Gadgetron::KaiserKernel>& kernel);

template Gadgetron::ConvInternal::ConvolutionMatrix<float>
Gadgetron::ConvInternal::make_conv_matrix<float, 3, Gadgetron::KaiserKernel>(
    const Gadgetron::hoNDArray<Gadgetron::vector_td<float, 3>> trajectory,
    const Gadgetron::vector_td<size_t, 3> &matrix_size,
    const ConvolutionKernel<float, 3, Gadgetron::KaiserKernel>& kernel);

template Gadgetron::ConvInternal::ConvolutionMatrix<float>
Gadgetron::ConvInternal::make_conv_matrix<float, 4, Gadgetron::KaiserKernel>(
    const Gadgetron::hoNDArray<Gadgetron::vector_td<float, 4>> trajectory,
    const Gadgetron::vector_td<size_t, 4> &matrix_size,
    const ConvolutionKernel<float, 4, Gadgetron::KaiserKernel>& kernel);

template Gadgetron::ConvInternal::ConvolutionMatrix<double>
Gadgetron::ConvInternal::make_conv_matrix<double, 1, Gadgetron::KaiserKernel>(
    const Gadgetron::hoNDArray<Gadgetron::vector_td<double, 1>> trajectory,
    const Gadgetron::vector_td<size_t, 1> &matrix_size,
    const ConvolutionKernel<double, 1, Gadgetron::KaiserKernel>& kernel);

template Gadgetron::ConvInternal::ConvolutionMatrix<double>
Gadgetron::ConvInternal::make_conv_matrix<double, 2, Gadgetron::KaiserKernel>(
    const Gadgetron::hoNDArray<Gadgetron::vector_td<double, 2>> trajectory,
    const Gadgetron::vector_td<size_t, 2> &matrix_size,
    const ConvolutionKernel<double, 2, Gadgetron::KaiserKernel>& kernel);

template Gadgetron::ConvInternal::ConvolutionMatrix<double>
Gadgetron::ConvInternal::make_conv_matrix<double, 3, Gadgetron::KaiserKernel>(
    const Gadgetron::hoNDArray<Gadgetron::vector_td<double, 3>> trajectory,
    const Gadgetron::vector_td<size_t, 3> &matrix_size,
    const ConvolutionKernel<double, 3, Gadgetron::KaiserKernel>& kernel);

template Gadgetron::ConvInternal::ConvolutionMatrix<double>
Gadgetron::ConvInternal::make_conv_matrix<double, 4, Gadgetron::KaiserKernel>(
    const Gadgetron::hoNDArray<Gadgetron::vector_td<double, 4>> trajectory,
    const Gadgetron::vector_td<size_t, 4> &matrix_size,
    const ConvolutionKernel<double, 4, Gadgetron::KaiserKernel>& kernel);

template Gadgetron::ConvInternal::ConvolutionMatrix<float>
Gadgetron::ConvInternal::make_conv_matrix<float, 1, Gadgetron::JincKernel>(
    const Gadgetron::hoNDArray<Gadgetron::vector_td<float, 1>> trajectory,
    const Gadgetron::vector_td<size_t, 1> &matrix_size,
    const ConvolutionKernel<float, 1, Gadgetron::JincKernel>& kernel);

template Gadgetron::ConvInternal::ConvolutionMatrix<float>
Gadgetron::ConvInternal::make_conv_matrix<float, 2, Gadgetron::JincKernel>(
    const Gadgetron::hoNDArray<Gadgetron::vector_td<float, 2>> trajectory,
    const Gadgetron::vector_td<size_t, 2> &matrix_size,
    const ConvolutionKernel<float, 2, Gadgetron::JincKernel>& kernel);

template Gadgetron::ConvInternal::ConvolutionMatrix<float>
Gadgetron::ConvInternal::make_conv_matrix<float, 3, Gadgetron::JincKernel>(
    const Gadgetron::hoNDArray<Gadgetron::vector_td<float, 3>> trajectory,
    const Gadgetron::vector_td<size_t, 3> &matrix_size,
    const ConvolutionKernel<float, 3, Gadgetron::JincKernel>& kernel);

template Gadgetron::ConvInternal::ConvolutionMatrix<float>
Gadgetron::ConvInternal::make_conv_matrix<float, 4, Gadgetron::JincKernel>(
    const Gadgetron::hoNDArray<Gadgetron::vector_td<float, 4>> trajectory,
    const Gadgetron::vector_td<size_t, 4> &matrix_size,
    const ConvolutionKernel<float, 4, Gadgetron::JincKernel>& kernel);

template Gadgetron::ConvInternal::ConvolutionMatrix<double>
Gadgetron::ConvInternal::make_conv_matrix<double, 1, Gadgetron::JincKernel>(
    const Gadgetron::hoNDArray<Gadgetron::vector_td<double, 1>> trajectory,
    const Gadgetron::vector_td<size_t, 1> &matrix_size,
    const ConvolutionKernel<double, 1, Gadgetron::JincKernel>& kernel);

template Gadgetron::ConvInternal::ConvolutionMatrix<double>
Gadgetron::ConvInternal::make_conv_matrix<double, 2, Gadgetron::JincKernel>(
    const Gadgetron::hoNDArray<Gadgetron::vector_td<double, 2>> trajectory,
    const Gadgetron::vector_td<size_t, 2> &matrix_size,
    const ConvolutionKernel<double, 2, Gadgetron::JincKernel>& kernel);

template Gadgetron::ConvInternal::ConvolutionMatrix<double>
Gadgetron::ConvInternal::make_conv_matrix<double, 3, Gadgetron::JincKernel>(
    const Gadgetron::hoNDArray<Gadgetron::vector_td<double, 3>> trajectory,
    const Gadgetron::vector_td<size_t, 3> &matrix_size,
    const ConvolutionKernel<double, 3, Gadgetron::JincKernel>& kernel);

template Gadgetron::ConvInternal::ConvolutionMatrix<double>
Gadgetron::ConvInternal::make_conv_matrix<double, 4, Gadgetron::JincKernel>(
    const Gadgetron::hoNDArray<Gadgetron::vector_td<double, 4>> trajectory,
    const Gadgetron::vector_td<size_t, 4> &matrix_size,
    const ConvolutionKernel<double, 4, Gadgetron::JincKernel>& kernel);

template Gadgetron::ConvInternal::ConvolutionMatrix<float>
Gadgetron::ConvInternal::transpose(
    const Gadgetron::ConvInternal::ConvolutionMatrix<float> &matrix);

template Gadgetron::ConvInternal::ConvolutionMatrix<double>
Gadgetron::ConvInternal::transpose(
        const Gadgetron::ConvInternal::ConvolutionMatrix<double> &matrix);

