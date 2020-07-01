//
// Created by dchansen on 10/12/18.
//

#include <GadgetronTimer.h>
#include <numeric>
#include "hoNFFT_sparseMatrix.h"
#include "KaiserBessel_kernel.h"
#include "vector_td_utilities.h"
#include <boost/range/irange.hpp>

namespace {
    using namespace Gadgetron;





    template<int N>
    struct iteration_counter {
    };


    template<class REAL, unsigned int D>
    void iterate_body(const vector_td<REAL, D> &point,
                      const vector_td<size_t, D> &matrix_size, REAL W, const vector_td<REAL, D> &beta,
                      std::vector<size_t> &indices,
                      std::vector<REAL> &weights, vector_td<REAL, D> &image_point, size_t index, iteration_counter<-1>) {

        indices.push_back(index);
        weights.push_back(KaiserBessel(abs(image_point - point), vector_td<REAL,D>(matrix_size), REAL(1) / W, beta));

    }

    template<class REAL, unsigned int D, int N>
    void iterate_body(const vector_td<REAL, D> &point,
                      const vector_td<size_t, D> &matrix_size, REAL W, const vector_td<REAL, D> &beta,
                      std::vector<size_t> &indices,
                      std::vector<REAL> &weights, vector_td<REAL, D> &image_point, size_t index, iteration_counter<N>) {

        size_t frame_offset = std::accumulate(&matrix_size[0], &matrix_size[N], 1, std::multiplies<size_t>());

        for (int i = std::ceil(point[N] - W * 0.5); i <= std::floor(point[N] + W * 0.5); i++) {
            auto wrapped_i = (i + matrix_size[N]) % matrix_size[N];
            size_t index2 = index + frame_offset * wrapped_i;
            image_point[N] = i;
            iterate_body(point, matrix_size, W, beta, indices, weights, image_point, index2,
                         iteration_counter<N - 1>());
        }
    }

    template<class REAL, unsigned int D>
    std::tuple<std::vector<size_t>, std::vector<REAL>>
    get_indices(const vector_td<REAL, D> &point, const vector_td<size_t, D> &matrix_size,
                REAL W, const vector_td<REAL, D> &beta) {

        std::vector<size_t> indices;
        indices.reserve(size_t(std::pow(std::ceil(W), D)));

        std::vector<REAL> weights;
        weights.reserve(size_t(std::pow(std::ceil(W), D)));

        vector_td<REAL, D> image_point;
        size_t index = 0;
        iterate_body(point, matrix_size, W, beta, indices, weights, image_point, index, iteration_counter<D - 1>());

        return std::make_tuple(std::move(indices), std::move(weights));
    }


}


template<class REAL, unsigned int D>
Gadgetron::NFFT_internal::NFFT_Matrix<REAL>
Gadgetron::NFFT_internal::make_NFFT_matrix(const Gadgetron::hoNDArray<Gadgetron::vector_td<REAL, D>> trajectories,
                                  const Gadgetron::vector_td<size_t, D> &image_dims, REAL W,
                                  const Gadgetron::vector_td<REAL, D> &beta) {

    NFFT_Matrix<REAL> matrix(trajectories.get_number_of_elements(),prod(image_dims));
#pragma omp parallel for 
    for (int i = 0; i < (int)trajectories.get_number_of_elements(); i++) {
        std::tie(matrix.indices[i], matrix.weights[i]) = get_indices(trajectories[i], image_dims, W, beta);
    }
    return matrix;
}



template<class REAL>
Gadgetron::NFFT_internal::NFFT_Matrix<REAL>
Gadgetron::NFFT_internal::transpose(const Gadgetron::NFFT_internal::NFFT_Matrix<REAL> &matrix) {

    auto work_index = boost::irange(size_t(0),matrix.n_cols);

    NFFT_Matrix<REAL> transposed(matrix.n_rows, matrix.n_cols);

    std::vector<int> counts(matrix.n_rows, 0);


    for (size_t i : work_index) {
        auto &rows = matrix.indices[i];
        for (auto &row : rows) {
            counts[row]++;
        }
    }

    for (size_t i = 0; i < counts.size(); i++) {
        transposed.indices[i].reserve(counts[i]);
        transposed.weights[i].reserve(counts[i]);
    }


    for (size_t i : work_index ) {
        auto &rows = matrix.indices[i];
        auto &weights = matrix.weights[i];
        for (size_t n = 0; n < rows.size(); n++) {
            transposed.indices[rows[n]].push_back(i);
            transposed.weights[rows[n]].push_back(weights[n]);
        }
    }

    return transposed;
}



template Gadgetron::NFFT_internal::NFFT_Matrix<float> Gadgetron::NFFT_internal::make_NFFT_matrix<float,1>(
        const Gadgetron::hoNDArray<Gadgetron::vector_td<float, 1>> trajectories,
        const Gadgetron::vector_td<size_t, 1> &image_dims, float W, const Gadgetron::vector_td<float, 1> &beta);

template Gadgetron::NFFT_internal::NFFT_Matrix<float> Gadgetron::NFFT_internal::make_NFFT_matrix<float,2>(
        const Gadgetron::hoNDArray<Gadgetron::vector_td<float, 2>> trajectories,
        const Gadgetron::vector_td<size_t, 2> &image_dims, float W, const Gadgetron::vector_td<float, 2> &beta);

template Gadgetron::NFFT_internal::NFFT_Matrix<float> Gadgetron::NFFT_internal::make_NFFT_matrix<float,3>(
        const Gadgetron::hoNDArray<Gadgetron::vector_td<float, 3>> trajectories,
        const Gadgetron::vector_td<size_t, 3> &image_dims, float W, const Gadgetron::vector_td<float, 3> &beta);


template Gadgetron::NFFT_internal::NFFT_Matrix<double> Gadgetron::NFFT_internal::make_NFFT_matrix<double,1>(
        const Gadgetron::hoNDArray<Gadgetron::vector_td<double, 1>> trajectories,
        const Gadgetron::vector_td<size_t, 1> &image_dims, double W, const Gadgetron::vector_td<double, 1> &beta);

template Gadgetron::NFFT_internal::NFFT_Matrix<double> Gadgetron::NFFT_internal::make_NFFT_matrix<double,2>(
        const Gadgetron::hoNDArray<Gadgetron::vector_td<double, 2>> trajectories,
        const Gadgetron::vector_td<size_t, 2> &image_dims, double W, const Gadgetron::vector_td<double, 2> &beta);

template Gadgetron::NFFT_internal::NFFT_Matrix<double> Gadgetron::NFFT_internal::make_NFFT_matrix<double,3>(
        const Gadgetron::hoNDArray<Gadgetron::vector_td<double, 3>> trajectories,
        const Gadgetron::vector_td<size_t, 3> &image_dims, double W, const Gadgetron::vector_td<double, 3> &beta);

template Gadgetron::NFFT_internal::NFFT_Matrix<float> Gadgetron::NFFT_internal::transpose(
        const Gadgetron::NFFT_internal::NFFT_Matrix<float> &matrix);
template Gadgetron::NFFT_internal::NFFT_Matrix<double> Gadgetron::NFFT_internal::transpose(
        const Gadgetron::NFFT_internal::NFFT_Matrix<double> &matrix);