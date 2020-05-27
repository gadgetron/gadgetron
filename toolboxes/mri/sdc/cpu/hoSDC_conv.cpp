/**
 * \file hoSDC_conv.cpp
 * \brief Convolution for sampling density compensation (host).
 */

#include "hoSDC_conv.h"

#include "vector_td_utilities.h"

#include <numeric>
#include <boost/range/irange.hpp>

namespace
{
    using namespace Gadgetron;

    template<int N>
    struct iteration_counter
    { };

    template<class REAL, unsigned int D>
    void iterate_body(const vector_td<REAL, D> &point, const hoSDC_kernel<REAL, D>& kernel,
                      std::vector<size_t> &indices, std::vector<REAL> &weights,
                      vector_td<REAL, D> &grid_point, size_t index, iteration_counter<-1>)
    {
        indices.push_back(index);
        vector_td<REAL, D> lut_point = abs(grid_point - point);
        weights.push_back(kernel.lookup(lut_point));
    }

    template<class REAL, unsigned int D, int N>
    void iterate_body(const vector_td<REAL, D> &point, const hoSDC_kernel<REAL, D>& kernel,
                      std::vector<size_t> &indices, std::vector<REAL> &weights,
                      vector_td<REAL, D> &grid_point, size_t index, iteration_counter<N>)
    {
        size_t frame_offset = std::accumulate(&kernel.get_grid_size()[0], &kernel.get_grid_size()[N], 1, std::multiplies<size_t>());

        for (int i = std::ceil(point[N] - kernel.get_rfp()); i <= std::floor(point[N] + kernel.get_rfp()); i++)
        {
            auto wrapped_i = (i + kernel.get_grid_size()[N]) % kernel.get_grid_size()[N];
            size_t index2 = index + frame_offset * wrapped_i;
            grid_point[N] = i;
            iterate_body(point, kernel, indices, weights, grid_point, index2, iteration_counter<N - 1>());
        }
    }

    template<class REAL, unsigned int D>
    std::tuple<std::vector<size_t>, std::vector<REAL>> get_indices(const vector_td<REAL, D> &point, const hoSDC_kernel<REAL, D>& kernel)
    {
        // kernel width
        REAL W = 2.0 * kernel.get_rfp();

        std::vector<size_t> indices;
        indices.reserve(size_t(std::pow(std::ceil(W), D)));

        std::vector<REAL> weights;
        weights.reserve(size_t(std::pow(std::ceil(W), D)));

        vector_td<REAL, D> grid_point;
        size_t index = 0;
        iterate_body(point, kernel, indices, weights, grid_point, index, iteration_counter<D - 1>());

        return std::make_tuple(std::move(indices), std::move(weights));
    }
}

namespace Gadgetron
{
    namespace SDC_internal
    {

        template<class REAL, unsigned int D>
        hoConvMatrix<REAL> make_conv_matrix(const hoNDArray<vector_td<REAL, D>> traj, const hoSDC_kernel<REAL, D>& kernel)
        {
            hoConvMatrix<REAL> matrix(traj.get_number_of_elements(), prod(kernel.get_grid_size()));

            #pragma omp parallel for 
            for (int i = 0; i < (int)traj.get_number_of_elements(); i++)
            {
                std::tie(matrix.indices[i], matrix.weights[i]) = get_indices(traj[i], kernel);
            }

            return matrix;
        }


        template<class REAL>
        hoConvMatrix<REAL> transpose(const hoConvMatrix<REAL> &matrix)
        {
            auto work_index = boost::irange(size_t(0), matrix.n_cols);

            hoConvMatrix<REAL> transposed(matrix.n_rows, matrix.n_cols);

            std::vector<int> counts(matrix.n_rows, 0);

            for (size_t i : work_index)
            {
                auto &rows = matrix.indices[i];
                for (auto &row : rows)
                {
                    counts[row]++;
                }
            }

            for (size_t i = 0; i < counts.size(); i++)
            {
                transposed.indices[i].reserve(counts[i]);
                transposed.weights[i].reserve(counts[i]);
            }

            for (size_t i : work_index)
            {
                // GDEBUG("i: %d\n", i);
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

    }   // namespace SDC_internal
}   // namespace Gadgetron


template Gadgetron::SDC_internal::hoConvMatrix<float> Gadgetron::SDC_internal::make_conv_matrix(const hoNDArray<vector_td<float, 1>> traj, const hoSDC_kernel<float, 1>& kernel);
template Gadgetron::SDC_internal::hoConvMatrix<float> Gadgetron::SDC_internal::make_conv_matrix(const hoNDArray<vector_td<float, 2>> traj, const hoSDC_kernel<float, 2>& kernel);
template Gadgetron::SDC_internal::hoConvMatrix<float> Gadgetron::SDC_internal::make_conv_matrix(const hoNDArray<vector_td<float, 3>> traj, const hoSDC_kernel<float, 3>& kernel);

template Gadgetron::SDC_internal::hoConvMatrix<double> Gadgetron::SDC_internal::make_conv_matrix(const hoNDArray<vector_td<double, 1>> traj, const hoSDC_kernel<double, 1>& kernel);
template Gadgetron::SDC_internal::hoConvMatrix<double> Gadgetron::SDC_internal::make_conv_matrix(const hoNDArray<vector_td<double, 2>> traj, const hoSDC_kernel<double, 2>& kernel);
template Gadgetron::SDC_internal::hoConvMatrix<double> Gadgetron::SDC_internal::make_conv_matrix(const hoNDArray<vector_td<double, 3>> traj, const hoSDC_kernel<double, 3>& kernel);

template Gadgetron::SDC_internal::hoConvMatrix<float> Gadgetron::SDC_internal::transpose(const Gadgetron::SDC_internal::hoConvMatrix<float>& matrix);
template Gadgetron::SDC_internal::hoConvMatrix<double> Gadgetron::SDC_internal::transpose(const Gadgetron::SDC_internal::hoConvMatrix<double>& matrix);