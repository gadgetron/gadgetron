/**
 * \file hoSDC.cpp
 * \brief Sampling density compensation (host specialization).
 */

#include "hoSDC.h"
#include "hoSDC_kernel.h"

#include "hoNDArray_elemwise.h"
#include "hoNDArray_reductions.h"
#include "NDArray_utils.h" // NDArrayViewRange
#include "vector_td_operators.h"

namespace Gadgetron
{
    namespace
    {
        template<class REAL>
        void matrix_vector_multiply(const Gadgetron::SDC_internal::hoConvMatrix<REAL>& matrix, const REAL* vector, REAL* result)
        {
            for (size_t i = 0; i < matrix.n_cols; i++)
            {
                auto &row_indices = matrix.indices[i];
                auto &weights = matrix.weights[i];

                #ifndef WIN32
                    #pragma omp simd
                #endif // WIN32
                for (size_t n = 0; n < row_indices.size(); n++)
                {
                    result[i] += vector[row_indices[n]] * weights[n];
                }
            }
        }
    }   // namespace

    template<class REAL, unsigned int D>
    hoSDC_impl<REAL, D>::hoSDC_impl(const vector_td<size_t, D>& matrix_size, REAL os_factor, size_t num_iterations)
      : SDC_impl<hoNDArray, REAL, D>(matrix_size, os_factor, num_iterations)
    {
        this->kernel_ = hoSDC_kernel<REAL, D>(this->matrix_size_, this->grid_size_);
    }

    template<class REAL, unsigned int D>
    void hoSDC_impl<REAL, D>::preprocess(const hoNDArray<vector_td<REAL, D>>& traj)
    {
        // call base class
        SDC_impl<hoNDArray, REAL, D>::preprocess(traj);

        // scale trajectories
        vector_td<REAL, D> grid_size = vector_td<REAL, D>(this->grid_size_);
        vector_td<REAL, D> win_length = (grid_size + this->kernel_.get_rfp() * REAL(2)) / grid_size;
        auto traj_scaled = traj;
        std::transform(traj_scaled.begin(), traj_scaled.end(), traj_scaled.begin(), [win_length](auto point) { return point * win_length; });
        std::transform(traj_scaled.begin(), traj_scaled.end(), traj_scaled.begin(), [grid_size](auto point) { return (point + REAL(0.5)) * grid_size; });

        // compute convolution matrices
        conv_matrix_.reserve(this->num_frames_);
        conv_matrix_T_.reserve(this->num_frames_);
        for (auto traj : NDArrayViewRange<hoNDArray<vector_td<REAL, D>>>(traj_scaled, 0))
        {
            conv_matrix_.push_back(SDC_internal::make_conv_matrix(traj, this->kernel_));
            conv_matrix_T_.push_back(SDC_internal::transpose(conv_matrix_.back()));
        }
    }

    template<class REAL, unsigned int D>
    void hoSDC_impl<REAL, D>::convolve(const hoNDArray<REAL>& in, hoNDArray<REAL>& out, SDC_conv_mode mode)
    {
        switch (mode)
        {
            case SDC_conv_mode::NC2C:   convolve_NC2C(in, out);     break;
            case SDC_conv_mode::C2NC:   convolve_C2NC(in, out);     break;
        }
    }

    template<class REAL, unsigned int D>
    void hoSDC_impl<REAL, D>::update(const hoNDArray<REAL>& in, hoNDArray<REAL>& out)
    {
        std::transform(out.begin(), out.end(), in.begin(), out.begin(), SDC_internal::safe_divides<REAL>());
    }

    template<class REAL, unsigned int D>
    void hoSDC_impl<REAL, D>::convolve_C2NC(const hoNDArray<REAL>& in, hoNDArray<REAL>& out)
    {
        size_t num_batches = in.get_number_of_elements() / conv_matrix_.front().n_rows;
        assert(num_batches == out.get_number_of_elements() / conv_matrix_.front().n_cols);

        clear(&out);

        #pragma omp parallel for
        for (int b = 0; b < (int)num_batches; b++)
        {
            const REAL* in_view = in.get_data_ptr() + b * conv_matrix_.front().n_rows;
            REAL* out_view = out.get_data_ptr() + b * conv_matrix_.front().n_cols;
            size_t matrix_index = b % conv_matrix_.size();
            matrix_vector_multiply(conv_matrix_[matrix_index], in_view, out_view);
        }
    }

    template<class REAL, unsigned int D>
    void hoSDC_impl<REAL, D>::convolve_NC2C(const hoNDArray<REAL>& in, hoNDArray<REAL>& out)
    {
        size_t num_batches = out.get_number_of_elements() / conv_matrix_.front().n_rows;
        assert(num_batches == in.get_number_of_elements() / conv_matrix_.front().n_cols);

        clear(&out);

        #pragma omp parallel for
        for (int b = 0; b < (int)num_batches; b++)
        {
            REAL *out_view = out.get_data_ptr() + b * conv_matrix_.front().n_rows;
            const REAL *in_view = in.get_data_ptr() + b * conv_matrix_.front().n_cols;
            size_t matrix_index = b % conv_matrix_.size();
            matrix_vector_multiply(conv_matrix_T_[matrix_index], in_view, out_view);
        }
    }

    template<class REAL, unsigned int D>
    boost::shared_ptr<hoNDArray<REAL>> estimate_dcw(
        const hoNDArray<vector_td<REAL, D>>& traj,
        const vector_td<size_t, D>& matrix_size,
        REAL os_factor,
        size_t num_iterations)
    { 
        hoSDC_impl<REAL, D> impl(matrix_size, os_factor, num_iterations);
        return impl.compute(traj);
    }

    template<class REAL, unsigned int D>
    boost::shared_ptr<hoNDArray<REAL>> estimate_dcw(
        const hoNDArray<vector_td<REAL, D>>& traj,
        const hoNDArray<REAL>& initial_dcw,
        const vector_td<size_t, D>& matrix_size,
        REAL os_factor,
        size_t num_iterations)
    { 
        hoSDC_impl<REAL, D> impl(matrix_size, os_factor, num_iterations);
        return impl.compute(traj, initial_dcw);
    }

}   // namespace Gadgetron


template EXPORTSDC boost::shared_ptr<Gadgetron::hoNDArray<float>> Gadgetron::estimate_dcw<float, 2>(
    const Gadgetron::hoNDArray<Gadgetron::vector_td<float, 2>>& traj,
    const Gadgetron::vector_td<size_t, 2>& matrix_size,
    float os_factor,
    size_t num_iterations);

template EXPORTSDC boost::shared_ptr<Gadgetron::hoNDArray<float>> Gadgetron::estimate_dcw<float, 3>(
    const Gadgetron::hoNDArray<Gadgetron::vector_td<float, 3>>& traj,
    const Gadgetron::vector_td<size_t, 3>& matrix_size,
    float os_factor,
    size_t num_iterations);

template EXPORTSDC boost::shared_ptr<Gadgetron::hoNDArray<double>> Gadgetron::estimate_dcw<double, 2>(
    const Gadgetron::hoNDArray<Gadgetron::vector_td<double, 2>>& traj,
    const Gadgetron::vector_td<size_t, 2>& matrix_size,
    double os_factor,
    size_t num_iterations);

template EXPORTSDC boost::shared_ptr<Gadgetron::hoNDArray<double>> Gadgetron::estimate_dcw<double, 3>(
    const Gadgetron::hoNDArray<Gadgetron::vector_td<double, 3>>& traj,
    const Gadgetron::vector_td<size_t, 3>& matrix_size,
    double os_factor,
    size_t num_iterations);

template EXPORTSDC boost::shared_ptr<Gadgetron::hoNDArray<float>> Gadgetron::estimate_dcw<float, 2>(
    const Gadgetron::hoNDArray<Gadgetron::vector_td<float, 2>>& traj,
    const Gadgetron::hoNDArray<float>& initial_dcw,
    const Gadgetron::vector_td<size_t, 2>& matrix_size,
    float os_factor,
    size_t num_iterations);

template EXPORTSDC boost::shared_ptr<Gadgetron::hoNDArray<float>> Gadgetron::estimate_dcw<float, 3>(
    const Gadgetron::hoNDArray<Gadgetron::vector_td<float, 3>>& traj,
    const Gadgetron::hoNDArray<float>& initial_dcw,
    const Gadgetron::vector_td<size_t, 3>& matrix_size,
    float os_factor,
    size_t num_iterations);

template EXPORTSDC boost::shared_ptr<Gadgetron::hoNDArray<double>> Gadgetron::estimate_dcw<double, 2>(
    const Gadgetron::hoNDArray<Gadgetron::vector_td<double, 2>>& traj,
    const Gadgetron::hoNDArray<double>& initial_dcw,
    const Gadgetron::vector_td<size_t, 2>& matrix_size,
    double os_factor,
    size_t num_iterations);

template EXPORTSDC boost::shared_ptr<Gadgetron::hoNDArray<double>> Gadgetron::estimate_dcw<double, 3>(
    const Gadgetron::hoNDArray<Gadgetron::vector_td<double, 3>>& traj,
    const Gadgetron::hoNDArray<double>& initial_dcw,
    const Gadgetron::vector_td<size_t, 3>& matrix_size,
    double os_factor,
    size_t num_iterations);
