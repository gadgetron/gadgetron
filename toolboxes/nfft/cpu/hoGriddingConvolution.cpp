
#include "hoGriddingConvolution.h"

#include "hoNDArray_elemwise.h"
#include "NDArray_utils.h"

#include "ConvolutionMatrix.h"

namespace Gadgetron
{
    template<class T, unsigned int D, template<class, unsigned int> class K>
    hoGriddingConvolution<T, D, K>::hoGriddingConvolution(
        const vector_td<size_t, D>& matrix_size,
        const vector_td<size_t, D>& matrix_size_os,
        const K<REAL, D>& kernel)
      : GriddingConvolutionBase<hoNDArray, T, D, K>(
          matrix_size, matrix_size_os, kernel)
    {

    }


    template<class T, unsigned int D, template<class, unsigned int> class K>
    hoGriddingConvolution<T, D, K>::hoGriddingConvolution(
        const vector_td<size_t, D>& matrix_size,
        const REAL os_factor,
        const K<REAL, D>& kernel)
      : GriddingConvolutionBase<hoNDArray, T, D, K>(
          matrix_size, os_factor, kernel)
    {

    }


    template<class T, unsigned int D, template<class, unsigned int> class K>
    hoGriddingConvolution<T, D, K>::~hoGriddingConvolution()
    {
        
    }


    template<class T, unsigned int D, template<class, unsigned int> class K>
    void hoGriddingConvolution<T, D, K>::preprocess(
        const hoNDArray<vector_td<REAL, D>> &trajectory,
        GriddingConvolutionPrepMode prep_mode)
    {       
        GriddingConvolutionBase<hoNDArray, T, D, K>::preprocess(
            trajectory, prep_mode);

        auto scaled_trajectory = trajectory;
        auto matrix_size_os_real = vector_td<REAL,D>(this->matrix_size_os_);
        std::transform(scaled_trajectory.begin(),
                       scaled_trajectory.end(),
                       scaled_trajectory.begin(),
                       [matrix_size_os_real](auto point)
                       { return (point + REAL(0.5)) * matrix_size_os_real; });

        conv_matrix_.reserve(this->num_frames_);
        conv_matrix_T_.reserve(this->num_frames_);

        for (auto traj : NDArrayViewRange<hoNDArray<vector_td<REAL,D>>>(
                            scaled_trajectory, 0))
        {
            conv_matrix_.push_back(ConvInternal::make_conv_matrix(
                traj, this->matrix_size_os_, this->kernel_));

            if (prep_mode == GriddingConvolutionPrepMode::NC2C ||
                prep_mode == GriddingConvolutionPrepMode::ALL)
            {
                conv_matrix_T_.push_back(ConvInternal::transpose(conv_matrix_.back()));
            }
        }
    }


    namespace
    {   
        /**
         * \brief Matrix-vector multiplication.
         * 
         * \tparam T Value type. Can be real or complex.
         * \param[in] matrix Convolution matrix.
         * \param[in] vector Vector.
         * \param[out] result Operation result.
         */
        template<class T>
        void mvm(
            const ConvInternal::ConvolutionMatrix<realType_t<T>>& matrix,
            const T* vector,
            T* result)
        {
            for (size_t i = 0; i < matrix.n_cols; i++)
            {
                auto &row_indices = matrix.indices[i];
                auto &weights = matrix.weights[i];

                for (size_t n = 0; n < row_indices.size(); n++)
                {
                    result[i] += vector[row_indices[n]] * weights[n];
                }
            }
        }
    }


    template<class T, unsigned int D, template<class, unsigned int> class K>
    void hoGriddingConvolution<T, D, K>::compute_C2NC(
        const hoNDArray<T> &image,
        hoNDArray<T> &samples,
        bool accumulate)
    {
        size_t nbatches = image.get_number_of_elements() / conv_matrix_.front().n_rows;
        assert(nbatches == samples.get_number_of_elements() / conv_matrix_.front().n_cols);

        if (!accumulate) clear(&samples);

        #pragma omp parallel for
        for (int b = 0; b < (int)nbatches; b++)
        {
            const T* image_view = image.get_data_ptr() + b * conv_matrix_.front().n_rows;
            T* samples_view = samples.get_data_ptr() + b * conv_matrix_.front().n_cols;
            size_t matrix_index = b % conv_matrix_.size();
            mvm(conv_matrix_[matrix_index], image_view, samples_view);
        }
    }


    template<class T, unsigned int D, template<class, unsigned int> class K>
    void hoGriddingConvolution<T, D, K>::compute_NC2C(
        const hoNDArray<T> &samples,
        hoNDArray<T> &image,
        bool accumulate)
    {
        size_t nbatches = image.get_number_of_elements() / conv_matrix_.front().n_rows;
        assert(nbatches == samples.get_number_of_elements() / conv_matrix_.front().n_cols);

        if (!accumulate) clear(&image);

        #pragma omp parallel for
        for (int b = 0; b < (int)nbatches; b++)
        {
            T* image_view = image.get_data_ptr() + b * conv_matrix_.front().n_rows;
            const T* samples_view = samples.get_data_ptr() + b * conv_matrix_.front().n_cols;
            size_t matrix_index = b % conv_matrix_.size();
            mvm(conv_matrix_T_[matrix_index], samples_view, image_view);
        }
    }
}

template class Gadgetron::hoGriddingConvolution<float, 1, Gadgetron::KaiserKernel>;
template class Gadgetron::hoGriddingConvolution<float, 2, Gadgetron::KaiserKernel>;
template class Gadgetron::hoGriddingConvolution<float, 3, Gadgetron::KaiserKernel>;
template class Gadgetron::hoGriddingConvolution<float, 4, Gadgetron::KaiserKernel>;

template class Gadgetron::hoGriddingConvolution<double, 1, Gadgetron::KaiserKernel>;
template class Gadgetron::hoGriddingConvolution<double, 2, Gadgetron::KaiserKernel>;
template class Gadgetron::hoGriddingConvolution<double, 3, Gadgetron::KaiserKernel>;
template class Gadgetron::hoGriddingConvolution<double, 4, Gadgetron::KaiserKernel>;

template class Gadgetron::hoGriddingConvolution<Gadgetron::complext<float>, 1, Gadgetron::KaiserKernel>;
template class Gadgetron::hoGriddingConvolution<Gadgetron::complext<float>, 2, Gadgetron::KaiserKernel>;
template class Gadgetron::hoGriddingConvolution<Gadgetron::complext<float>, 3, Gadgetron::KaiserKernel>;
template class Gadgetron::hoGriddingConvolution<Gadgetron::complext<float>, 4, Gadgetron::KaiserKernel>;

template class Gadgetron::hoGriddingConvolution<Gadgetron::complext<double>, 1, Gadgetron::KaiserKernel>;
template class Gadgetron::hoGriddingConvolution<Gadgetron::complext<double>, 2, Gadgetron::KaiserKernel>;
template class Gadgetron::hoGriddingConvolution<Gadgetron::complext<double>, 3, Gadgetron::KaiserKernel>;
template class Gadgetron::hoGriddingConvolution<Gadgetron::complext<double>, 4, Gadgetron::KaiserKernel>;

template class Gadgetron::hoGriddingConvolution<float, 1, Gadgetron::JincKernel>;
template class Gadgetron::hoGriddingConvolution<float, 2, Gadgetron::JincKernel>;
template class Gadgetron::hoGriddingConvolution<float, 3, Gadgetron::JincKernel>;
template class Gadgetron::hoGriddingConvolution<float, 4, Gadgetron::JincKernel>;

template class Gadgetron::hoGriddingConvolution<double, 1, Gadgetron::JincKernel>;
template class Gadgetron::hoGriddingConvolution<double, 2, Gadgetron::JincKernel>;
template class Gadgetron::hoGriddingConvolution<double, 3, Gadgetron::JincKernel>;
template class Gadgetron::hoGriddingConvolution<double, 4, Gadgetron::JincKernel>;

template class Gadgetron::hoGriddingConvolution<Gadgetron::complext<float>, 1, Gadgetron::JincKernel>;
template class Gadgetron::hoGriddingConvolution<Gadgetron::complext<float>, 2, Gadgetron::JincKernel>;
template class Gadgetron::hoGriddingConvolution<Gadgetron::complext<float>, 3, Gadgetron::JincKernel>;
template class Gadgetron::hoGriddingConvolution<Gadgetron::complext<float>, 4, Gadgetron::JincKernel>;

template class Gadgetron::hoGriddingConvolution<Gadgetron::complext<double>, 1, Gadgetron::JincKernel>;
template class Gadgetron::hoGriddingConvolution<Gadgetron::complext<double>, 2, Gadgetron::JincKernel>;
template class Gadgetron::hoGriddingConvolution<Gadgetron::complext<double>, 3, Gadgetron::JincKernel>;
template class Gadgetron::hoGriddingConvolution<Gadgetron::complext<double>, 4, Gadgetron::JincKernel>;

template class Gadgetron::GriddingConvolution<Gadgetron::hoNDArray, float, 1, Gadgetron::KaiserKernel>;
template class Gadgetron::GriddingConvolution<Gadgetron::hoNDArray, float, 2, Gadgetron::KaiserKernel>;
template class Gadgetron::GriddingConvolution<Gadgetron::hoNDArray, float, 3, Gadgetron::KaiserKernel>;
template class Gadgetron::GriddingConvolution<Gadgetron::hoNDArray, float, 4, Gadgetron::KaiserKernel>;

template class Gadgetron::GriddingConvolution<Gadgetron::hoNDArray, double, 1, Gadgetron::KaiserKernel>;
template class Gadgetron::GriddingConvolution<Gadgetron::hoNDArray, double, 2, Gadgetron::KaiserKernel>;
template class Gadgetron::GriddingConvolution<Gadgetron::hoNDArray, double, 3, Gadgetron::KaiserKernel>;
template class Gadgetron::GriddingConvolution<Gadgetron::hoNDArray, double, 4, Gadgetron::KaiserKernel>;

template class Gadgetron::GriddingConvolution<Gadgetron::hoNDArray, Gadgetron::complext<float>, 1, Gadgetron::KaiserKernel>;
template class Gadgetron::GriddingConvolution<Gadgetron::hoNDArray, Gadgetron::complext<float>, 2, Gadgetron::KaiserKernel>;
template class Gadgetron::GriddingConvolution<Gadgetron::hoNDArray, Gadgetron::complext<float>, 3, Gadgetron::KaiserKernel>;
template class Gadgetron::GriddingConvolution<Gadgetron::hoNDArray, Gadgetron::complext<float>, 4, Gadgetron::KaiserKernel>;

template class Gadgetron::GriddingConvolution<Gadgetron::hoNDArray, Gadgetron::complext<double>, 1, Gadgetron::KaiserKernel>;
template class Gadgetron::GriddingConvolution<Gadgetron::hoNDArray, Gadgetron::complext<double>, 2, Gadgetron::KaiserKernel>;
template class Gadgetron::GriddingConvolution<Gadgetron::hoNDArray, Gadgetron::complext<double>, 3, Gadgetron::KaiserKernel>;
template class Gadgetron::GriddingConvolution<Gadgetron::hoNDArray, Gadgetron::complext<double>, 4, Gadgetron::KaiserKernel>;

template class Gadgetron::GriddingConvolution<Gadgetron::hoNDArray, float, 1, Gadgetron::JincKernel>;
template class Gadgetron::GriddingConvolution<Gadgetron::hoNDArray, float, 2, Gadgetron::JincKernel>;
template class Gadgetron::GriddingConvolution<Gadgetron::hoNDArray, float, 3, Gadgetron::JincKernel>;
template class Gadgetron::GriddingConvolution<Gadgetron::hoNDArray, float, 4, Gadgetron::JincKernel>;

template class Gadgetron::GriddingConvolution<Gadgetron::hoNDArray, double, 1, Gadgetron::JincKernel>;
template class Gadgetron::GriddingConvolution<Gadgetron::hoNDArray, double, 2, Gadgetron::JincKernel>;
template class Gadgetron::GriddingConvolution<Gadgetron::hoNDArray, double, 3, Gadgetron::JincKernel>;
template class Gadgetron::GriddingConvolution<Gadgetron::hoNDArray, double, 4, Gadgetron::JincKernel>;

template class Gadgetron::GriddingConvolution<Gadgetron::hoNDArray, Gadgetron::complext<float>, 1, Gadgetron::JincKernel>;
template class Gadgetron::GriddingConvolution<Gadgetron::hoNDArray, Gadgetron::complext<float>, 2, Gadgetron::JincKernel>;
template class Gadgetron::GriddingConvolution<Gadgetron::hoNDArray, Gadgetron::complext<float>, 3, Gadgetron::JincKernel>;
template class Gadgetron::GriddingConvolution<Gadgetron::hoNDArray, Gadgetron::complext<float>, 4, Gadgetron::JincKernel>;

template class Gadgetron::GriddingConvolution<Gadgetron::hoNDArray, Gadgetron::complext<double>, 1, Gadgetron::JincKernel>;
template class Gadgetron::GriddingConvolution<Gadgetron::hoNDArray, Gadgetron::complext<double>, 2, Gadgetron::JincKernel>;
template class Gadgetron::GriddingConvolution<Gadgetron::hoNDArray, Gadgetron::complext<double>, 3, Gadgetron::JincKernel>;
template class Gadgetron::GriddingConvolution<Gadgetron::hoNDArray, Gadgetron::complext<double>, 4, Gadgetron::JincKernel>;
