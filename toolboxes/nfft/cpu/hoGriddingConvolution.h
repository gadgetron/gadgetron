#pragma once 

#include "GriddingConvolution.h"

#include "hoNDArray.h"

#include "ConvolutionMatrix.h"

namespace Gadgetron
{
    /**
     * \brief Gridding convolution (CPU implementation).
     * 
     * \tparam T Value type. Can be real or complex.
     * \tparam D Number of dimensions.
     * \tparam K Convolution kernel type.
     */
    template<class T, unsigned int D, template<class, unsigned int> class K>
    class hoGriddingConvolution : public GriddingConvolutionBase<hoNDArray, T, D, K>
    {
        typedef typename realType<T>::Type REAL;

    public:

        /**
         * \brief Constructor.
         * 
         * \param matrix_size Matrix size.
         * \param matrix_size_os Matrix size with oversampling.
         * \param kernel Convolution kernel.
         */
        hoGriddingConvolution(const vector_td<size_t, D>& matrix_size,
                              const vector_td<size_t, D>& matrix_size_os,
                              const K<REAL, D>& kernel);

        /**
         * \brief Constructor.
         * 
         * \param matrix_size Matrix size.
         * \param os_factor Matrix size with oversampling.
         * \param kernel Convolution kernel.
         */
        hoGriddingConvolution(const vector_td<size_t, D>& matrix_size,
                              REAL os_factor,
                              const K<REAL, D>& kernel);

        ~hoGriddingConvolution();

        /**
         * \brief Prepare gridding convolution.
         * 
         * \param trajectory Trajectory, normalized to [-0.5, 0.5].
         * \param prep_mode Preparation mode.
         */
        virtual void preprocess(
            const hoNDArray<vector_td<REAL, D>>& trajectory, 
            GriddingConvolutionPrepMode prep_mode = GriddingConvolutionPrepMode::ALL) override;

    private:

        /**
         * \brief Compute gridding convolution (Cartesian to non-Cartesian).
         * 
         * \param[in] image Image.
         * \param[out] samples Non-Cartesian samples.
         * \param[in] accumulate If true, accumulate the result to output array.
         *                       If false, overwrite output array.
         */
        void compute_C2NC(const hoNDArray<T>& image,
                          hoNDArray<T>& samples,
                          bool accumulate) override;

        /**
         * \brief Compute gridding convolution (non-Cartesian to Cartesian).
         * 
         * \param[in] samples Non-Cartesian samples.
         * \param[out] image Image.
         * \param[in] accumulate If true, accumulate the result to output array.
         *                       If false, overwrite output array.
         */
        void compute_NC2C(const hoNDArray<T> &samples,
                           hoNDArray<T> &image,
                           bool accumulate) override;

        std::vector<ConvInternal::ConvolutionMatrix<REAL>> conv_matrix_;
        std::vector<ConvInternal::ConvolutionMatrix<REAL>> conv_matrix_T_;
    };

    /**
     * \brief Gridding convolution factory (CPU specialization).
     * 
     * \tparam T Value type. Can be real or complex.
     * \tparam D Number of dimensions.
     * \tparam K Convolution kernel type.
     */
    template<class T, unsigned int D, template<class, unsigned int> class K>
    struct GriddingConvolution<hoNDArray, T, D, K>
    {
        /**
         * \brief Make a new gridding convolution object.
         * 
         * \param matrix_size Matrix size.
         * \param matrix_size_os Matrix size with oversampling.
         * \param kernel Convolution kernel.
         * \return std::unique_ptr<hoGriddingConvolution<T, D, K>> 
         */
        static std::unique_ptr<hoGriddingConvolution<T, D, K>> make(
            const vector_td<size_t, D>& matrix_size,
            const vector_td<size_t, D>& matrix_size_os,
            const K<realType_t<T>, D>& kernel)
        {
            return std::make_unique<hoGriddingConvolution<T, D, K>>(
                matrix_size, matrix_size_os, kernel);
        }

        /**
         * \brief Make a new gridding convolution object.
         * 
         * \param matrix_size Matrix size.
         * \param os_factor Oversampling factor.
         * \param kernel Convolution kernel.
         * \return std::unique_ptr<hoGriddingConvolution<T, D, K>> 
         */
        static std::unique_ptr<hoGriddingConvolution<T, D, K>> make(
            const vector_td<size_t, D>& matrix_size,
            realType_t<T> os_factor,
            const K<realType_t<T>, D>& kernel)
        {
            return std::make_unique<hoGriddingConvolution<T, D, K>>(
                matrix_size, os_factor, kernel);
        }
    };
}
