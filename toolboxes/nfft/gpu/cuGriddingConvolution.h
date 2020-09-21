#pragma once 

#include "GriddingConvolution.h"

#include "cuSparseMatrix.h"
#include "cuNDArray.h"

/**
 * \brief Convolution type.
 * 
 * This is outside namespace `Gadgetron` for compatibility reasons.
 */
enum class ConvolutionType
{
    STANDARD,
    ATOMIC,
    SPARSE_MATRIX
};

namespace Gadgetron
{
    // Forward declarations.
    template<class T, unsigned int D, template<class, unsigned int> class K>
    struct Convolver;

    template<class T, unsigned int D, template<class, unsigned int> class K, ConvolutionType C>
    struct ConvolverC2NC;

    template<class T, unsigned int D, template<class, unsigned int> class K, ConvolutionType C>
    struct ConvolverNC2C;

    /**
     * \brief Gridding convolution (GPU implementation).
     * 
     * \tparam T Value type. Can be real or complex.
     * \tparam D Number of dimensions.
     * \tparam K Convolution kernel type.
     */
    template<class T, unsigned int D, template<class, unsigned int> class K>
    class cuGriddingConvolution : public GriddingConvolutionBase<cuNDArray, T, D, K>
    {
        using REAL = realType_t<T>;

    public:

        /**
         * \brief Constructor.
         * 
         * \param matrix_size Matrix size.
         * \param matrix_size_os Matrix size with oversampling.
         * \param kernel Convolution kernel.
         * \param conv_type Convolution type.
         */
        cuGriddingConvolution(const vector_td<size_t, D>& matrix_size,
                              const vector_td<size_t, D>& matrix_size_os,
                              const K<REAL, D>& kernel,
                              ConvolutionType conv_type = ConvolutionType::STANDARD);

        /**
         * \brief Constructor.
         * 
         * \param matrix_size Matrix size.
         * \param os_factor Oversampling factor.
         * \param kernel Convolution kernel.
         * \param conv_type Convolution type.
         */
        cuGriddingConvolution(const vector_td<size_t, D>& matrix_size,
                              REAL os_factor,
                              const K<REAL, D>& kernel,
                              ConvolutionType conv_type = ConvolutionType::STANDARD);

        ~cuGriddingConvolution();

        /**
         * \brief Prepare gridding convolution.
         * 
         * \param trajectory Trajectory, normalized to [-0.5, 0.5].
         * \param prep_mode Preparation mode.
         */
        virtual void preprocess(
            const cuNDArray<vector_td<REAL, D>>& trajectory, 
            GriddingConvolutionPrepMode prep_mode = GriddingConvolutionPrepMode::ALL) override;


        /**
         * \brief Get pointer to kernel on device.
         *
         * \return const K<REAL, D>* Device pointer to kernel.
         * \warning This is a device pointer which must not be dereferenced in
         *          host code.
         */
        const K<REAL, D>* get_kernel_d() const;


    private:      

        /**
         * \brief Compute gridding convolution (Cartesian to non-Cartesian).
         * 
         * \param[in] image Image.
         * \param[out] samples Non-Cartesian samples.
         * \param[in] accumulate If true, accumulate the result to output array.
         *                       If false, overwrite output array.
         */
        void compute_C2NC(const cuNDArray<T>& image,
                          cuNDArray<T>& samples,
                          bool accumulate) override;

        /**
         * \brief Compute gridding convolution (non-Cartesian to Cartesian).
         * 
         * \param[in] samples Non-Cartesian samples.
         * \param[out] image Image.
         * \param[in] accumulate If true, accumulate the result to output array.
         *                       If false, overwrite output array.
         */
        void compute_NC2C(const cuNDArray<T> &samples,
                          cuNDArray<T> &image,
                          bool accumulate) override;

        /**
         * \brief Basic object initialization.
         * 
         * \param conv_type Convolution type.
         */
        void initialize(ConvolutionType conv_type);

        /**
         * \brief Check inputs to "compute" functions.
         * 
         * \param samples Non-Cartesian samples.
         * \param image Cartesian image.
         */
        void check_inputs(const cuNDArray<T>& samples,
                          const cuNDArray<T>& image);

        /**
         * \brief Pointer to kernel.
         * 
         * This is a device pointer and must not be dereferenced in host code.
         */
        K<REAL, D>* d_kernel_;

        /**
         * \brief Convolution implementation (Cartesian to Non-Cartesian).
         */
        std::unique_ptr<Convolver<T, D, K>> conv_C2NC_;

        /**
         * \brief Convolution implementation (Non-Cartesian to Cartesian).
         */
        std::unique_ptr<Convolver<T, D, K>> conv_NC2C_;

        /**
         * \brief Non-Cartesian sample coordinates.
         */
        thrust::device_vector<vector_td<REAL, D>> trajectory_;

        /**
         * \brief Active device.
         */
        int device_;

        /**
         * \brief Warp size as a power of 2.
         */
        unsigned int warp_size_power_;

        /**
         * \brief Matrix padding.
         */
        vector_td<size_t, D> matrix_padding_;       

        friend ConvolverC2NC<T, D, K, ConvolutionType::STANDARD>;
        friend ConvolverNC2C<T, D, K, ConvolutionType::STANDARD>;

        friend ConvolverC2NC<T, D, K, ConvolutionType::ATOMIC>;
        friend ConvolverNC2C<T, D, K, ConvolutionType::ATOMIC>;

        friend ConvolverC2NC<T, D, K, ConvolutionType::SPARSE_MATRIX>;
        friend ConvolverNC2C<T, D, K, ConvolutionType::SPARSE_MATRIX>;
    };

    template<class T, unsigned int D, template<class, unsigned int> class K>
    const K<realType_t<T>, D>* cuGriddingConvolution<T, D, K>::get_kernel_d() const
    {
        return d_kernel_;
    }

    template<class T, unsigned int D, template<class, unsigned int> class K>
    struct Convolver
    {
        using REAL = realType_t<T>;

    public:

        Convolver(cuGriddingConvolution<T, D, K>& plan)
          : plan_(plan) { };
        
        virtual ~Convolver() { };

        virtual void prepare(
            const thrust::device_vector<vector_td<REAL, D>>& trajectory) = 0;

        virtual void compute(
            const cuNDArray<T>& input,
            cuNDArray<T>& output,
            bool accumulate) = 0;

    protected:

        cuGriddingConvolution<T, D, K>& plan_;
    };


    template<class T, unsigned int D, template<class, unsigned int> class K, ConvolutionType C>
    struct ConvolverC2NC
      : public Convolver<T, D, K>
    {
        using REAL = realType_t<T>;

    public:

        ConvolverC2NC(cuGriddingConvolution<T, D, K>& plan)
          : Convolver<T, D, K>(plan) { };

        void prepare(
            const thrust::device_vector<vector_td<REAL, D>>& trajectory) override;

        void compute(
            const cuNDArray<T>& image,
            cuNDArray<T>& samples,
            bool accumulate) override;
    };


    template<class T, unsigned int D, template<class, unsigned int> class K>
    class ConvolverNC2C<T, D, K, ConvolutionType::STANDARD>
      : public Convolver<T, D, K>
    {
        using REAL = realType_t<T>;

    public:      

        ConvolverNC2C(cuGriddingConvolution<T, D, K>& plan)
          : Convolver<T, D, K>(plan) { };

        void prepare(
            const thrust::device_vector<vector_td<REAL, D>>& trajectory) override;

        void compute(
            const cuNDArray<T>& samples,
            cuNDArray<T>& image,
            bool accumulate) override;

    private:

        void wrap_image(
            cuNDArray<T>& source,
            cuNDArray<T>& target,
            bool accumulate);

        thrust::device_vector<unsigned int> tuples_first, tuples_last;
        thrust::device_vector<unsigned int> bucket_begin, bucket_end;
    };


    template<class T, unsigned int D, template<class, unsigned int> class K>
    class ConvolverNC2C<T, D, K, ConvolutionType::ATOMIC>
      : public Convolver<T, D, K>
    {
        using REAL = realType_t<T>;

    public:

        ConvolverNC2C(cuGriddingConvolution<T, D, K>& plan)
          : Convolver<T, D, K>(plan) { };

        void prepare(
            const thrust::device_vector<vector_td<REAL, D>>& trajectory) override;

        void compute(
            const cuNDArray<T>& samples,
            cuNDArray<T>& image,
            bool accumulate) override;
    };


    template<class T, unsigned int D, template<class, unsigned int> class K>
    class ConvolverNC2C<T, D, K, ConvolutionType::SPARSE_MATRIX>
      : public Convolver<T, D, K>
    {
        using REAL = realType_t<T>;

    public:

        ConvolverNC2C(cuGriddingConvolution<T, D, K>& plan)
          : Convolver<T, D, K>(plan) { };

        void prepare(
            const thrust::device_vector<vector_td<REAL, D>>& trajectory) override;

        void compute(
            const cuNDArray<T>& samples,
            cuNDArray<T>& image,
            bool accumulate) override;

    private:

        std::unique_ptr<cuCsrMatrix<T>> conv_matrix_;
    };

    /**
     * \brief Gridding convolution factory (GPU specialization).
     * 
     * \tparam T Value type. Can be real or complex.
     * \tparam D Number of dimensions.
     * \tparam K Convolution kernel type.
     */
    template<class T, unsigned int D, template<class, unsigned int> class K>
    struct GriddingConvolution<cuNDArray, T, D, K>
    {
        /**
         * \brief Make a new gridding convolution object.
         * 
         * \param matrix_size Matrix size.
         * \param matrix_size_os Matrix size with oversampling.
         * \param kernel Convolution kernel.
         * \param conv_type Convolution type.
         * \return std::unique_ptr<cuGriddingConvolution<T, D, K>> 
         */
        static std::unique_ptr<cuGriddingConvolution<T, D, K>> make(
            const vector_td<size_t, D>& matrix_size,
            const vector_td<size_t, D>& matrix_size_os,
            const K<realType_t<T>, D>& kernel,
            ConvolutionType conv_type = ConvolutionType::STANDARD)
        {
            return std::make_unique<cuGriddingConvolution<T, D, K>>(
                matrix_size, matrix_size_os, kernel, conv_type);
        }

        /**
         * \brief Make a new gridding convolution object.
         * 
         * \param matrix_size Matrix size.
         * \param matrix_size_os Matrix size with oversampling.
         * \param kernel Convolution kernel.
         * \param conv_type Convolution type.
         * \return std::unique_ptr<cuGriddingConvolution<T, D, K>> 
         */
        static std::unique_ptr<cuGriddingConvolution<T, D, K>> make(
            const vector_td<size_t, D>& matrix_size,
            realType_t<T> os_factor,
            const K<realType_t<T>, D>& kernel,
            ConvolutionType conv_type = ConvolutionType::STANDARD)
        {
            return std::make_unique<cuGriddingConvolution<T, D, K>>(
                matrix_size, os_factor, kernel, conv_type);
        }
    };

} // namespace Gadgetron
