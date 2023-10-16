
#pragma once

#include "core_defines.h"
#include "vector_td.h"

namespace Gadgetron
{
    /**
     * \brief Convolution kernel.
     *
     * This object is usable in both host and device code. 
     *
     * This is an abstract base class and should not be used directly.
     * All kernels for use in a `GriddingConvolution` should derive from this
     * class.
     *
     * A CRTP-based static dispatch mechanism is used for maximum performance.
     * 
     * \tparam REAL Value type. Must be a real type.
     * \tparam D Number of dimensions.
     * \tparam K Kernel type.
     */
    template<class REAL, unsigned int D, template<class, unsigned int> class K>
    class ConvolutionKernel
    {
    public:
        
        /**
         * \brief Constructor.
         * 
         * \param width Kernel width. 
         */
        __host__ __device__
        ConvolutionKernel(REAL width);

        /**
         * \brief Destructor.
         */
        __host__ __device__
        ~ConvolutionKernel();

        /**
         * \brief Get kernel value at given coordinates.
         * 
         * \param u Coordinates.
         * \return REAL Kernel value.
         */
        __host__ __device__
        REAL get(const vector_td<REAL, D>& u) const;

        /**
         * \brief Get kernel value at given radius along axis.
         * 
         * \param r Radius.
         * \param ax Axis.
         * \return REAL Kernel value.
         */
        __host__ __device__
        REAL get(REAL r, size_t ax = 0) const;

        /**
         * \brief Get kernel width.
         */
        __host__ __device__
        REAL get_width() const;

        /**
         * \brief Get kernel radius.
         */
        __host__ __device__
        REAL get_radius() const;

    protected:

        /**
         * \brief Compute kernel value at given coordinates.
         *
         * Calls derived class implementation.
         * 
         * \param u Coordinates.
         * \return REAL Kernel value. 
         */
        __host__ __device__
        REAL compute(const vector_td<REAL, D>& u) const;

        /**
         * \brief Compute kernel value at given radius along axis.
         *
         * Calls derived class implementation.
         * 
         * \param r Radius.
         * \param ax Axis.
         * \return REAL Kernel value. 
         */
        __host__ __device__
        REAL compute(REAL r, size_t ax = 0) const;

        /**
         * \brief Look up kernel value at given coordinates.
         *
         * Calls derived class implementation.
         * 
         * \param u Coordinates.
         * \return REAL Kernel value. 
         */
        __host__ __device__
        REAL lookup(const vector_td<REAL, D>& u) const;

        /**
         * \brief Look up kernel value at given radius along axis.
         *
         * Calls derived class implementation.
         * 
         * \param r Radius.
         * \param ax Axis.
         * \return REAL Kernel value. 
         */
        __host__ __device__
        REAL lookup(REAL r, size_t ax = 0) const;

        /**
         * \brief Kernel width.
         * 
         * This is equal to 2 * radius_.
         */
        REAL width_;
        
        /**
         * \brief Kernel radius.
         *
         * This is equal to width_ / 2.
         */
        REAL radius_;
    };


    /**
     * \brief Kaiser-Bessel kernel.
     * 
     * \tparam REAL Value type. Must be a real type.
     * \tparam D Number of dimensions.
     */
    template<class REAL, unsigned int D>
    class KaiserKernel : public ConvolutionKernel<REAL, D, KaiserKernel>
    {
    public:

        /**
         * \brief Constructor.
         * 
         * \param matrix_size Matrix size.
         * \param matrix_size_os Matrix size with oversampling.
         * \param width Kernel width.
         */
        __host__ __device__
        KaiserKernel(const vector_td<unsigned int, D>& matrix_size,
                     const vector_td<unsigned int, D>& matrix_size_os,
                     REAL width);

        /**
         * \brief Constructor.
         * 
         * \param matrix_size Matrix size.
         * \param os_factor Matrix size with oversampling.
         * \param width Kernel width.
         */
        __host__ __device__
        KaiserKernel(const vector_td<unsigned int, D>& matrix_size,
                     REAL os_factor,
                     REAL width);

        /**
         * \brief Get the beta factor.
         */
        __host__ __device__
        vector_td<REAL, D> get_beta() const;

    private:

        /**
         * \brief Compute kernel value at given coordinates.
         * 
         * \param u Coordinates.
         * \return REAL Kernel value. 
         */
        __host__ __device__
        REAL compute(const vector_td<REAL, D>& u) const;

        /**
         * \brief Compute kernel value at given radius along axis.
         * 
         * \param r Radius.
         * \param ax Axis.
         * \return REAL Kernel value. 
         */
        __host__ __device__
        REAL compute(REAL r, size_t ax = 0) const;

        /**
         * \brief Look up kernel value at given coordinates.
         *
         * \param u Coordinates.
         * \return REAL Kernel value.
         * \warning Not implemented yet.
         */
        __host__ __device__
        REAL lookup(const vector_td<REAL, D>& u) const;

        /**
         * \brief Look up kernel value at given radius along axis.
         *
         * \param r Radius.
         * \param ax Axis.
         * \return REAL Kernel value.
         * \warning Not implemented yet.
         */
        __host__ __device__
        REAL lookup(REAL r, size_t ax = 0) const;

        /**
         * \brief Compute beta factor.
         * 
         * Compute beta factor based on the oversampling ratio.
         *
         * \return vector_td<REAL, D> Beta factor. 
         */
        __host__ __device__
        vector_td<REAL, D> compute_beta() const;

        /**
         * \brief Matrix size (gridding convolution).
         */
        vector_td<unsigned int, D> matrix_size_;

        /**
         * \brief Matrix size with oversampling (gridding convolution).
         *
         * This is equal to matrix_size_ * os_factor_.
         */
        vector_td<unsigned int, D> matrix_size_os_;

        /**
         * \brief Oversampling factor (gridding convolution).
         */
        vector_td<REAL, D> os_factor_;

        /**
         * \brief Beta factor.
         */
        vector_td<REAL, D> beta_;

        /**
         * \brief Width inverse.
         *
         * This is equal to 1 / width_.
         */
        REAL width_inv_;

        friend class ConvolutionKernel<REAL, D, KaiserKernel>;
    };


    /**
     * \brief Jinc kernel.
     * 
     * \tparam REAL Value type. Must be a real type.
     * \tparam D Number of dimensions.
     */
    template<class REAL, unsigned int D>
    class JincKernel : public ConvolutionKernel<REAL, D, JincKernel>
    {
    public:

        /**
         * \brief Constructor.
         * 
         * \param kernelWidth Matrix size.
         */
        __host__ __device__
        JincKernel(float kernelWidth);

    private:

        /**
         * \brief Compute kernel value at given coordinates.
         * 
         * \param u Coordinates.
         * \return REAL Kernel value. 
         */
        __host__ __device__
        REAL compute(const vector_td<REAL, D>& u) const;

        /**
         * \brief Compute kernel value at given radius along axis.
         * 
         * Becase this is a circularly symmetric kernel, the `axis` parameter is
         * irrelevant.
         *
         * \param r Radius.
         * \param ax Axis.
         * \return REAL Kernel value. 
         */
        __host__ __device__
        REAL compute(REAL r, size_t ax = 0) const;

        /**
         * \brief Look up kernel value at given coordinates.
         *
         * \param u Coordinates.
         * \return REAL Kernel value.
         * \warning Not implemented yet.
         */
        __host__ __device__
        REAL lookup(const vector_td<REAL, D>& u) const;

        /**
         * \brief Look up kernel value at given radius along axis.
         *
         * \param r Radius.
         * \param ax Axis.
         * \return REAL Kernel value.
         * \warning Not implemented yet.
         */
        __host__ __device__
        REAL lookup(REAL r, size_t ax = 0) const;

        /**
         * \brief Matrix size (gridding convolution).
         */
        vector_td<unsigned int, D> matrix_size_;

        /**
         * \brief Matrix size with oversampling (gridding convolution).
         *
         * This is equal to matrix_size_ * os_factor_.
         */
        vector_td<unsigned int, D> matrix_size_os_;

        /**
         * \brief Oversampling factor (gridding convolution).
         */
        vector_td<REAL, D> os_factor_;

        /**
         * \brief Order of polynomial.
         */
        static constexpr unsigned int poly_order = 5U;

        // This is what we would like to do, but direct use in device code of 
        // const-qualified arrays is not currently supported by CUDA.
        // static constexpr REAL poly_coefs[] = {
        //     0.999923599661866,
        //     0.014023570134078,
        //     -2.969947516849635,
        //     0.891173008769047,
        //     2.326043034495162,
        //     -1.270644932677047
        // };

        // An ugly solution to the above problem, but it keeps nvcc happy.
        static constexpr REAL poly_coef_0 = 0.999923599661866;
        static constexpr REAL poly_coef_1 = 0.014023570134078;
        static constexpr REAL poly_coef_2 = -2.969947516849635;
        static constexpr REAL poly_coef_3 = 0.891173008769047;
        static constexpr REAL poly_coef_4 = 2.326043034495162;
        static constexpr REAL poly_coef_5 = -1.270644932677047;

        /**
         * \brief Normalized kernel radius.
         */
        static constexpr REAL norm_radius = 0.96960938;

        friend class ConvolutionKernel<REAL, D, JincKernel>;
    };

} // namespace Gadgetron

#include "ConvolutionKernel.hpp"