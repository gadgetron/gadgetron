#pragma once

#include <vector>

#include "complext.h"
#include "vector_td.h"

#include "GriddingConvolution.h"

/**
  Enum specifying the direction of the NFFT standalone FFT.
*/
namespace Gadgetron {
    enum class NFFT_fft_mode {
        FORWARDS, /**< forwards FFT. */
        BACKWARDS /**< backwards FFT. */
    };

    /**
      Enum specifying the direction of the NFFT standalone convolution
   */
    enum class NFFT_conv_mode {
        C2NC, /**< convolution: Cartesian to non-Cartesian. */
        NC2C /**< convolution: non-Cartesian to Cartesian. */
    };


/**
     Enum defining the desired NFFT operation
  */
    enum class NFFT_comp_mode {
        FORWARDS_C2NC, /**< forwards NFFT Cartesian to non-Cartesian. */
        FORWARDS_NC2C, /**< forwards NFFT non-Cartesian to Cartesian. */
        BACKWARDS_C2NC, /**< backwards NFFT Cartesian to non-Cartesian. */
        BACKWARDS_NC2C /**< backwards NFFT non-Cartesian to Cartesian. */
    };

/**
   Enum to specify the preprocessing mode.
*/
    enum class NFFT_prep_mode {
        C2NC, /**< preprocess to perform a Cartesian to non-Cartesian NFFT. */
        NC2C, /**< preprocess to perform a non-Cartesian to Cartesian NFFT. */
        ALL /**< preprocess to perform NFFTs in both directions. */
    };


    template<template<class> class ARRAY, class REAL, unsigned int D>
    struct NFFT {

    };


    template<template<class> class ARRAY, class REAL, unsigned int D>
    class NFFT_plan {
    public:


        NFFT_plan(const vector_td<size_t,D>& matrix_size, const vector_td<size_t,D>& matrix_size_os, REAL W);
        NFFT_plan(const vector_td<size_t,D>& matrix_size, REAL oversampling_factor, REAL W);


        void reconfigure(const vector_td<size_t,D>& matrix_size, const vector_td<size_t,D>& matrix_size_os, REAL W);
        void reconfigure(const vector_td<size_t,D>& matrix_size, REAL oversampling_Factor,  REAL W);
        /**
           Perform NFFT preprocessing for a given trajectory.
           \param trajectory the NFFT non-Cartesian trajectory normalized to the range [-1/2;1/2].
           \param mode enum class specifying the preprocessing mode
        */
        virtual void preprocess(const ARRAY<vector_td<REAL,D>> &trajectory,
                                NFFT_prep_mode mode = NFFT_prep_mode::ALL);


        /**
           Execute the NFFT.
           \param[in] in the input array.
           \param[out] out the output array.
           \param mode enum class specifying the mode of operation.
        */
        virtual void compute(const ARRAY<complext<REAL>>& in, ARRAY<complext<REAL>>&out,
        const ARRAY<REAL> *dcw, NFFT_comp_mode mode);

        /**
           Execute an NFFT iteraion (from Cartesian image space to non-Cartesian Fourier space and back to Cartesian image space).
           \param[in] in the input array.
           \param[out] out the output array.
           \param[in] dcw optional density compensation weights weighing the input samples according to the sampling density.
           If an 0x0-pointer is provided no density compensation is used.
           \param[in] halfway_dims specifies the dimensions of the intermediate Fourier space (codomain).
        */
        virtual void mult_MH_M(const ARRAY<complext<REAL>>& in, ARRAY<complext < REAL>>& out,
        const ARRAY<REAL> *dcw
        );

    public: // Utilities


        /**
           Perform "standalone" convolution
           \param[in] in the input array.
           \param[out] out the output array.
           \param[in] dcw optional density compensation weights.
           \param[in] mode enum class specifying the mode of the convolution
           \param[in] accumulate specifies whether the result is added to the output (accumulation) or if the output is overwritten.
        */
        virtual void convolve(const ARRAY<complext<REAL>>& in, ARRAY<complext < REAL>>& out,
        NFFT_conv_mode mode,
        bool accumulate = false);


        /**
           Cartesian FFT. For completeness, just invokes the cuNDFFT class.
           \param[in,out] data the data for the inplace FFT.
           \param mode enum class specifying the direction of the FFT.
           \param do_scale boolean specifying whether FFT normalization is desired.
        */
        virtual void fft(ARRAY<complext < REAL>>& data,
        NFFT_fft_mode mode,
        bool do_scale = true
        ) = 0;

        /**
           NFFT deapodization.
           \param[in,out] image the image to be deapodized (inplace).
        */
        virtual void deapodize(ARRAY<complext < REAL>>

        & image,
        bool fourier_domain = false
        ) = 0;


    public: // Setup queries

        /**
           Get the matrix size.
        */
        inline typename uint64d<D>::Type get_matrix_size() {
            return matrix_size_;
        }

        /**
           Get the oversampled matrix size.
        */
        inline typename uint64d<D>::Type get_matrix_size_os() {
            return matrix_size_os_;
        }

        /**
           Get the convolution kernel size
        */
        inline REAL get_W() {
            return width_;
        }


        virtual ~NFFT_plan() = default;


    protected:

        // Dedicated computes
        void compute_NFFT_C2NC(ARRAY <complext<REAL>>&in, ARRAY <complext<REAL>>& out);

        void compute_NFFT_NC2C(const ARRAY <complext<REAL>>& in, ARRAY <complext<REAL>>& out);

        void compute_NFFTH_NC2C(const ARRAY <complext<REAL>>& in, ARRAY <complext<REAL>>& out);

        void compute_NFFTH_C2NC(ARRAY <complext<REAL>>& in, ARRAY <complext<REAL>>& out);

        typename uint64d<D>::Type matrix_size_;
        typename uint64d<D>::Type matrix_size_os_;
        REAL width_;
        size_t number_of_frames;
        size_t number_of_samples;

        std::unique_ptr<GriddingConvolutionBase<ARRAY, complext<REAL>, D, KaiserKernel>> conv_;
    };




}
#include "NFFT.hpp"