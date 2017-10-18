/**
    \brief CPU implementation of the non-cartesian FFT

    Comparisions were made to the gridkb function provided in the 
    Stanford Medical Image Reconstruction course lecture notes and 
    the Cuda version of the NFFT (cuNFFT)

    Uses the Kaiser Bessel for convolution
*/

#pragma once 

#include "hoNDArray.h"
#include "vector_td.h"
#include <complex>

#include <boost/shared_ptr.hpp>

#include "cpunfft_export.h"

namespace Gadgetron{

    /**
        NFFT class declaration
        ----------------------

        Real: desired precision : float or double only 
        D: dimensionality: 1D, 2D, and 3D supported
    */

    template<class Real, unsigned int D>
    class EXPORTCPUNFFT hoNFFT_plan
    {
        typedef std::complex<Real> ComplexType;

        /**
            Main interface
        */

        public:

            /**
                Default constructor
            */

            hoNFFT_plan();

            /**
                Constructor defining the required NFFT parameters

                Note: Dimensions of n must be equal.

                /param n: the non-oversampled matrix size to use for NFFT
                /param osf: the oversampling factor
                /param wg: the width of the oversampled kernel
            */

            hoNFFT_plan(
                typename uint64d<D>::Type n,
                Real osf,
                Real wg
            );

            /** 
                Destructor
            */

            ~hoNFFT_plan();

            /** 
                Perform NFFT preprocessing for a given trajectory

                \param k: the NFFT non cartesian trajectory
                \param mode: enum specifying the preprocessing mode
            */

            void preprocess(
                const hoNDArray<typename reald<Real, D>::Type>& k
            );

            /**
                Enum defining the desired NFFT operation
            */

            enum NFFT_comp_mode{
                NFFT_FORWARDS_C2NC, /** forwards cartesian to non */
                NFFT_BACKWARDS_NC2C, /** backwards non to cartesian */
                NFFT_BACKWARDS_C2NC, /** backwards cartesian to non */
                NFFT_FORWARDS_NC2C /** forwards non to cartesian */
            };

            /**
                Execute the NFFT

                \param d: the input data array
                \param m: the output matrix
                \param w: optional density compensation if not iterative 
                    provide a 0x0 if non density compensation
                \param mode: enum specifyiing the mode of operation
            */

            void compute(
                hoNDArray<ComplexType> &d,
                hoNDArray<ComplexType> &m,
                hoNDArray<Real>& w,
                NFFT_comp_mode mode
            );

            /**
                To be used by an operator for iterative reconstruction 

                \param in: the input data
                \param out: the data after MH_H has been applied

                Note: dimensions of in and out should be the same
            */
            void mult_MH_M(
                hoNDArray<ComplexType> &in,
                hoNDArray<ComplexType> &out
            );

        /**
            Utilities
        */

        public:

            /**
                Enum specifying the direction of the NFFT convolution 
            */

            enum NFFT_conv_mode{
                NFFT_CONV_C2NC, /** cartesian to non convolution */
                NFFT_CONV_NC2C /** non to cartesian convolution */
            };

            /** 
                Perform standalone convolution

                \param d: input array
                \param m: output array
                \param mode: enum specifying the mode of the convolution
            */

            void convolve(
                hoNDArray<ComplexType> &d,
                hoNDArray<ComplexType> &m,
                NFFT_conv_mode mode
            );

            /** 
                enum specifying the direction of the NFFT fft
            */

            enum NFFT_fft_mode{
                NFFT_FORWARDS,
                NFFT_BACKWARDS
            };

            /**
                Cartesian fft. Making use of the hoNDFFT class.

                \param d: input and output for he fft 
                \param mode: enum specifying the mode of the fft 
            */

            void fft(
                hoNDArray<ComplexType> &d,
                NFFT_fft_mode mode
            );

            /**
                NFFT deapodization

                \param d: input and output image to be deapodized 
                \param fourierDomain: has data been ffted
            */

            void deapodize(
                hoNDArray<ComplexType> &d,
                bool fourierDomain = false
            );

        /**
            Private implementation methods
        */

        private:

            /**
                Initialize variables and compute tables
            */

            void initialize();

            /**
                Dedicated convolutions

                The two methods below are entirely symmetric in 
                thier implementation. They could probably be
                combined for conciseness.
            */

            void convolve_NFFT_C2NC(
                hoNDArray<ComplexType> &d,
                hoNDArray<ComplexType> &m
            );

            void convolve_NFFT_NC2C(
                hoNDArray<ComplexType> &d,
                hoNDArray<ComplexType> &m
            );

            /**
                Bessel function
            */

            Real bessi0(Real x);

        /** 
            Implementation variables
        */

        private:

            typename uint64d<D>::Type n;

            Real wg, kw, kosf, kwidth, beta, osf;

            hoNDArray<Real> p, daf, nx, ny, nz;

            hoNDArray<ComplexType> da;

            hoNDArray<typename reald<Real, D>::Type> k;

    };

}
