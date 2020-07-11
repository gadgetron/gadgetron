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
#include "complext.h"
#include <complex>
#include "NFFT.h"

#include <boost/shared_ptr.hpp>
#include "hoArmadillo.h"

#include "hoNFFT_sparseMatrix.h"

namespace Gadgetron{

    /**
        NFFT class declaration
        ----------------------

        REAL: desired precision : float or double only
        D: dimensionality: 1D, 2D, and 3D supported
    */

    template<class REAL, unsigned int D>
    class hoNFFT_plan : public NFFT_plan<hoNDArray,REAL,D>
    {
        using ComplexType = std::complex<REAL>;

        /**
            Main interface
        */

        public:





            hoNFFT_plan(
                    const vector_td<size_t,D>& matrix_size,
                    const vector_td<size_t,D>& matrix_size_os,
                    REAL W
            );

             hoNFFT_plan(
                    const vector_td<size_t,D>& matrix_size,
                    REAL oversampling_factor = 1.5f,
                    REAL W = 5.5f
            );



            /** 
                Perform NFFT preprocessing for a given trajectory

                \param k: the NFFT non cartesian trajectory
                \param mode: enum specifying the preprocessing mode
            */

            virtual void preprocess(
                const hoNDArray<vector_td<REAL, D>>& k, NFFT_prep_mode prep_mode = NFFT_prep_mode::ALL
            ) override;


            void compute(
                const hoNDArray<ComplexType > &d,
                hoNDArray<ComplexType > &m,
                const hoNDArray<REAL>* dcw,
                NFFT_comp_mode mode
            );
            virtual void compute(
                const hoNDArray<complext<REAL>> &d,
                hoNDArray<complext<REAL>> &m,
                const hoNDArray<REAL>* dcw,
                NFFT_comp_mode mode
            ) override;

            /**
                To be used by an operator for iterative reconstruction 

                \param in: the input data
                \param out: the data after MH_H has been applied

                Note: dimensions of in and out should be the same
            */
            void mult_MH_M(
                const hoNDArray<ComplexType> &in,
                hoNDArray<ComplexType> &out,
                const hoNDArray<REAL>* dcw
            );

            virtual void mult_MH_M(
                const hoNDArray<complext<REAL>> &in,
                hoNDArray<complext<REAL>> &out,
                const hoNDArray<REAL>* dcw
            ) override;

        /**
            Utilities
        */

        public:


            /** 
                Perform standalone convolution

                \param d: input array
                \param m: output array
                \param mode: enum specifying the mode of the convolution
            */

            void convolve(
                const hoNDArray<ComplexType> &d,
                hoNDArray<ComplexType> &m,
                NFFT_conv_mode mode,
                bool accumulate = false
            );

            virtual void convolve(
                const hoNDArray<complext<REAL>> &d,
                hoNDArray<complext<REAL>> &m,
                NFFT_conv_mode mode,
                bool accumulate = false
            ) override;
            /**
                Cartesian fft. Making use of the hoNDFFT class.

                \param d: input and output for he fft 
                \param mode: enum specifying the mode of the fft 
            */

            void fft(
                hoNDArray<ComplexType> &d,
                NFFT_fft_mode mode,
                bool do_scale=true
            );

            virtual void fft(
                hoNDArray<complext<REAL>> &d,
                NFFT_fft_mode mode,
                bool do_scale=true
            ) override;

            /**
                NFFT deapodization

                \param d: input and output image to be deapodized 
                \param fourierDomain: has data been ffted
            */

            void deapodize(
                hoNDArray<ComplexType> &d,
                bool fourierDomain = false
            );

            virtual void deapodize(
                hoNDArray<complext<REAL>> &d,
                bool fourierDomain = false
            ) override;

        /**
            Private implementation methods
        */

        private:

            /**
                Dedicated convolutions

                The two methods below are entirely symmetric in 
                thier implementation. They could probably be
                combined for conciseness.
            */

            void convolve_NFFT_C2NC(
                const hoNDArray<ComplexType> &d,
                hoNDArray<ComplexType> &m, bool accumulate
            );

            void convolve_NFFT_NC2C(
                const hoNDArray<ComplexType> &d,
                hoNDArray<ComplexType> &m, bool accumulate
            );


            static vector_td<REAL,D> compute_beta(REAL W, const vector_td<size_t,D>& matrix_size, const vector_td<size_t,D>& matrix_size_os);


        /** 
            Implementation variables
        */

        private:

        vector_td<REAL,D> beta;
        std::vector<NFFT_internal::NFFT_Matrix<REAL>> convolution_matrix;
        std::vector<NFFT_internal::NFFT_Matrix<REAL>> convolution_matrix_T;

        hoNDArray<ComplexType> deapodization_filter_IFFT;
        hoNDArray<ComplexType> deapodization_filter_FFT;

    };


    template<class REAL, unsigned int D> struct NFFT<hoNDArray,REAL,D>{

        using NFFT_plan = hoNFFT_plan<REAL,D>;
         static boost::shared_ptr<hoNFFT_plan<REAL,D>> make_plan(const vector_td<size_t,D>& matrix_size, const vector_td<size_t,D>& matrix_size_os,
                    REAL W);
    };

}
