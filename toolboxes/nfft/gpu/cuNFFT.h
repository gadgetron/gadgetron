/** \file cuNFFT.h
    \brief Cuda implementation of the non-Cartesian FFT

    Reference information on the CUDA/GPU implementation of the NFFT can be found in the papers
    
    Accelerating the Non-equispaced Fast Fourier Transform on Commodity Graphics Hardware.
    T.S. Sørensen, T. Schaeffter, K.Ø. Noe, M.S. Hansen. 
    IEEE Transactions on Medical Imaging 2008; 27(4):538-547.
    
    Real-time Reconstruction of Sensitivity Encoded Radial Magnetic Resonance Imaging Using a Graphics Processing Unit.
    T.S. Sørensen, D. Atkinson, T. Schaeffter, M.S. Hansen.
    IEEE Transactions on Medical Imaging 2009; 28(12):1974-1985. 
*/

#pragma once

#include <boost/shared_ptr.hpp>
#include <thrust/device_vector.h>

#include "cuNDArray.h"
#include "vector_td.h"
#include "complext.h"
#include "cuSparseMatrix.h"

#include "NFFT.h"
#include "cuGriddingConvolution.h"


namespace Gadgetron
{

    template <class REAL, unsigned int D>
    using cuNFFT_plan = NFFT_plan<cuNDArray,REAL,D>;


    /** \class cuNFFT_impl
        \brief Cuda implementation of the non-Cartesian FFT

        ------------------------------
        --- NFFT class declaration ---
        ------------------------------
        REAL:  desired precision : float or double
        D:  dimensionality : { 1,2,3,4 }
        ATOMICS: use atomic device memory transactions : { true, false }

        For the tested hardware the implementation using atomic operations is slower as its non-atomic counterpart.
        However, using atomic operations has the advantage of not requiring any pre-processing.
        As the preprocessing step can be quite costly in terms of memory usage,
        the atomic mode can be necessary for very large images or for 3D/4D volumes.
        Notice: currently no devices support atomics operations in double precision.
    */
    template<class REAL, unsigned int D, ConvolutionType CONV = ConvolutionType::STANDARD>
    class cuNFFT_impl : public cuNFFT_plan<REAL, D> {

    public: // Main interface

        /**
           Constructor defining the required NFFT parameters.
           \param matrix_size the matrix size to use for the NFFT. Define as a multiple of 32.
           \param matrix_size_os intermediate oversampled matrix size. Define as a multiple of 32.
           The ratio between matrix_size_os and matrix_size define the oversampling ratio for the NFFT implementation.
           Use an oversampling ratio between 1 and 2. The higher ratio the better quality results,
           however at the cost of increased execution times.
           \param W the concolution window size used in the NFFT implementation.
           The larger W the better quality at the cost of increased runtime.
           \param device the device (GPU id) to use for the NFFT computation.
           The default value of -1 indicates that the currently active device is used.
        */
        cuNFFT_impl(const vector_td<size_t,D>& matrix_size, const vector_td<size_t,D>& matrix_size_os,
                    REAL W, int device = -1);

        /**
           Destructor
        */
        virtual ~cuNFFT_impl();

        /**
           Cartesian FFT. For completeness, just invokes the cuNDFFT class.
           \param[in,out] data the data for the inplace FFT.
           \param mode enum class specifying the direction of the FFT.
           \param do_scale boolean specifying whether FFT normalization is desired.
        */
        virtual void fft(cuNDArray <complext<REAL>> &data, NFFT_fft_mode mode, bool do_scale = true) override;

        /**
           NFFT deapodization.
           \param[in,out] image the image to be deapodized (inplace).
        */
        virtual void deapodize(cuNDArray <complext<REAL>> &image, bool fourier_domain = false);

    private:

        void initialize(int device);

        boost::shared_ptr<cuNDArray<complext<REAL>>> compute_deapodization_filter(bool FFTed = false);

        // Inverse fourier transformed deapodization filter.
        boost::shared_ptr<cuNDArray<complext<REAL>>> deapodization_filter; 

        // Fourier transformed deapodization filter.
        boost::shared_ptr<cuNDArray<complext<REAL>>> deapodization_filterFFT; 

        int device_;
    };

    // Pure virtual class to cause compile errors if you try to use NFFT with double and atomics
    // - since this is not supported on the device
    template<unsigned int D>
    class cuNFFT_impl<double, D, ConvolutionType::ATOMIC> {
        virtual void atomics_not_supported_for_type_double() = 0;
    };

    template<class REAL, unsigned int D>
    struct NFFT<cuNDArray,REAL,D> 
    {
        static boost::shared_ptr<cuNFFT_plan<REAL,D>> make_plan(
            const vector_td<size_t,D>& matrix_size,
            const vector_td<size_t,D>& matrix_size_os,
            REAL W,
            ConvolutionType conv = ConvolutionType::STANDARD);
    };

    template<unsigned int D>
    struct NFFT<cuNDArray,double,D>
    {
        static boost::shared_ptr<cuNFFT_plan<double,D>> make_plan(
            const vector_td<size_t,D>& matrix_size,
            const vector_td<size_t,D>& matrix_size_os,
            double W,
            ConvolutionType conv = ConvolutionType::STANDARD);
    };
}
