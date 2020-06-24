/** \file hoNDFFT.h
    \brief Wrappers for FFTW for ndarrays of type std::complex.
*/

#ifndef hoNDFFT_H
#define hoNDFFT_H

#include "cpufft_export.h"
#include "hoNDArray.h"

#include "complext.h"
#include <complex>
#include <fftw3.h>
#include <iostream>
#include <mutex>

#ifdef USE_OMP
#include "omp.h"
#endif // USE_OMP

namespace Gadgetron {

  namespace FFT {

/**
         * Performs a standard in-place FFT over the specified dimensions
         * @tparam ComplexType Complex type, such as std::complex<float> or complex<double>
         * @param data Complex input data to be transformed
         * @param dimensions List of dimensions along which to perform the transform
 */
template <class ComplexType,
          class ENABLER = std::enable_if_t<is_complex_type_v<ComplexType>>>
void fft(hoNDArray<ComplexType> &data, std::vector<size_t> dimensions);

/**
         * \overload
 */
template <class ComplexType,
          class ENABLER = std::enable_if_t<is_complex_type_v<ComplexType>>>
void fft(hoNDArray<ComplexType> &data, size_t dimensions);

/**
         * Performs a standard in-place inverse FFT over the specified dimensions
         * @tparam ComplexType Complex type, such as std::complex<float> or complex<double>
         * @param data Complex input data to be transformed
         * @param dimensions List of dimensions along which to perform the transform
 */

template <class ComplexType,
          class ENABLER = std::enable_if_t<is_complex_type_v<ComplexType>>>
void ifft(hoNDArray<ComplexType> &data, std::vector<size_t> dimensions);

/**
         * \overload
 */
template <class ComplexType,
          class ENABLER = std::enable_if_t<is_complex_type_v<ComplexType>>>
void ifft(hoNDArray<ComplexType> &data, size_t dimensions);

/**
     * Performs a centered FFT over the 1st dimension
     * @tparam ComplexType Complex type, such as std::complex<float> or complex<double>
     * @param data Complex input data to be transformed
 */
template <class ComplexType,
          class ENABLER = std::enable_if_t<is_complex_type_v<ComplexType>>>
hoNDArray<ComplexType> fft1c(const hoNDArray<ComplexType> &data);
/**
      * Performs a centered FFT over the 1st and 2nd dimensions
      * @tparam ComplexType Complex type, such as std::complex<float> or complex<double>
      * @param data Complex input data to be transformed
 */
template <class ComplexType,
          class ENABLER = std::enable_if_t<is_complex_type_v<ComplexType>>>
hoNDArray<ComplexType> fft2c(const hoNDArray<ComplexType> &data);
/**
  * Performs a centered FFT over the 1st, 2nd and 3rd dimension
  * @tparam ComplexType Complex type, such as std::complex<float> or complex<double>
  * @param data Complex input data to be transformed
 */
template <class ComplexType,
          class ENABLER = std::enable_if_t<is_complex_type_v<ComplexType>>>
hoNDArray<ComplexType> fft3c(const hoNDArray<ComplexType> &data);

/**
    * Performs a centered FFT over the 1st dimension
    * @tparam ComplexType Complex type, such as std::complex<float> or complex<double>
    * @param data Complex input data to be transformed
 */
template <class ComplexType,
          class ENABLER = std::enable_if_t<is_complex_type_v<ComplexType>>>
hoNDArray<ComplexType> ifft1c(const hoNDArray<ComplexType> &data);
/**
      * Performs a centered FFT over the 1st and 2nd dimensions
      * @tparam ComplexType Complex type, such as std::complex<float> or complex<double>
      * @param data Complex input data to be transformed
 */
template <class ComplexType,
          class ENABLER = std::enable_if_t<is_complex_type_v<ComplexType>>>
hoNDArray<ComplexType> ifft2c(const hoNDArray<ComplexType> &data);
/**
  * Performs a centered FFT over the 1st, 2nd and 3rd dimension
  * @tparam ComplexType Complex type, such as std::complex<float> or complex<double>
  * @param data Complex input data to be transformed
 */
template <class ComplexType,
          class ENABLER = std::enable_if_t<is_complex_type_v<ComplexType>>>
hoNDArray<ComplexType> ifft3c(const hoNDArray<ComplexType> &data);

}


    /**
Generic class for Fast Fourier Transforms using FFTW on the hoNDArray class.
This class is a singleton because the planning and memory allocation routines of FFTW are NOT threadsafe.
The class' template type is a REAL, ie. float or double.

            Note that scaling is 1/sqrt(N) fir both FFT and IFFT, where N is the number of elements along the FFT
dimensions Access using e.g. FFT<float>::instance()
*/
    template <typename T> class EXPORTCPUFFT hoNDFFT {
    public:
        typedef std::complex<T> ComplexType;

        static hoNDFFT<T>* instance();

        void fft(hoNDArray<ComplexType>* input, unsigned int dim_to_transform);
        void ifft(hoNDArray<ComplexType>* input, unsigned int dim_to_transform);

        void fft(hoNDArray<ComplexType>* input);

        void ifft(hoNDArray<ComplexType>* input);

        void fft(hoNDArray<complext<T>>* input, unsigned int dim_to_transform);

        void ifft(hoNDArray<complext<T>>* input, unsigned int dim_to_transform);

        void fft(hoNDArray<complext<T>>* input);

        void ifft(hoNDArray<complext<T>>* input);

        // 1D
        void fftshift1D(hoNDArray<ComplexType>& a);
        void fftshift1D(const hoNDArray<ComplexType>& a, hoNDArray<ComplexType>& r);

        void ifftshift1D(hoNDArray<ComplexType>& a);
        void ifftshift1D(const hoNDArray<ComplexType>& a, hoNDArray<ComplexType>& r);

        // 2D
        void fftshift2D(hoNDArray<ComplexType>& a);
        void fftshift2D(const hoNDArray<ComplexType>& a, hoNDArray<ComplexType>& r);

        void ifftshift2D(hoNDArray<ComplexType>& a);
        void ifftshift2D(const hoNDArray<ComplexType>& a, hoNDArray<ComplexType>& r);

        // 3D
        void fftshift3D(hoNDArray<ComplexType>& a);
        void fftshift3D(const hoNDArray<ComplexType>& a, hoNDArray<ComplexType>& r);

        void ifftshift3D(hoNDArray<ComplexType>& a);
        void ifftshift3D(const hoNDArray<ComplexType>& a, hoNDArray<ComplexType>& r);

        // 1D fft, in-place and out-of-place
        // the first dimension will be transformed
        void fft1(hoNDArray<ComplexType>& a);
        void ifft1(hoNDArray<ComplexType>& a);

        void fft1(const hoNDArray<ComplexType>& a, hoNDArray<ComplexType>& r);
        void ifft1(const hoNDArray<ComplexType>& a, hoNDArray<ComplexType>& r);

        // centered 1D fft
        void fft1c(hoNDArray<ComplexType>& a);
        void ifft1c(hoNDArray<ComplexType>& a);

        void fft1c(const hoNDArray<ComplexType>& a, hoNDArray<ComplexType>& r);
        void ifft1c(const hoNDArray<ComplexType>& a, hoNDArray<ComplexType>& r);

        void fft1c(const hoNDArray<ComplexType>& a, hoNDArray<ComplexType>& r, hoNDArray<ComplexType>& buf);
        void ifft1c(const hoNDArray<ComplexType>& a, hoNDArray<ComplexType>& r, hoNDArray<ComplexType>& buf);

        // 2D fft, in-place and out-of-place
        // the first and second dimensions will be transformed
        void fft2(hoNDArray<ComplexType>& a);
        void ifft2(hoNDArray<ComplexType>& a);

        void fft2(const hoNDArray<ComplexType>& a, hoNDArray<ComplexType>& r);
        void ifft2(const hoNDArray<ComplexType>& a, hoNDArray<ComplexType>& r);

        // centered 2D fft
        void fft2c(hoNDArray<ComplexType>& a);
        void ifft2c(hoNDArray<ComplexType>& a);

        void fft2c(const hoNDArray<ComplexType>& a, hoNDArray<ComplexType>& r);
        void ifft2c(const hoNDArray<ComplexType>& a, hoNDArray<ComplexType>& r);

        void fft2c(const hoNDArray<ComplexType>& a, hoNDArray<ComplexType>& r, hoNDArray<ComplexType>& buf);
        void ifft2c(const hoNDArray<ComplexType>& a, hoNDArray<ComplexType>& r, hoNDArray<ComplexType>& buf);

        // 3D fft, in-place and out-of-place
        // the first, second and third dimensions will be transformed
        void fft3(hoNDArray<ComplexType>& a);
        void ifft3(hoNDArray<ComplexType>& a);

        void fft3(const hoNDArray<ComplexType>& a, hoNDArray<ComplexType>& r);
        void ifft3(const hoNDArray<ComplexType>& a, hoNDArray<ComplexType>& r);

        // centered 3D fft
        void fft3c(hoNDArray<ComplexType>& a);
        void ifft3c(hoNDArray<ComplexType>& a);

        void fft3c(const hoNDArray<ComplexType>& a, hoNDArray<ComplexType>& r);
        void ifft3c(const hoNDArray<ComplexType>& a, hoNDArray<ComplexType>& r);

        void fft3c(const hoNDArray<ComplexType>& a, hoNDArray<ComplexType>& r, hoNDArray<ComplexType>& buf);
        void ifft3c(const hoNDArray<ComplexType>& a, hoNDArray<ComplexType>& r, hoNDArray<ComplexType>& buf);

    protected:
        // We are making these protected since this class is a singleton

        hoNDFFT() = default;
        static hoNDFFT* instance_;
    };

    namespace FFT {


    }

}

#endif // hoNDFFT_H
