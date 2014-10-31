/** \file hoNDFFT.h
    \brief Wrappers for FFTW for ndarrays of type std::complex.
*/

#ifndef hoNDFFT_H
#define hoNDFFT_H

#include "hoNDArray.h"
#include "cpufft_export.h"

#include <boost/thread/mutex.hpp>
#include <iostream>
#include <fftw3.h>
#include <complex>

#ifdef USE_OMP
    #include "omp.h"
#endif // USE_OMP

namespace Gadgetron{

    /** 
    Generic class for Fast Fourier Transforms using FFTW on the hoNDArray class.
    This class is a singleton because the planning and memory allocation routines of FFTW are NOT threadsafe.
    The class' template type is a REAL, ie. float or double.

    Access using e.g.
    FFT<float>::instance()
    */
    template <typename T> class EXPORTCPUFFT hoNDFFT
    {
    public:

        typedef std::complex<T> ComplexType;

        static hoNDFFT<T>* instance(); 

        void fft(hoNDArray< ComplexType >* input, unsigned int dim_to_transform)
        {
            //-1 refers to the sign of the transform, -1 for FFTW_FORWARD
            fft_int(input,dim_to_transform,-1);
        }

        void ifft(hoNDArray< ComplexType >* input, unsigned int dim_to_transform)
        {
            //1 refers to the sign of the transform, +1 for FFTW_BACKWARD
            fft_int(input,dim_to_transform,1);
        }

        void fft(hoNDArray< ComplexType >* input)
        {
            for (size_t i = 0; i < input->get_number_of_dimensions(); i++) {
                //-1 refers to the sign of the transform, -1 for FFTW_FORWARD
                fft_int(input,i,-1);
            }
        }

        void ifft(hoNDArray< ComplexType >* input)
        {
            for (size_t i = 0; i < input->get_number_of_dimensions(); i++) {
                //1 refers to the sign of the transform, +1 for FFTW_BACKWARD
                fft_int(input,i,1);
            }
        }


        void fft(hoNDArray< complext<T> >* input, unsigned int dim_to_transform)
        {
            fft((hoNDArray<ComplexType>*) input, dim_to_transform);
        }

        void ifft(hoNDArray< complext<T> >* input, unsigned int dim_to_transform)
        {
            ifft((hoNDArray<ComplexType>*) input, dim_to_transform);
        }

        void fft(hoNDArray< complext<T> >* input)
        {
            fft((hoNDArray<ComplexType>*) input);
        }

        void ifft(hoNDArray< complext<T> >* input)
        {
        	ifft((hoNDArray<ComplexType>*) input);
        }


        // 1D
        bool fftshift1D(hoNDArray< ComplexType >& a);
        bool fftshift1D(const hoNDArray< ComplexType >& a, hoNDArray< ComplexType >& r);

        bool ifftshift1D(hoNDArray< ComplexType >& a);
        bool ifftshift1D(const hoNDArray< ComplexType >& a, hoNDArray< ComplexType >& r);

        // 2D
        bool fftshift2D(hoNDArray< ComplexType >& a);
        bool fftshift2D(const hoNDArray< ComplexType >& a, hoNDArray< ComplexType >& r);

        bool ifftshift2D(hoNDArray< ComplexType >& a);
        bool ifftshift2D(const hoNDArray< ComplexType >& a, hoNDArray< ComplexType >& r);

        // 3D
        bool fftshift3D(hoNDArray< ComplexType >& a);
        bool fftshift3D(const hoNDArray< ComplexType >& a, hoNDArray< ComplexType >& r);

        bool ifftshift3D(hoNDArray< ComplexType >& a);
        bool ifftshift3D(const hoNDArray< ComplexType >& a, hoNDArray< ComplexType >& r);

        // 1D fft, in-place and out-of-place
        // the first dimension will be transformed
        bool fft1(hoNDArray< ComplexType >& a);
        bool ifft1(hoNDArray< ComplexType >& a);

        bool fft1(const hoNDArray< ComplexType >& a, hoNDArray< ComplexType >& r);
        bool ifft1(const hoNDArray< ComplexType >& a, hoNDArray< ComplexType >& r);

        // centered 1D fft
        bool fft1c(hoNDArray< ComplexType >& a);
        bool ifft1c(hoNDArray< ComplexType >& a);

        bool fft1c(const hoNDArray< ComplexType >& a, hoNDArray< ComplexType >& r);
        bool ifft1c(const hoNDArray< ComplexType >& a, hoNDArray< ComplexType >& r);

        bool fft1c(const hoNDArray< ComplexType >& a, hoNDArray< ComplexType >& r, hoNDArray< ComplexType >& buf);
        bool ifft1c(const hoNDArray< ComplexType >& a, hoNDArray< ComplexType >& r, hoNDArray< ComplexType >& buf);

        // 2D fft, in-place and out-of-place
        // the first and second dimensions will be transformed
        bool fft2(hoNDArray< ComplexType >& a);
        bool ifft2(hoNDArray< ComplexType >& a);

        bool fft2(const hoNDArray< ComplexType >& a, hoNDArray< ComplexType >& r);
        bool ifft2(const hoNDArray< ComplexType >& a, hoNDArray< ComplexType >& r);

        // centered 2D fft
        bool fft2c(hoNDArray< ComplexType >& a);
        bool ifft2c(hoNDArray< ComplexType >& a);

        bool fft2c(const hoNDArray< ComplexType >& a, hoNDArray< ComplexType >& r);
        bool ifft2c(const hoNDArray< ComplexType >& a, hoNDArray< ComplexType >& r);

        bool fft2c(const hoNDArray< ComplexType >& a, hoNDArray< ComplexType >& r, hoNDArray< ComplexType >& buf);
        bool ifft2c(const hoNDArray< ComplexType >& a, hoNDArray< ComplexType >& r, hoNDArray< ComplexType >& buf);

        // 3D fft, in-place and out-of-place
        // the first, second and third dimensions will be transformed
        bool fft3(hoNDArray< ComplexType >& a);
        bool ifft3(hoNDArray< ComplexType >& a);

        bool fft3(const hoNDArray< ComplexType >& a, hoNDArray< ComplexType >& r);
        bool ifft3(const hoNDArray< ComplexType >& a, hoNDArray< ComplexType >& r);

        // centered 3D fft
        bool fft3c(hoNDArray< ComplexType >& a);
        bool ifft3c(hoNDArray< ComplexType >& a);

        bool fft3c(const hoNDArray< ComplexType >& a, hoNDArray< ComplexType >& r);
        bool ifft3c(const hoNDArray< ComplexType >& a, hoNDArray< ComplexType >& r);

        bool fft3c(const hoNDArray< ComplexType >& a, hoNDArray< ComplexType >& r, hoNDArray< ComplexType >& buf);
        bool ifft3c(const hoNDArray< ComplexType >& a, hoNDArray< ComplexType >& r, hoNDArray< ComplexType >& buf);

    protected:

        //We are making these protected since this class is a singleton

        hoNDFFT() {
            set_function_pointers();

#ifdef USE_OMP
            num_of_max_threads_ = omp_get_num_procs();
#else
            num_of_max_threads_ = 1;
#endif // USE_OMP
        }

        virtual ~hoNDFFT() { fftw_cleanup_ptr_(); }

        void fft_int(hoNDArray< ComplexType >* input, size_t dim_to_transform, int sign);

        void set_function_pointers();

        int   (*fftw_import_wisdom_from_file_ptr_)(FILE*);
        void  (*fftw_export_wisdom_to_file_ptr_)(FILE*);
        void  (*fftw_cleanup_ptr_)(void);
        void* (*fftw_malloc_ptr_)(size_t);
        void  (*fftw_free_ptr_)(void* p);
        void  (*fftw_execute_ptr_)(void*);
        void* (*fftw_plan_dft_1d_ptr_)(int, void*, void*, int, unsigned);
        void  (*fftw_destroy_plan_ptr_)(void*);

        static hoNDFFT<T>* instance_;
        boost::mutex mutex_;

        int num_of_max_threads_;

        // the fft and ifft shift pivot for a certain length
        // [0 .. pivot-1] will be shifted to the right end
        size_t fftshiftPivot(size_t len);
        size_t ifftshiftPivot(size_t len);

        // 1D
        bool fftshift1D(const ComplexType* a, ComplexType* r, size_t x, size_t pivot);
        bool ifftshift1D(const ComplexType* a, ComplexType* r, size_t x, size_t pivot);

        bool fftshiftPivot1D(ComplexType* a, size_t x, size_t n, size_t pivot);
        bool fftshiftPivot1D(const ComplexType* a, ComplexType* r, size_t x, size_t n, size_t pivot);

        // 2D
        bool fftshiftPivot2D(const ComplexType* a, ComplexType* r, size_t x, size_t y, size_t n, size_t pivotx, size_t pivoty);
        bool fftshiftPivot2D(ComplexType* a, size_t x, size_t y, size_t n, size_t pivotx, size_t pivoty);

        bool fftshift2D(const ComplexType* a, ComplexType* r, size_t x, size_t y, size_t n);
        bool ifftshift2D(const ComplexType* a, ComplexType* r, size_t x, size_t y, size_t n);

        bool fftshift2D(ComplexType* a, size_t x, size_t y, size_t n);
        bool ifftshift2D(ComplexType* a, size_t x, size_t y, size_t n);

        // 3D
        bool fftshiftPivot3D(const ComplexType* a, ComplexType* r, size_t x, size_t y, size_t z, size_t n, size_t pivotx, size_t pivoty, size_t pivotz);
        bool fftshiftPivot3D(ComplexType* a, size_t x, size_t y, size_t z, size_t n, size_t pivotx, size_t pivoty, size_t pivotz);

        bool fftshift3D(const ComplexType* a, ComplexType* r, size_t x, size_t y, size_t z, size_t n);
        bool ifftshift3D(const ComplexType* a, ComplexType* r, size_t x, size_t y, size_t z, size_t n);

        bool fftshift3D(ComplexType* a, size_t x, size_t y, size_t z, size_t n);
        bool ifftshift3D(ComplexType* a, size_t x, size_t y, size_t z, size_t n);

        // forward: true, fft; false, inverse fft
        bool fft1(hoNDArray< ComplexType >& a, bool forward);
        bool fft2(hoNDArray< ComplexType >& a, bool forward);
        bool fft3(hoNDArray< ComplexType >& a, bool forward);

        bool fft1(hoNDArray< ComplexType >& a, hoNDArray< ComplexType >& r, bool forward);
        bool fft2(hoNDArray< ComplexType >& a, hoNDArray< ComplexType >& r, bool forward);
        bool fft3(hoNDArray< ComplexType >& a, hoNDArray< ComplexType >& r, bool forward);

        // get the number of threads used for fft
        int get_num_threads_fft1(size_t n0, size_t num);
        int get_num_threads_fft2(size_t n0, size_t n1, size_t num);
        int get_num_threads_fft3(size_t n0, size_t n1, size_t n2, size_t num);
    };
}

#endif //hoNDFFT_H
