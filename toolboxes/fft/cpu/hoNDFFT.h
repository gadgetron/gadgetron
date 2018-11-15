/** \file hoNDFFT.h
    \brief Wrappers for FFTW for ndarrays of type std::complex.
*/

#ifndef hoNDFFT_H
#define hoNDFFT_H

#include "hoNDArray.h"
#include "cpufft_export.h"

#include <mutex>
#include <iostream>
#include <fftw3.h>
#include <complex>

#ifdef USE_OMP
    #include "omp.h"
#endif // USE_OMP

namespace Gadgetron{

template<class T> struct fftw_types{};

template<> struct fftw_types<float>{
	typedef fftwf_complex complex;
	typedef fftwf_plan_s plan;
};

template<> struct fftw_types<double>{
	typedef fftw_complex complex;
	typedef fftw_plan_s plan;
};


    /** 
    Generic class for Fast Fourier Transforms using FFTW on the hoNDArray class.
    This class is a singleton because the planning and memory allocation routines of FFTW are NOT threadsafe.
    The class' template type is a REAL, ie. float or double.

		Note that scaling is 1/sqrt(N) fir both FFT and IFFT, where N is the number of elements along the FFT dimensions
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
        void fftshift1D(hoNDArray< ComplexType >& a);
        void fftshift1D(const hoNDArray< ComplexType >& a, hoNDArray< ComplexType >& r);

        void ifftshift1D(hoNDArray< ComplexType >& a);
        void ifftshift1D(const hoNDArray< ComplexType >& a, hoNDArray< ComplexType >& r);

        // 2D
        void fftshift2D(hoNDArray< ComplexType >& a);
        void fftshift2D(const hoNDArray< ComplexType >& a, hoNDArray< ComplexType >& r);

        void ifftshift2D(hoNDArray< ComplexType >& a);
        void ifftshift2D(const hoNDArray< ComplexType >& a, hoNDArray< ComplexType >& r);

        // 3D
        void fftshift3D(hoNDArray< ComplexType >& a);
        void fftshift3D(const hoNDArray< ComplexType >& a, hoNDArray< ComplexType >& r);

        void ifftshift3D(hoNDArray< ComplexType >& a);
        void ifftshift3D(const hoNDArray< ComplexType >& a, hoNDArray< ComplexType >& r);

        // 1D fft, in-place and out-of-place
        // the first dimension will be transformed
        void fft1(hoNDArray< ComplexType >& a);
        void ifft1(hoNDArray< ComplexType >& a);

        void fft1(const hoNDArray< ComplexType >& a, hoNDArray< ComplexType >& r);
        void ifft1(const hoNDArray< ComplexType >& a, hoNDArray< ComplexType >& r);

        // centered 1D fft
        void fft1c(hoNDArray< ComplexType >& a);
        void ifft1c(hoNDArray< ComplexType >& a);

        void fft1c(const hoNDArray< ComplexType >& a, hoNDArray< ComplexType >& r);
        void ifft1c(const hoNDArray< ComplexType >& a, hoNDArray< ComplexType >& r);

        void fft1c(const hoNDArray< ComplexType >& a, hoNDArray< ComplexType >& r, hoNDArray< ComplexType >& buf);
        void ifft1c(const hoNDArray< ComplexType >& a, hoNDArray< ComplexType >& r, hoNDArray< ComplexType >& buf);

        // 2D fft, in-place and out-of-place
        // the first and second dimensions will be transformed
        void fft2(hoNDArray< ComplexType >& a);
        void ifft2(hoNDArray< ComplexType >& a);

        void fft2(const hoNDArray< ComplexType >& a, hoNDArray< ComplexType >& r);
        void ifft2(const hoNDArray< ComplexType >& a, hoNDArray< ComplexType >& r);

        // centered 2D fft
        void fft2c(hoNDArray< ComplexType >& a);
        void ifft2c(hoNDArray< ComplexType >& a);

        void fft2c(const hoNDArray< ComplexType >& a, hoNDArray< ComplexType >& r);
        void ifft2c(const hoNDArray< ComplexType >& a, hoNDArray< ComplexType >& r);

        void fft2c(const hoNDArray< ComplexType >& a, hoNDArray< ComplexType >& r, hoNDArray< ComplexType >& buf);
        void ifft2c(const hoNDArray< ComplexType >& a, hoNDArray< ComplexType >& r, hoNDArray< ComplexType >& buf);

        // 3D fft, in-place and out-of-place
        // the first, second and third dimensions will be transformed
        void fft3(hoNDArray< ComplexType >& a);
        void ifft3(hoNDArray< ComplexType >& a);

        void fft3(const hoNDArray< ComplexType >& a, hoNDArray< ComplexType >& r);
        void ifft3(const hoNDArray< ComplexType >& a, hoNDArray< ComplexType >& r);

        // centered 3D fft
        void fft3c(hoNDArray< ComplexType >& a);
        void ifft3c(hoNDArray< ComplexType >& a);

        void fft3c(const hoNDArray< ComplexType >& a, hoNDArray< ComplexType >& r);
        void ifft3c(const hoNDArray< ComplexType >& a, hoNDArray< ComplexType >& r);

        void fft3c(const hoNDArray< ComplexType >& a, hoNDArray< ComplexType >& r, hoNDArray< ComplexType >& buf);
        void ifft3c(const hoNDArray< ComplexType >& a, hoNDArray< ComplexType >& r, hoNDArray< ComplexType >& buf);


        void timeswitch(hoNDArray<ComplexType>* a, int transform_dim);
    protected:

        //We are making these protected since this class is a singleton

        hoNDFFT() {


#ifdef USE_OMP
            num_of_max_threads_ = omp_get_num_procs();
#else
            num_of_max_threads_ = 1;
#endif // USE_OMP
        }

        virtual ~hoNDFFT() { fftw_cleanup_(); }

        void fft_int(hoNDArray< ComplexType >* input, size_t dim_to_transform, int sign);


        int   fftw_import_wisdom_from_file_(FILE*);
        void  fftw_export_wisdom_to_file_(FILE*);
        void  fftw_cleanup_();
        void* fftw_malloc_(size_t);
        void  fftw_free_(void* p);
        void  fftw_execute_dft_(typename fftw_types<T>::plan * p, ComplexType*, ComplexType*);
        void  fftw_execute_(typename fftw_types<T>::plan * p);
        void fftw_print_plan_(typename fftw_types<T>::plan *p);

        typename fftw_types<T>::plan * fftw_plan_dft_1d_(int rank, ComplexType*, ComplexType*, int, unsigned);
        typename fftw_types<T>::plan * fftw_plan_dft_2d_(int dim0,int dim1, ComplexType*, ComplexType*, int, unsigned);
        typename fftw_types<T>::plan * fftw_plan_dft_3d_(int dim0,int dim1,int dim2, ComplexType*, ComplexType*, int, unsigned);

        typename fftw_types<T>::plan * fftw_plan_many_dft_(int rank, const int *n, int howmany,
                                          ComplexType *in, const int *inembed,
                                          int istride, int idist,
                                          ComplexType *out, const int *onembed,
                                          int ostride, int odist,
                                          int sign, unsigned flags);
        typename fftw_types<T>::plan * fftw_plan_dft_(int rank, ComplexType*, ComplexType*, int, unsigned);

        void  fftw_destroy_plan_(typename fftw_types<T>::plan *);



        static hoNDFFT<T>* instance_;
        std::mutex mutex_;

        int num_of_max_threads_;

        // the fft and ifft shift pivot for a certain length
        // [0 .. pivot-1] will be shifted to the right end
        size_t fftshiftPivot(size_t len);
        size_t ifftshiftPivot(size_t len);

        // 1D
        void fftshift1D(const ComplexType* a, ComplexType* r, size_t x, size_t pivot);
        void ifftshift1D(const ComplexType* a, ComplexType* r, size_t x, size_t pivot);

        void fftshiftPivot1D(ComplexType* a, size_t x, size_t n, size_t pivot);
        void fftshiftPivot1D(const ComplexType* a, ComplexType* r, size_t x, size_t n, size_t pivot);

        // 2D
        void fftshiftPivot2D(const ComplexType* a, ComplexType* r, size_t x, size_t y, size_t n, size_t pivotx, size_t pivoty);
        void fftshiftPivot2D(ComplexType* a, size_t x, size_t y, size_t n, size_t pivotx, size_t pivoty);

        void fftshift2D(const ComplexType* a, ComplexType* r, size_t x, size_t y, size_t n);
        void ifftshift2D(const ComplexType* a, ComplexType* r, size_t x, size_t y, size_t n);

        void fftshift2D(ComplexType* a, size_t x, size_t y, size_t n);
        void ifftshift2D(ComplexType* a, size_t x, size_t y, size_t n);

        // 3D
        void fftshiftPivot3D(const ComplexType* a, ComplexType* r, size_t x, size_t y, size_t z, size_t n, size_t pivotx, size_t pivoty, size_t pivotz);
        void fftshiftPivot3D(ComplexType* a, size_t x, size_t y, size_t z, size_t n, size_t pivotx, size_t pivoty, size_t pivotz);

        void fftshift3D(const ComplexType* a, ComplexType* r, size_t x, size_t y, size_t z, size_t n);
        void ifftshift3D(const ComplexType* a, ComplexType* r, size_t x, size_t y, size_t z, size_t n);

        void fftshift3D(ComplexType* a, size_t x, size_t y, size_t z, size_t n);
        void ifftshift3D(ComplexType* a, size_t x, size_t y, size_t z, size_t n);

        // forward: true, fft; false, inverse fft
        void fft1(hoNDArray< ComplexType >& a, bool forward);
        void fft2(hoNDArray< ComplexType >& a, bool forward);
        void fft3(hoNDArray< ComplexType >& a, bool forward);

        void fft1(hoNDArray< ComplexType >& a, hoNDArray< ComplexType >& r, bool forward);
        void fft2(hoNDArray< ComplexType >& a, hoNDArray< ComplexType >& r, bool forward);
        void fft3(hoNDArray< ComplexType >& a, hoNDArray< ComplexType >& r, bool forward);

        // get the number of threads used for fft
        int get_num_threads_fft1(size_t n0, size_t num);
        int get_num_threads_fft2(size_t n0, size_t n1, size_t num);
        int get_num_threads_fft3(size_t n0, size_t n1, size_t n2, size_t num);

        /**
         * Multiplies array k'th element by -1^k, causing the result of an FFT to have the frequencies centered
         */

        void phaseshift(hoNDArray<ComplexType>* a, T phase, int transform_dim);
    };
}

#endif //hoNDFFT_H
