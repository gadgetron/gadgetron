/** \file  hoNDMath_util.cpp
\brief math function utility
*/

#include "hoNDMath_util.h"
#include <stdexcept>
#include <iostream>
#include "GadgetronCommon.h"

/// uncomment this to disable the explicit MKL calls
// #undef USE_MKL

#ifdef USE_MKL
    #include "mkl.h"

    #ifdef MKL_ILP64
        #define ILP_MODE_ON 1
    #endif
#endif // USE_MKL

#ifndef lapack_int
    #ifdef USE_MKL
        #ifdef MKL_ILP64
            #define lapack_int __int64
        #else
            #define lapack_int int
        #endif // MKL_ILP64
    #else
        #define lapack_int int
    #endif // USE_MKL
#endif // lapack_int

#ifdef USE_MKL

    #define saxpy_ saxpy
    #define daxpy_ daxpy
    #define caxpy_ caxpy
    #define zaxpy_ zaxpy

    #define snrm2_ snrm2
    #define scnrm2_ scnrm2
    #define dnrm2_ dnrm2
    #define dznrm2_ dznrm2

    #define sasum_ sasum
    #define scasum_ scasum
    #define dasum_ dasum
    #define dzasum_ dzasum

    #define sscal_ sscal
    #define dscal_ dscal
    #define cscal_ cscal
    #define zscal_ zscal
    #define csscal_ csscal
    #define zdscal_ zdscal

    #define cdotc_ cdotc
    #define zdotc_ zdotc
    #define cdotu_ cdotu
    #define zdotu_ zdotu

    #define isamax_ isamax
    #define idamax_ idamax
    #define icamax_ icamax
    #define izamax_ izamax

#else

    #ifndef lapack_complex_float
        #define lapack_complex_float GT_Complex8
    #endif // lapack_complex_float

    #ifndef lapack_complex_double
        #define lapack_complex_double GT_Complex16
    #endif // #ifndef lapack_complex_double

    //Declaration of BLAS and LAPACK routines
    extern "C"
    {
        /// Computes a vector-scalar product and adds the result to a vector.
        /// y = a*x + y
        void saxpy_(lapack_int* N, float* a, float* x, lapack_int* incx, float* y, lapack_int* incy);
        void daxpy_(lapack_int* N, double* a, double* x, lapack_int* incx, double* y, lapack_int* incy);
        void caxpy_(lapack_int* N, lapack_complex_float* a, lapack_complex_float* x, lapack_int* incx, lapack_complex_float* y, lapack_int* incy);
        void zaxpy_(lapack_int* N, lapack_complex_double* a, lapack_complex_double* x, lapack_int* incx, lapack_complex_double* y, lapack_int* incy);

        /// Computes the Euclidean norm of a vector.
        float snrm2_(lapack_int* N, float* x, lapack_int* incx);
        float scnrm2_(lapack_int* N, lapack_complex_float* x, lapack_int* incx);
        double dnrm2_(lapack_int* N, double* x, lapack_int* incx);
        double dznrm2_(lapack_int* N, lapack_complex_double* x, lapack_int* incx);

        /// Computes the sum of magnitudes of the vector elements.
        float sasum_(lapack_int* N, float* x, lapack_int* incx);
        float scasum_(lapack_int* N, lapack_complex_float* x, lapack_int* incx);
        double dasum_(lapack_int* N, double* x, lapack_int* incx);
        double dzasum_(lapack_int* N, lapack_complex_double* x, lapack_int* incx);

        /// Computes the product of a vector by a scalar.
        void sscal_(lapack_int* N, float* a, float* x, lapack_int* incx);
        void dscal_(lapack_int* N, double* a, double* x, lapack_int* incx);
        void cscal_(lapack_int* N, lapack_complex_float* a, lapack_complex_float* x, lapack_int* incx);
        void zscal_(lapack_int* N, lapack_complex_double* a, lapack_complex_double* x, lapack_int* incx);
        void csscal_(lapack_int* N, float* a, lapack_complex_float* x, lapack_int* incx);
        void zdscal_(lapack_int* N, double* a, lapack_complex_double* x, lapack_int* incx);

        /// Computes a dot product of a conjugated vector with another vector.
        void cdotc_(void* r, lapack_int* N, lapack_complex_float* x, lapack_int* incx, lapack_complex_float* y, lapack_int* incy);
        void zdotc_(void* r, lapack_int* N, lapack_complex_double* x, lapack_int* incx, lapack_complex_double* y, lapack_int* incy);

        /// Computes a vector-vector dot product.
        void cdotu_(void* r, lapack_int* N, lapack_complex_float* x, lapack_int* incx, lapack_complex_float* y, lapack_int* incy);
        void zdotu_(void* r, lapack_int* N, lapack_complex_double* x, lapack_int* incx, lapack_complex_double* y, lapack_int* incy);

        /// Finds the index of the element with the maximal absolute value.
        lapack_int isamax_(lapack_int* N, float* x, lapack_int* incx);
        lapack_int idamax_(lapack_int* N, double* x, lapack_int* incx);
        lapack_int icamax_(lapack_int* N, lapack_complex_float* x, lapack_int* incx);
        lapack_int izamax_(lapack_int* N, lapack_complex_double* x, lapack_int* incx);
        }
#endif // USE_MKL

#define NumElementsUseThreading 128*128*6

const unsigned long long FourGBLimit = (unsigned long long)(1024.0*1024*1024*4);

namespace Gadgetron { namespace math {

    /// --------------------------------------------------------------------

    template <typename T> inline void scal_64bit_mode(size_t N, T a, T* x)
    {
        long long n;

#pragma omp parallel for private(n) if (N>NumElementsUseThreading)
        for (n = 0; n < (long long)N; n++)
        {
            const T& c = x[n];
            const typename realType<T>::Type re = c.real();
            const typename realType<T>::Type im = c.imag();

            const typename realType<T>::Type ar = a.real();
            const typename realType<T>::Type ai = a.imag();

            x[n].real(re*ar-im*ai);
            x[n].imag(re*ai+im*ar);
        }
    }

    template <typename T> inline void scal_64bit_mode(size_t N, typename realType<T>::Type a, T* x)
    {
        long long n;

#pragma omp parallel for private(n) if (N>NumElementsUseThreading)
        for (n = 0; n < (long long)N; n++)
        {
            const T& c = x[n];
            const typename realType<T>::Type re = c.real();
            const typename realType<T>::Type im = c.imag();

            x[n].real(re*a);
            x[n].imag(im*a);
        }
    }

    // -----------------------------------

    template <> EXPORTCPUCOREMATH void scal(size_t N, float a, float* x)
    {
        lapack_int num = (lapack_int)N;
        lapack_int incx = 1;

#ifdef ILP_MODE_ON
        sscal_(&num, &a, x, &incx);
#else
        if ( N < FourGBLimit )
        {
            sscal_(&num, &a, x, &incx);
        }
        else
        {
            long long n;
#pragma omp parallel for private(n) if (N>NumElementsUseThreading)
            for (n = 0; n < (long long)N; n++)
            {
                x[n] *= a;
            }
        }
#endif // ILP_MODE_ON
    }

    // -----------------------------------

    template <> EXPORTCPUCOREMATH void scal(size_t N, double a, double* x)
    {
        lapack_int num = (lapack_int)N;
        lapack_int incx = 1;

#ifdef ILP_MODE_ON
        dscal_(&num, &a, x, &incx);
#else
        if ( N < FourGBLimit )
        {
            dscal_(&num, &a, x, &incx);
        }
        else
        {
            long long n;
#pragma omp parallel for private(n) if (N>NumElementsUseThreading)
            for (n = 0; n < (long long)N; n++)
            {
                x[n] *= a;
            }
        }
#endif // ILP_MODE_ON
    }

    // -----------------------------------

    template <> EXPORTCPUCOREMATH void scal(size_t N, GT_Complex8 a, GT_Complex8* x)
    {
        lapack_int num = (lapack_int)N;
        lapack_int incx = 1;

#ifdef ILP_MODE_ON
        cscal_(&num, (lapack_complex_float*)(&a), (lapack_complex_float*)(x), &incx);
#else
        if ( N < FourGBLimit )
        {
            cscal_(&num, (lapack_complex_float*)(&a), (lapack_complex_float*)(x), &incx);
        }
        else
        {
            scal_64bit_mode(N, a, x);
        }
#endif // ILP_MODE_ON
    }

    // -----------------------------------

    template <> EXPORTCPUCOREMATH void scal(size_t N, GT_Complex16 a, GT_Complex16* x)
    {
        lapack_int num = (lapack_int)N;
        lapack_int incx = 1;

#ifdef ILP_MODE_ON
        zscal_(&num, (lapack_complex_double*)(&a), (lapack_complex_double*)(x), &incx);
#else
        if ( N < FourGBLimit )
        {
            zscal_(&num, (lapack_complex_double*)(&a), (lapack_complex_double*)(x), &incx);
        }
        else
        {
            scal_64bit_mode(N, a, x);
        }
#endif // ILP_MODE_ON
    }

    // -----------------------------------

    template <> EXPORTCPUCOREMATH void scal(size_t N, float a, GT_Complex8* x)
    {
        lapack_int num = (lapack_int)N;
        lapack_int incx = 1;

#ifdef ILP_MODE_ON
        csscal_(&num, &a, (lapack_complex_float*)(x), &incx);
#else
        if ( N < FourGBLimit )
        {
            csscal_(&num, &a, (lapack_complex_float*)(x), &incx);
        }
        else
        {
            scal_64bit_mode(N, a, x);
        }
#endif // ILP_MODE_ON
    }

    // -----------------------------------

    template <> EXPORTCPUCOREMATH void scal(size_t N, double a, GT_Complex16* x)
    {
        lapack_int num = (lapack_int)N;
        lapack_int incx = 1;

#ifdef ILP_MODE_ON
        zdscal_(&num, &a, (lapack_complex_double*)(x), &incx);
#else
        if ( N < FourGBLimit )
        {
            zdscal_(&num, &a, (lapack_complex_double*)(x), &incx);
        }
        else
        {
            scal_64bit_mode(N, a, x);
        }
#endif // ILP_MODE_ON
    }

    /// --------------------------------------------------------------------

    template <typename T> inline void axpy_64bit_mode(T a, size_t N, const T* x, const T* y, T* r)
    {
        long long n;

        #pragma omp parallel for private(n) shared(N, r, a, x, y) if(N>NumElementsUseThreading)
        for ( n=0; n<(long long)N; ++n)
        {
            const T& vx = x[n];
            const typename realType<T>::Type re1 = vx.real();
            const typename realType<T>::Type im1 = vx.imag();

            const T& vy = y[n];
            const typename realType<T>::Type re2 = vy.real();
            const typename realType<T>::Type im2 = vy.imag();

            const typename realType<T>::Type ar = a.real();
            const typename realType<T>::Type ai = a.imag();

            r[n].real(re2 + ar*re1 - ai*im1);
            r[n].imag(im2 + ar*im1 + ai*re1);
        }
    }

    template <> EXPORTCPUCOREMATH void axpy(float a, size_t N, const float* x, const float* y, float* r)
    {
        lapack_int num = (lapack_int)N;
        lapack_int incx = 1;
        lapack_int incy = 1;

        if ( y != r )
        {
            memcpy(r, y, sizeof(float)*N);
        }

#ifdef ILP_MODE_ON
        saxpy_(&num, &a, const_cast<float*>(x), &incx, r, &incy);
#else
        if ( N < FourGBLimit )
        {
            saxpy_(&num, &a, const_cast<float*>(x), &incx, r, &incy);
        }
        else
        {
            long long n;

            #pragma omp parallel for default(none) private(n) shared(N, r, a , x, y) if(N>NumElementsUseThreading)
            for ( n=0; n<(long long)N; ++n)
            {
                r[n] = a*x[n] + y[n];
            }
        }
#endif // ILP_MODE_ON
    }

    template <> EXPORTCPUCOREMATH void axpy(double a, size_t N, const double* x, const double* y, double* r)
    {
        lapack_int num = (lapack_int)N;
        lapack_int incx = 1;
        lapack_int incy = 1;

        if ( y != r )
        {
            memcpy(r, y, sizeof(double)*N);
        }

#ifdef ILP_MODE_ON
        daxpy_(&num, &a, const_cast<double*>(x), &incx, r, &incy);
#else
        if ( N < FourGBLimit )
        {
            daxpy_(&num, &a, const_cast<double*>(x), &incx, r, &incy);
        }
        else
        {
            long long n;

            #pragma omp parallel for default(none) private(n) shared(N, r, a , x, y) if(N>NumElementsUseThreading)
            for ( n=0; n<(long long)N; ++n)
            {
                r[n] = a*x[n] + y[n];
            }
        }
#endif // ILP_MODE_ON
    }

    template <> EXPORTCPUCOREMATH void axpy(GT_Complex8 a, size_t N, const GT_Complex8* x, const GT_Complex8* y, GT_Complex8* r)
    {
        lapack_int num = (lapack_int)N;
        lapack_int incx = 1;
        lapack_int incy = 1;

        if ( y != r )
        {
            memcpy(r, y, sizeof(GT_Complex8)*N);
        }

#ifdef ILP_MODE_ON
        caxpy_(&num, (lapack_complex_float*)(&a), (lapack_complex_float*)(x), &incx, (lapack_complex_float*)(r), &incy);
#else
        if ( N < FourGBLimit )
        {
            caxpy_(&num, (lapack_complex_float*)(&a), (lapack_complex_float*)(x), &incx, (lapack_complex_float*)(r), &incy);
        }
        else
        {
            axpy_64bit_mode(a, N, x, y, r);
        }
#endif // ILP_MODE_ON
    }

    template <> EXPORTCPUCOREMATH void axpy(GT_Complex16 a, size_t N, const GT_Complex16* x, const GT_Complex16* y, GT_Complex16* r)
    {
        lapack_int num = (lapack_int)N;
        lapack_int incx = 1;
        lapack_int incy = 1;

        if ( y != r )
        {
            memcpy(r, y, sizeof(GT_Complex16)*N);
        }

#ifdef ILP_MODE_ON
        zaxpy_(&num, (lapack_complex_double*)(&a), (lapack_complex_double*)(x), &incx, (lapack_complex_double*)(r), &incy);
#else
        if ( N < FourGBLimit )
        {
            zaxpy_(&num, (lapack_complex_double*)(&a), (lapack_complex_double*)(x), &incx, (lapack_complex_double*)(r), &incy);
        }
        else
        {
            axpy_64bit_mode(a, N, x, y, r);
        }
#endif // ILP_MODE_ON
    }

    /// --------------------------------------------------------------------

#ifdef USE_MKL
    template <> EXPORTCPUCOREMATH void add(size_t N, const float* x, const float* y, float* r)
    {
        vsAdd((lapack_int)N, x, y, r);
        GADGET_CHECK_THROW(vmlGetErrStatus()==0);
    }

    template <> EXPORTCPUCOREMATH void add(size_t N, const double* x, const double* y, double* r)
    {
        vdAdd((lapack_int)N, x, y, r);
        GADGET_CHECK_THROW(vmlGetErrStatus()==0);
    }

    template <> EXPORTCPUCOREMATH void add(size_t N, const GT_Complex8* x, const GT_Complex8* y, GT_Complex8* r)
    {
        vcAdd((lapack_int)N, (MKL_Complex8*)(x), (MKL_Complex8*)(y), (MKL_Complex8*)(r));
        GADGET_CHECK_THROW(vmlGetErrStatus()==0);
    }

    template <> EXPORTCPUCOREMATH void add(size_t N, const GT_Complex16* x, const GT_Complex16* y, GT_Complex16* r)
    {
        vzAdd((lapack_int)N, (MKL_Complex16*)(x), (MKL_Complex16*)(y), (MKL_Complex16*)(r));
        GADGET_CHECK_THROW(vmlGetErrStatus()==0);
    }

#else

    template <typename T> void add(size_t N, const T* x, const T* y, T* r)
    {
        T a = 1;
        if ( x == r && y!=r )
        {
            axpy(a, N, y, x, r);
        }
        else if (x != r && y==r )
        {
            axpy(a, N, x, y, r);
        }
        else if ( x==r && y==r )
        {
            scal( N, (typename realType<T>::Type)(2), r);
        }
        else
        {
            axpy(a, N, x, y, r);
        }
    }

    template EXPORTCPUCOREMATH void add(size_t N, const float* x, const float* y, float* r);
    template EXPORTCPUCOREMATH void add(size_t N, const double* x, const double* y, double* r);
    template EXPORTCPUCOREMATH void add(size_t N, const GT_Complex8* x, const GT_Complex8* y, GT_Complex8* r);
    template EXPORTCPUCOREMATH void add(size_t N, const GT_Complex16* x, const GT_Complex16* y, GT_Complex16* r);

#endif // USE_MKL

    /// --------------------------------------------------------------------

#ifdef USE_MKL
    template <> EXPORTCPUCOREMATH void subtract(size_t N, const float* x, const float* y, float* r)
    {
        vsSub((lapack_int)N, x, y, r);
        GADGET_CHECK_THROW(vmlGetErrStatus()==0);
    }

    template <> EXPORTCPUCOREMATH void subtract(size_t N, const double* x, const double* y, double* r)
    {
        vdSub((lapack_int)N, x, y, r);
        GADGET_CHECK_THROW(vmlGetErrStatus()==0);
    }

    template <> EXPORTCPUCOREMATH void subtract(size_t N, const GT_Complex8* x, const GT_Complex8* y, GT_Complex8* r)
    {
        vcSub((lapack_int)N, (MKL_Complex8*)(x), (MKL_Complex8*)(y), (MKL_Complex8*)(r));
        GADGET_CHECK_THROW(vmlGetErrStatus()==0);
    }

    template <> EXPORTCPUCOREMATH void subtract(size_t N, const GT_Complex16* x, const GT_Complex16* y, GT_Complex16* r)
    {
        vzSub((lapack_int)N, (MKL_Complex16*)(x), (MKL_Complex16*)(y), (MKL_Complex16*)(r));
        GADGET_CHECK_THROW(vmlGetErrStatus()==0);
    }

#else

    template <typename T> void subtract(size_t N, const T* x, const T* y, T* r)
    {
        long long n;
#pragma omp parallel for default(none) private(n) shared(N, x, y, r) if (N>NumElementsUseThreading)
        for ( n=0; n<(long long)N; n++ )
        {
            r[n] = x[n] - y[n];
        }
    }

    template EXPORTCPUCOREMATH void subtract(size_t N, const float* x, const float* y, float* r);
    template EXPORTCPUCOREMATH void subtract(size_t N, const double* x, const double* y, double* r);
    template EXPORTCPUCOREMATH void subtract(size_t N, const GT_Complex8* x, const GT_Complex8* y, GT_Complex8* r);
    template EXPORTCPUCOREMATH void subtract(size_t N, const GT_Complex16* x, const GT_Complex16* y, GT_Complex16* r);
#endif // USE_MKL

    /// --------------------------------------------------------------------

#ifdef USE_MKL
    template <> EXPORTCPUCOREMATH void multiply(size_t N, const float* x, const float* y, float* r)
    {
        vsMul((lapack_int)N, x, y, r);
        GADGET_CHECK_THROW(vmlGetErrStatus()==0);
    }

    template <> EXPORTCPUCOREMATH void multiply(size_t N, const double* x, const double* y, double* r)
    {
        vdMul((lapack_int)N, x, y, r);
        GADGET_CHECK_THROW(vmlGetErrStatus()==0);
    }

    template <> EXPORTCPUCOREMATH void multiply(size_t N, const GT_Complex8* x, const GT_Complex8* y, GT_Complex8* r)
    {
        vcMul((lapack_int)N, (MKL_Complex8*)(x), (MKL_Complex8*)(y), (MKL_Complex8*)(r));
        GADGET_CHECK_THROW(vmlGetErrStatus()==0);
    }

    template <> EXPORTCPUCOREMATH void multiply(size_t N, const GT_Complex16* x, const GT_Complex16* y, GT_Complex16* r)
    {
        vzMul((lapack_int)N, (MKL_Complex16*)(x), (MKL_Complex16*)(y), (MKL_Complex16*)(r));
        GADGET_CHECK_THROW(vmlGetErrStatus()==0);
    }

#else

    template <typename T> 
    void multiply(size_t N, const T* x, const T* y, T* r)
    {
        long long n;
#pragma omp parallel for private(n) if (N>NumElementsUseThreading)
        for ( n=0; n<(long long)N; n++ )
        {
            const T& a = x[n];
            const T& b = y[n];
            r[n] = a*b;
        }
    }

    template EXPORTCPUCOREMATH void multiply(size_t N, const float* x, const float* y, float* r);
    template EXPORTCPUCOREMATH void multiply(size_t N, const double* x, const double* y, double* r);

    template <> EXPORTCPUCOREMATH void multiply(size_t N, const GT_Complex8* x, const GT_Complex8* y, GT_Complex8* r)
    {
        long long i;
        #pragma omp parallel for private(i) if (N>NumElementsUseThreading)
        for (i = 0; i < (long long)N; i++)
        {
            const std::complex<float>& a1 = x[i];
            const std::complex<float>& b1 = y[i];
            const float a = a1.real();
            const float b = a1.imag();
            const float c = b1.real();
            const float d = b1.imag();

            reinterpret_cast<float(&)[2]>(r[i])[0] = a*c-b*d;
            reinterpret_cast<float(&)[2]>(r[i])[1] = a*d+b*c;
        }
    }

    template <> EXPORTCPUCOREMATH void multiply(size_t N, const GT_Complex16* x, const GT_Complex16* y, GT_Complex16* r)
    {
        long long i;
        #pragma omp parallel for private(i) if (N>NumElementsUseThreading)
        for (i = 0; i < (long long)N; i++)
        {
            const std::complex<double>& a1 = x[i];
            const std::complex<double>& b1 = y[i];
            const double a = a1.real();
            const double b = a1.imag();
            const double c = b1.real();
            const double d = b1.imag();

            reinterpret_cast<double(&)[2]>(r[i])[0] = a*c-b*d;
            reinterpret_cast<double(&)[2]>(r[i])[1] = a*d+b*c;
        }
    }

#endif // USE_MKL

    /// --------------------------------------------------------------------

#ifdef USE_MKL
    template <> EXPORTCPUCOREMATH void divide(size_t N, const float* x, const float* y, float* r)
    {
        vsDiv((lapack_int)N, x, y, r);
        GADGET_CHECK_THROW(vmlGetErrStatus()==0);
    }

    template <> EXPORTCPUCOREMATH void divide(size_t N, const double* x, const double* y, double* r)
    {
        vdDiv((lapack_int)N, x, y, r);
        GADGET_CHECK_THROW(vmlGetErrStatus()==0);
    }

    template <> EXPORTCPUCOREMATH void divide(size_t N, const GT_Complex8* x, const GT_Complex8* y, GT_Complex8* r)
    {
        vcDiv((lapack_int)N, (MKL_Complex8*)(x), (MKL_Complex8*)(y), (MKL_Complex8*)(r));
        GADGET_CHECK_THROW(vmlGetErrStatus()==0);
    }

    template <> EXPORTCPUCOREMATH void divide(size_t N, const GT_Complex16* x, const GT_Complex16* y, GT_Complex16* r)
    {
        vzDiv((lapack_int)N, (MKL_Complex16*)(x), (MKL_Complex16*)(y), (MKL_Complex16*)(r));
        GADGET_CHECK_THROW(vmlGetErrStatus()==0);
    }

#else

    template <typename T> 
    void divide(size_t N, const T* x, const T* y, T* r)
    {
        long long n;
#pragma omp parallel for private(n) if (N>NumElementsUseThreading)
        for ( n=0; n<(long long)N; n++ )
        {
            const T& a = x[n];
            const T& b = y[n];
            r[n] = a/b;
        }
    }

    template EXPORTCPUCOREMATH void divide(size_t N, const float* x, const float* y, float* r);
    template EXPORTCPUCOREMATH void divide(size_t N, const double* x, const double* y, double* r);

    template <> EXPORTCPUCOREMATH void divide(size_t N, const GT_Complex8* x, const GT_Complex8* y, GT_Complex8* r)
    {
        long long i;
        #pragma omp parallel for private(i) if (N>NumElementsUseThreading)
        for (i = 0; i < (long long)N; i++)
        {
            const std::complex<float>& a1 = x[i];
            const std::complex<float>& b1 = y[i];
            const float a = a1.real();
            const float b = a1.imag();
            const float c = b1.real();
            const float d = b1.imag();

            const float m = 1/(c*c+d*d);

            reinterpret_cast<float(&)[2]>(r[i])[0] = (a*c+b*d)*m;
            reinterpret_cast<float(&)[2]>(r[i])[1] = (b*c-a*d)*m;
        }
    }

    template <> EXPORTCPUCOREMATH void divide(size_t N, const GT_Complex16* x, const GT_Complex16* y, GT_Complex16* r)
    {
        long long i;
        #pragma omp parallel for private(i) if (N>NumElementsUseThreading)
        for (i = 0; i < (long long)N; i++)
        {
            const std::complex<double>& a1 = x[i];
            const std::complex<double>& b1 = y[i];
            const double a = a1.real();
            const double b = a1.imag();
            const double c = b1.real();
            const double d = b1.imag();

            const double m = 1/(c*c+d*d);

            reinterpret_cast<double(&)[2]>(r[i])[0] = (a*c+b*d)*m;
            reinterpret_cast<double(&)[2]>(r[i])[1] = (b*c-a*d)*m;
        }
    }
#endif // USE_MKL

    /// --------------------------------------------------------------------

#ifdef USE_MKL
    template <> EXPORTCPUCOREMATH void sqrt(size_t N, const float* x, float* r)
    {
        vsSqrt((lapack_int)N, x, r);
        GADGET_CHECK_THROW(vmlGetErrStatus()==0);
    }

    template <> EXPORTCPUCOREMATH void sqrt(size_t N, const double* x, double* r)
    {
        vdSqrt((lapack_int)N, x, r);
        GADGET_CHECK_THROW(vmlGetErrStatus()==0);
    }

    template <> EXPORTCPUCOREMATH void sqrt(size_t N, const GT_Complex8* x, GT_Complex8* r)
    {
        vcSqrt((lapack_int)N, (MKL_Complex8*)(x), (MKL_Complex8*)(r));
        GADGET_CHECK_THROW(vmlGetErrStatus()==0);
    }

    template <> EXPORTCPUCOREMATH void sqrt(size_t N, const GT_Complex16* x, GT_Complex16* r)
    {
        vzSqrt((lapack_int)N, (MKL_Complex16*)(x), (MKL_Complex16*)(r));
        GADGET_CHECK_THROW(vmlGetErrStatus()==0);
    }

#else

    template <typename T> 
    void sqrt(size_t N, const T* x, T* r)
    {
        long long n;
#pragma omp parallel for default(none) private(n) shared(N, x, r) if (N>NumElementsUseThreading)
        for ( n=0; n<(long long)N; n++ )
        {
            r[n] = std::sqrt(x[n]);
        }
    }

    template EXPORTCPUCOREMATH void sqrt(size_t N, const float* x, float* r);
    template EXPORTCPUCOREMATH void sqrt(size_t N, const double* x, double* r);
    template EXPORTCPUCOREMATH void sqrt(size_t N, const GT_Complex8* x, GT_Complex8* r);
    template EXPORTCPUCOREMATH void sqrt(size_t N, const GT_Complex16* x, GT_Complex16* r);
#endif // USE_MKL

    /// --------------------------------------------------------------------

    template <typename T> 
    void minAbsolute(size_t N, const T* x, T& r, size_t& ind)
    {
        ind = 0;
        if ( N == 0 ) return;

        long long n;

        typename realType<T>::Type v = abs(x[0]);
        typename realType<T>::Type v2;

        ind = 0;
        for ( n=1; n<(long long)N; n++ )
        {
            v2 = abs(x[n]);
            if ( v2 < v )
            {
                v = v2;
                ind = n;
            }
        }

        r = x[ind];
    }

    template EXPORTCPUCOREMATH void minAbsolute(size_t N, const float* x, float& r, size_t& ind);
    template EXPORTCPUCOREMATH void minAbsolute(size_t N, const double* x, double& r, size_t& ind);
    template EXPORTCPUCOREMATH void minAbsolute(size_t N, const GT_Complex8* x, GT_Complex8& r, size_t& ind);
    template EXPORTCPUCOREMATH void minAbsolute(size_t N, const GT_Complex16* x, GT_Complex16& r, size_t& ind);

    /// --------------------------------------------------------------------

    template <typename T> 
    void maxAbsolute(size_t N, const T* x, T& r, size_t& ind)
    {
        try
        {
            ind = 0;
            if ( N == 0 ) return;

            long long n;

            typename realType<T>::Type v = abs(x[0]);
            typename realType<T>::Type v2;

            ind = 0;
            for ( n=1; n<(long long)N; n++ )
            {
                v2 = abs(x[n]);
                if ( v2 > v )
                {
                    v = v2;
                    ind = n;
                }
            }

            r = x[ind];
        }
        catch(...)
        {
            GADGET_THROW("Error happened in maxAbsolute(size_t N, const T* x, T& r, size_t& ind) ... ");
        }
    }

    template EXPORTCPUCOREMATH void maxAbsolute(size_t N, const float* x, float& r, size_t& ind);
    template EXPORTCPUCOREMATH void maxAbsolute(size_t N, const double* x, double& r, size_t& ind);
    template EXPORTCPUCOREMATH void maxAbsolute(size_t N, const GT_Complex8* x, GT_Complex8& r, size_t& ind);
    template EXPORTCPUCOREMATH void maxAbsolute(size_t N, const GT_Complex16* x, GT_Complex16& r, size_t& ind);

    /// --------------------------------------------------------------------

#ifdef USE_MKL

    template <> EXPORTCPUCOREMATH void multiplyConj(size_t N, const GT_Complex8* x, const GT_Complex8* y, GT_Complex8* r)
    {
        vcMulByConj((lapack_int)N, (MKL_Complex8*)(x), (MKL_Complex8*)(y), (MKL_Complex8*)(r));
        GADGET_CHECK_THROW(vmlGetErrStatus()==0);
    }

    template <> EXPORTCPUCOREMATH void multiplyConj(size_t N, const GT_Complex16* x, const GT_Complex16* y, GT_Complex16* r)
    {
        vzMulByConj((lapack_int)N, (MKL_Complex16*)(x), (MKL_Complex16*)(y), (MKL_Complex16*)(r));
        GADGET_CHECK_THROW(vmlGetErrStatus()==0);
    }

#else

    template <> EXPORTCPUCOREMATH 
    void multiplyConj(size_t N, const GT_Complex8* x, const GT_Complex8* y, GT_Complex8* r)
    {
        long long n;

#pragma omp parallel for private(n) shared(N, x, y, r) if (N>NumElementsUseThreading)
        for ( n=0; n<(long long)N; n++ )
        {
            const float a = x[n].real();
            const float b = x[n].imag();
            const float c = y[n].real();
            const float d = y[n].imag();

            reinterpret_cast<float(&)[2]>(r[n])[0] = (a*c + b*d);
            reinterpret_cast<float(&)[2]>(r[n])[1] = (c*b - a*d);
        }
    }

    template <> EXPORTCPUCOREMATH 
    void multiplyConj(size_t N, const GT_Complex16* x, const GT_Complex16* y, GT_Complex16* r)
    {
        long long n;

#pragma omp parallel for private(n) shared(N, x, y, r) if (N>NumElementsUseThreading)
        for ( n=0; n<(long long)N; n++ )
        {
            const double a = x[n].real();
            const double b = x[n].imag();
            const double c = y[n].real();
            const double d = y[n].imag();

            reinterpret_cast<double(&)[2]>(r[n])[0] = (a*c + b*d);
            reinterpret_cast<double(&)[2]>(r[n])[1] = (c*b - a*d);
        }
    }

#endif // USE_MKL
    /// --------------------------------------------------------------------

#ifdef USE_MKL

    template <> EXPORTCPUCOREMATH void conjugate(size_t N, const GT_Complex8* x, GT_Complex8* r)
    {
        vcConj((lapack_int)N, (MKL_Complex8*)(x), (MKL_Complex8*)(r));
        GADGET_CHECK_THROW(vmlGetErrStatus()==0);
    }

    template <> EXPORTCPUCOREMATH void conjugate(size_t N, const GT_Complex16* x, GT_Complex16* r)
    {
        vzConj((lapack_int)N, (MKL_Complex16*)(x), (MKL_Complex16*)(r));
        GADGET_CHECK_THROW(vmlGetErrStatus()==0);
    }

#else

    template <> EXPORTCPUCOREMATH 
    void conjugate(size_t N, const GT_Complex8* x, GT_Complex8* r)
    {
        long long n;

#pragma omp parallel for default(none) private(n) shared(N, x, r) if (N>NumElementsUseThreading)
        for ( n=0; n<(long long)N; n++ )
        {
            reinterpret_cast<float(&)[2]>(r[n])[0] = reinterpret_cast< const float(&)[2]>(x[n])[0];
            reinterpret_cast<float(&)[2]>(r[n])[1] = -(reinterpret_cast< const float(&)[2]>(x[n])[1]);
        }
    }

    template <> EXPORTCPUCOREMATH 
    void conjugate(size_t N, const GT_Complex16* x, GT_Complex16* r)
    {
        long long n;

#pragma omp parallel for default(none) private(n) shared(N, x, r) if (N>NumElementsUseThreading)
        for ( n=0; n<(long long)N; n++ )
        {
            reinterpret_cast<double(&)[2]>(r[n])[0] = reinterpret_cast< const double(&)[2]>(x[n])[0];
            reinterpret_cast<double(&)[2]>(r[n])[1] = -(reinterpret_cast<const double(&)[2]>(x[n])[1]);
        }
    }

#endif // USE_MKL

    /// --------------------------------------------------------------------

    template <typename T> 
    void addEpsilon(size_t N, T* x)
    {
        typename realType<T>::Type eps = std::numeric_limits<typename realType<T>::Type>::epsilon();

        long long n;

#pragma omp parallel for default(none) private(n) shared(N, x, eps) if (N>NumElementsUseThreading)
        for (n=0; n<(long long)N; n++ )
        {
            if ( abs(x[n]) < eps )
            {
                x[n] += eps;
            }
        }
    }

    template EXPORTCPUCOREMATH void addEpsilon(size_t N, float* x);
    template EXPORTCPUCOREMATH void addEpsilon(size_t N, double* x);

    template <> EXPORTCPUCOREMATH 
    void addEpsilon(size_t N, GT_Complex8* x)
    {
        const float eps = std::numeric_limits<float>::epsilon();

        long long n;

#pragma omp parallel for private(n) shared(N, x, eps) if (N>NumElementsUseThreading)
        for (n=0; n<(long long)N; n++ )
        {
            if ( abs(x[n]) < eps )
            {
                reinterpret_cast<float(&)[2]>(x[n])[0] += eps;
            }
        }
    }

    template <> EXPORTCPUCOREMATH 
    void addEpsilon(size_t N, GT_Complex16* x)
    {
        const float eps = std::numeric_limits<double>::epsilon();

        long long n;

#pragma omp parallel for private(n) shared(N, x, eps) if (N>NumElementsUseThreading)
        for (n=0; n<(long long)N; n++ )
        {
            if ( abs(x[n]) < eps )
            {
                reinterpret_cast<double(&)[2]>(x[n])[0] += eps;
            }
        }
    }

    /// --------------------------------------------------------------------

    template <> EXPORTCPUCOREMATH void norm2(size_t N, const float* x, float& r)
    {
        lapack_int num = (lapack_int)N;
        lapack_int incx = 1;

#ifdef ILP_MODE_ON
        r = snrm2_(&num, (float*)(x), &incx);
#else
        if ( N < FourGBLimit )
        {
            r = snrm2_(&num, (float*)(x), &incx);
        }
        else
        {
            long long i;

            float sum(0);

#pragma omp parallel for reduction(+:sum) if (N>NumElementsUseThreading)
            for (i = 0; i < (long long)N; i++)
            {
                const float& re = x[i];
                sum += ( re*re );
            }

            r = std::sqrt(sum);
        }
#endif // ILP_MODE_ON
    }

    template <> EXPORTCPUCOREMATH void norm2(size_t N, const double* x, double& r)
    {
        lapack_int num = (lapack_int)N;
        lapack_int incx = 1;

#ifdef ILP_MODE_ON
        r = dnrm2_(&num, (double*)(x), &incx);
#else
        if ( N < FourGBLimit )
        {
            r = dnrm2_(&num, (double*)(x), &incx);
        }
        else
        {
            long long i;

            double sum(0);

#pragma omp parallel for reduction(+:sum) if (N>NumElementsUseThreading)
            for (i = 0; i < (long long)N; i++)
            {
                const double& re = x[i];
                sum += ( re*re );
            }

            r = std::sqrt(sum);
        }
#endif // ILP_MODE_ON
    }

    template <> EXPORTCPUCOREMATH void norm2(size_t N, const GT_Complex8* x, float& r)
    {
        lapack_int num = (lapack_int)N;
        lapack_int incx = 1;

#ifdef ILP_MODE_ON
        r = scnrm2_(&num, (lapack_complex_float*)(x), &incx);
#else
        if ( N < FourGBLimit )
        {
            r = scnrm2_(&num, (lapack_complex_float*)(x), &incx);
        }
        else
        {
            long long i;

            float sum(0);

#pragma omp parallel for reduction(+:sum) if (N>NumElementsUseThreading)
            for (i = 0; i < (long long)N; i++)
            {
                const std::complex<float>& c = x[i];
                const float re = c.real();
                const float im = c.imag();
                sum += ( (re*re) + (im * im) );
            }

            r = std::sqrt(sum);
        }
#endif // ILP_MODE_ON
    }

    template <> EXPORTCPUCOREMATH void norm2(size_t N, const GT_Complex16* x, double& r)
    {
        lapack_int num = (lapack_int)N;
        lapack_int incx = 1;

#ifdef ILP_MODE_ON
        r = dznrm2_(&num, (lapack_complex_double*)(x), &incx);
#else
        if ( N < FourGBLimit )
        {
            r = dznrm2_(&num, (lapack_complex_double*)(x), &incx);
        }
        else
        {
            long long i;

            double sum(0);

#pragma omp parallel for reduction(+:sum) if (N>NumElementsUseThreading)
            for (i = 0; i < (long long)N; i++)
            {
                const std::complex<double>& c = x[i];
                const double re = c.real();
                const double im = c.imag();
                sum += ( (re*re) + (im * im) );
            }

            r = std::sqrt(sum);
        }
#endif // ILP_MODE_ON
    }

    template <typename T> inline 
    typename realType<T>::Type norm2(size_t N, const T* x)
    {
        typename realType<T>::Type r;
        norm2(N, x, r);
        return r;
    }

    template EXPORTCPUCOREMATH float norm2(size_t N, const float* x);
    template EXPORTCPUCOREMATH double norm2(size_t N, const double* x);
    template EXPORTCPUCOREMATH float norm2(size_t N, const GT_Complex8* x);
    template EXPORTCPUCOREMATH double norm2(size_t N, const GT_Complex16* x);

    /// --------------------------------------------------------------------

    template <typename T> 
    void norm1(size_t N, const T* x, typename realType<T>::Type& r)
    {
        long long n;

        typename realType<T>::Type norm1Sum(0);

        #pragma omp parallel for reduction(+:norm1Sum) if (N>NumElementsUseThreading)
        for (n=0; n<(long long)N; n++)
        {
            const T& c = x[n];
            norm1Sum += GT_ABS(c);
        }

        r = norm1Sum;
    }

    template EXPORTCPUCOREMATH void norm1(size_t N, const float* x, float& r);
    template EXPORTCPUCOREMATH void norm1(size_t N, const double* x, double& r);

    template <> EXPORTCPUCOREMATH void norm1(size_t N, const GT_Complex8* x, float& r)
    {
        long long i;
        float sum = 0.0f;
        #pragma omp parallel for reduction(+:sum) if (N>NumElementsUseThreading)
        for (i = 0; i < (long long)N; i++)
        {
            const std::complex<float>& c = x[i];
            const float re = c.real();
            const float im = c.imag();
            sum += std::sqrt( (re*re) + (im * im) );
        }

        r = sum;
    }

    template <> EXPORTCPUCOREMATH void norm1(size_t N, const GT_Complex16* x, double& r)
    {
        long long i;
        double sum = 0.0;
        #pragma omp parallel for reduction(+:sum) if (N>NumElementsUseThreading)
        for (i = 0; i < (long long)N; i++)
        {
            const std::complex<double>& c = x[i];
            const double re = c.real();
            const double im = c.imag();
            sum += std::sqrt( (re*re) + (im * im) );
        }

        r = sum;
    }

    template <typename T> inline 
    typename realType<T>::Type norm1(size_t N, const T* x)
    {
        typename realType<T>::Type r;
        norm1(N, x, r);
        return r;
    }

    template EXPORTCPUCOREMATH float norm1(size_t N, const float* x);
    template EXPORTCPUCOREMATH double norm1(size_t N, const double* x);
    template EXPORTCPUCOREMATH float norm1(size_t N, const GT_Complex8* x);
    template EXPORTCPUCOREMATH double norm1(size_t N, const GT_Complex16* x);

    /// --------------------------------------------------------------------

    template <typename T> void dotc_64bit_mode(size_t N, const T* x, const T* y, T& r)
    {
        long long n;

        T sum(0);

        typename realType<T>::Type sa(0), sb(0);

#pragma omp parallel for reduction(+:sa) if (N>NumElementsUseThreading)
        for (n = 0; n < (long long)N; n++)
        {
            const typename realType<T>::Type a = x[n].real();
            const typename realType<T>::Type b = x[n].imag();
            const typename realType<T>::Type c = y[n].real();
            const typename realType<T>::Type d = y[n].imag();

            sa += (a*c + b*d);
            sb += (c*b - a*d);
        }

        r = T(sa, sb);
    }

    template <> EXPORTCPUCOREMATH void dotc(size_t N, const GT_Complex8* x, const GT_Complex8* y, GT_Complex8& r)
    {
        lapack_int num = (lapack_int)N;
        lapack_int incx=1, incy=1;

#ifdef ILP_MODE_ON
        cdotc_((lapack_complex_float*)(&r), &num, (lapack_complex_float*)(x), &incx, (lapack_complex_float*)(y), &incy);
#else
        if ( N < FourGBLimit )
        {
            cdotc_((lapack_complex_float*)(&r), &num, (lapack_complex_float*)(x), &incx, (lapack_complex_float*)(y), &incy);
        }
        else
        {
            dotc_64bit_mode(N, x, y, r);
        }
#endif // ILP_MODE_ON
    }

    template <> EXPORTCPUCOREMATH void dotc(size_t N, const GT_Complex16* x, const GT_Complex16* y, GT_Complex16& r)
    {
        lapack_int num = (lapack_int)N;
        lapack_int incx=1, incy=1;

#ifdef ILP_MODE_ON
        zdotc_((lapack_complex_double*)(&r), &num, (lapack_complex_double*)(x), &incx, (lapack_complex_double*)(y), &incy);
#else
        if ( N < FourGBLimit )
        {
            zdotc_((lapack_complex_double*)(&r), &num, (lapack_complex_double*)(x), &incx, (lapack_complex_double*)(y), &incy);
        }
        else
        {
            dotc_64bit_mode(N, x, y, r);
        }
#endif // ILP_MODE_ON
    }

    template <typename T> T dotc(size_t N, const T* x, const T* y)
    {
        T r;
        dotc(N, x, y, r);
        return r;
    }

    template EXPORTCPUCOREMATH GT_Complex8 dotc(size_t N, const GT_Complex8* x, const GT_Complex8* y);
    template EXPORTCPUCOREMATH GT_Complex16 dotc(size_t N, const GT_Complex16* x, const GT_Complex16* y);

    /// --------------------------------------------------------------------

    template <> EXPORTCPUCOREMATH void dotu(size_t N, const float* x, const float* y, float& r)
    {
        long long n;

        float res(0);

        #pragma omp parallel for reduction(+:res) if (N>NumElementsUseThreading)
        for (n=0; n<(long long)N; n++)
        {
            res += x[n]*y[n];
        }

        r = res;
    }

    template <> EXPORTCPUCOREMATH void dotu(size_t N, const double* x, const double* y, double& r)
    {
        long long n;

        double res(0);

        #pragma omp parallel for reduction(+:res) if (N>NumElementsUseThreading)
        for (n=0; n<(long long)N; n++)
        {
            res += x[n]*y[n];
        }

        r = res;
    }

    template <typename T> inline void dotu_64bit_mode(size_t N, const T* x, const T* y, T& r)
    {
        long long n;

        T sum(0);

        typename realType<T>::Type sa(0), sb(0);
#pragma omp parallel for reduction(+:sa) if (N>NumElementsUseThreading)
        for (n = 0; n < (long long)N; n++)
        {
            const typename realType<T>::Type a = x[n].real();
            const typename realType<T>::Type b = x[n].imag();
            const typename realType<T>::Type c = y[n].real();
            const typename realType<T>::Type d = y[n].imag();

            sa += (a*c - b*d);
            sb += (c*b + a*d);
        }

        r = T(sa, sb);
    }

    template <> EXPORTCPUCOREMATH void dotu(size_t N, const GT_Complex8* x, const GT_Complex8* y, GT_Complex8& r)
    {
        lapack_int num = (lapack_int)N;
        lapack_int incx=1, incy=1;

#ifdef ILP_MODE_ON
        cdotu_((lapack_complex_float*)(&r), &num, (lapack_complex_float*)(x), &incx, (lapack_complex_float*)(y), &incy);
#else
        if ( N < FourGBLimit )
        {
            cdotu_((lapack_complex_float*)(&r), &num, (lapack_complex_float*)(x), &incx, (lapack_complex_float*)(y), &incy);
        }
        else
        {
            dotu_64bit_mode(N, x, y, r);
        }
#endif // ILP_MODE_ON
    }

    template <> EXPORTCPUCOREMATH void dotu(size_t N, const GT_Complex16* x, const GT_Complex16* y, GT_Complex16& r)
    {
        lapack_int num = (lapack_int)N;
        lapack_int incx=1, incy=1;

#ifdef ILP_MODE_ON
        zdotu_((lapack_complex_double*)&r, &num, (lapack_complex_double*)(x), &incx, (lapack_complex_double*)(y), &incy);
#else
        if ( N < FourGBLimit )
        {
            zdotu_((lapack_complex_double*)&r, &num, (lapack_complex_double*)(x), &incx, (lapack_complex_double*)(y), &incy);
        }
        else
        {
            dotu_64bit_mode(N, x, y, r);
        }
#endif // ILP_MODE_ON
    }

    template <typename T> inline T dotu(size_t N, const T* x, const T* y)
    {
        T r;
        dotu(N, x, y, r);
        return r;
    }

    template EXPORTCPUCOREMATH float dotu(size_t N, const float* x, const float* y);
    template EXPORTCPUCOREMATH double dotu(size_t N, const double* x, const double* y);
    template EXPORTCPUCOREMATH GT_Complex8 dotu(size_t N, const GT_Complex8* x, const GT_Complex8* y);
    template EXPORTCPUCOREMATH GT_Complex16 dotu(size_t N, const GT_Complex16* x, const GT_Complex16* y);

    /// --------------------------------------------------------------------

    template <> EXPORTCPUCOREMATH void asum(size_t N, const float* x, float& r)
    {
        lapack_int num = (lapack_int)(N);
        lapack_int incx = 1;

#ifdef ILP_MODE_ON
        r = sasum_(&num, (float*)(x), &incx);
#else
        if ( N < FourGBLimit )
        {
            r = sasum_(&num, (float*)(x), &incx);
        }
        else
        {
            norm1(N, x, r);
        }
#endif // ILP_MODE_ON
    }

    template <> EXPORTCPUCOREMATH void asum(size_t N, const double* x, double& r)
    {
        lapack_int num = (lapack_int)(N);
        lapack_int incx = 1;

#ifdef ILP_MODE_ON
        r = dasum_(&num, (double*)(x), &incx);
#else
        if ( N < FourGBLimit )
        {
            r = dasum_(&num, (double*)(x), &incx);
        }
        else
        {
            norm1(N, x, r);
        }
#endif // ILP_MODE_ON
    }

    template <typename T> inline void asum_64bit_mode(size_t N, const T* x, typename realType<T>::Type& r)
    {
        long long i;
        typename realType<T>::Type sum(0);
        #pragma omp parallel for reduction(+:sum) if (N>NumElementsUseThreading)
        for (i = 0; i < (long long)N; i++)
        {
            const T& c = x[i];
            const typename realType<T>::Type re = c.real();
            const typename realType<T>::Type im = c.imag();
            sum += ( GT_ABS(re) + GT_ABS(im) );
        }

        r = sum;
    }

    template <> EXPORTCPUCOREMATH void asum(size_t N, const GT_Complex8* x, float& r)
    {
        lapack_int num = (lapack_int)(N);
        lapack_int incx = 1;

#ifdef ILP_MODE_ON
        r = scasum_(&num, (lapack_complex_float*)(x), &incx);
#else
        if ( N < FourGBLimit )
        {
            r = scasum_(&num, (lapack_complex_float*)(x), &incx);
        }
        else
        {
            asum_64bit_mode(N, x, r);
        }
#endif // ILP_MODE_ON
    }

    template <> EXPORTCPUCOREMATH void asum(size_t N, const GT_Complex16* x, double& r)
    {
        lapack_int num = (lapack_int)(N);
        lapack_int incx = 1;

#ifdef ILP_MODE_ON
        r = dzasum_(&num, (lapack_complex_double*)(x), &incx);
#else
        if ( N < FourGBLimit )
        {
            r = dzasum_(&num, (lapack_complex_double*)(x), &incx);
        }
        else
        {
            asum_64bit_mode(N, x, r);
        }
#endif // ILP_MODE_ON
    }

    template <typename T> inline typename realType<T>::Type asum(size_t N, const T* x)
    {
        typename realType<T>::Type r;
        asum(N, x, r);
        return r;
    }

    template EXPORTCPUCOREMATH float asum(size_t N, const float* x);
    template EXPORTCPUCOREMATH double asum(size_t N, const double* x);
    template EXPORTCPUCOREMATH float asum(size_t N, const GT_Complex8* x);
    template EXPORTCPUCOREMATH double asum(size_t N, const GT_Complex16* x);

    /// --------------------------------------------------------------------

    template <> EXPORTCPUCOREMATH size_t amax(size_t N, const float* x)
    {
        lapack_int num = (lapack_int)(N);
        lapack_int incx = 1;

        return isamax_(&num, (float*)(x), &incx);
    }

    template <> EXPORTCPUCOREMATH size_t amax(size_t N, const double* x)
    {
        lapack_int num = (lapack_int)(N);
        lapack_int incx = 1;

        return idamax_(&num, (double*)(x), &incx);
    }

    template <> EXPORTCPUCOREMATH size_t amax(size_t N, const GT_Complex8* x)
    {
        lapack_int num = (lapack_int)(N);
        lapack_int incx = 1;

        return icamax_(&num, (lapack_complex_float*)(x), &incx);
    }

    template <> EXPORTCPUCOREMATH size_t amax(size_t N, const GT_Complex16* x)
    {
        lapack_int num = (lapack_int)(N);
        lapack_int incx = 1;

        return izamax_(&num, (lapack_complex_double*)(x), &incx);
    }

    /// --------------------------------------------------------------------

#ifdef USE_MKL
    template <> EXPORTCPUCOREMATH void absolute(size_t N, const float* x, float* r)
    {
        vsAbs((lapack_int)N, x, r);
        GADGET_CHECK_THROW(vmlGetErrStatus()==0);
    }

    template <> EXPORTCPUCOREMATH void absolute(size_t N, const double* x, double* r)
    {
        vdAbs((lapack_int)N, x, r);
        GADGET_CHECK_THROW(vmlGetErrStatus()==0);
    }

    template <> EXPORTCPUCOREMATH void absolute(size_t N, const GT_Complex8* x, float* r)
    {
        vcAbs((lapack_int)N, (MKL_Complex8*)(x), r);
        GADGET_CHECK_THROW(vmlGetErrStatus()==0);
    }

    template <> EXPORTCPUCOREMATH void absolute(size_t N, const GT_Complex16* x, double* r)
    {
        vzAbs((lapack_int)N, (MKL_Complex16*)(x), r);
        GADGET_CHECK_THROW(vmlGetErrStatus()==0);
    }

#else

    template <typename T> 
    void absolute(size_t N, const T* x, typename realType<T>::Type* r)
    {
        long long n;

        #pragma omp parallel for default(none) private(n) shared(N, x, r) if (N>NumElementsUseThreading)
        for ( n=0; n<(long long)N; n++ )
        {
            r[n]= GT_ABS(x[n]);
        }
    }

    template EXPORTCPUCOREMATH void absolute(size_t N, const float* x, float* r);
    template EXPORTCPUCOREMATH void absolute(size_t N, const double* x, double* r);

    template <> EXPORTCPUCOREMATH void absolute(size_t N, const GT_Complex8* x, float* r)
    {
        long long n;

        #pragma omp parallel for default(none) private(n) shared(N, x, r) if (N>NumElementsUseThreading)
        for ( n=0; n<(long long)N; n++ )
        {
            const GT_Complex8& c = x[n];
            const float re = c.real();
            const float im = c.imag();
            r[n]= std::sqrt( (re*re) + (im * im) );
        }
    }

    template <> EXPORTCPUCOREMATH void absolute(size_t N, const GT_Complex16* x, double* r)
    {
        long long n;

        #pragma omp parallel for default(none) private(n) shared(N, x, r) if (N>NumElementsUseThreading)
        for ( n=0; n<(long long)N; n++ )
        {
            const GT_Complex16& c = x[n];
            const double re = c.real();
            const double im = c.imag();
            r[n]= std::sqrt( (re*re) + (im * im) );
        }
    }
#endif // USE_MKL

    template <> EXPORTCPUCOREMATH 
    void absolute(size_t N, const std::complex<float>* x, std::complex<float>* r)
    {
        try
        {
            long long n;

            #pragma omp parallel for default(none) private(n) shared(N, x, r) if (N>NumElementsUseThreading)
            for ( n=0; n<(long long)N; n++ )
            {
                const std::complex<float>& c = x[n];
                const float re = c.real();
                const float im = c.imag();

                reinterpret_cast<float(&)[2]>(r[n])[0] = std::sqrt( (re*re) + (im * im) );
                reinterpret_cast<float(&)[2]>(r[n])[1] = 0;
            }
        }
        catch(...)
        {
            GADGET_THROW("Error happened in absolute(size_t N, const std::complex<float>* x, std::complex<float>* r) ... ");
        }
    }

    template <> EXPORTCPUCOREMATH 
    void absolute(size_t N, const std::complex<double>* x, std::complex<double>* r)
    {
        try
        {
            long long n;

            #pragma omp parallel for default(none) private(n) shared(N, x, r) if (N>NumElementsUseThreading)
            for ( n=0; n<(long long)N; n++ )
            {
                const std::complex<double>& c = x[n];
                const double re = c.real();
                const double im = c.imag();

                reinterpret_cast<double(&)[2]>(r[n])[0] = std::sqrt( (re*re) + (im * im) );
                reinterpret_cast<double(&)[2]>(r[n])[1] = 0;
            }
        }
        catch(...)
        {
            GADGET_THROW("Error happened in absolute(size_t N, const std::complex<double>* x, std::complex<double>* r) ... ");
        }
    }

    /// --------------------------------------------------------------------

#ifdef USE_MKL

    template <> EXPORTCPUCOREMATH void argument(size_t N, const GT_Complex8* x, float* r)
    {
        vcArg((lapack_int)N, (MKL_Complex8*)(x), r);
        GADGET_CHECK_THROW(vmlGetErrStatus()==0);
    }

    template <> EXPORTCPUCOREMATH void argument(size_t N, const GT_Complex16* x, double* r)
    {
        vzArg((lapack_int)N, (MKL_Complex16*)(x), r);
        GADGET_CHECK_THROW(vmlGetErrStatus()==0);
    }

#else

    template <typename T> 
    void argument(size_t N, const T* x, typename realType<T>::Type* r)
    {
        long long n;

        #pragma omp parallel for default(none) private(n) shared(N, x, r) if (N>NumElementsUseThreading)
        for ( n=0; n<(long long)N; n++ )
        {
            r[n] = std::arg( x[n] );
        }
    }

    template EXPORTCPUCOREMATH void argument(size_t N, const GT_Complex8* x, float* r);
    template EXPORTCPUCOREMATH void argument(size_t N, const GT_Complex16* x, double* r);
#endif // USE_MKL

    /// --------------------------------------------------------------------

    template <typename T> 
    void inv(size_t N, const T* x, T* r)
    {
        T v(1.0);
        long long n;

        #pragma omp parallel for default(none) private(n) shared(N, x, r, v) if (N>NumElementsUseThreading)
        for ( n=0; n<(long long)N; n++ )
        {
            r[n] = v/x[n];
        }
    }

#ifdef USE_MKL

    template <> EXPORTCPUCOREMATH void inv(size_t N, const float* x, float* r)
    {
        vsInv((lapack_int)N, x, r);
        GADGET_CHECK_THROW(vmlGetErrStatus()==0);
    }

    template <> EXPORTCPUCOREMATH void inv(size_t N, const double* x, double* r)
    {
        vdInv((lapack_int)N, x, r);
        GADGET_CHECK_THROW(vmlGetErrStatus()==0);
    }

#else
    template EXPORTCPUCOREMATH void inv(size_t N, const float* x, float* r);
    template EXPORTCPUCOREMATH void inv(size_t N, const double* x, double* r);
#endif // USE_MKL

    template EXPORTCPUCOREMATH void inv(size_t N, const GT_Complex8* x, GT_Complex8* r);
    template EXPORTCPUCOREMATH void inv(size_t N, const GT_Complex16* x, GT_Complex16* r);

    /// --------------------------------------------------------------------

    template<typename T> 
    void conv2(size_t RO, size_t E1, size_t num, const T* x, size_t kRO, size_t kE1, const T* y, T* z)
    {
        try
        {
            long long halfKRO = (long long)(kRO/2);
            long long halfKE1 = (long long)(kE1/2);

            hoNDArray<T> flipY(2*halfKRO+1, 2*halfKE1+1);
            T* pKer = flipY.begin();

            long long n;
            long long ro, e1;

            // flip the kernel
            for ( e1=0; e1<(long long)kE1; e1++ )
            {
                long long flip_e1 = 2*halfKE1 - e1;

                for ( ro=0; ro<(long long)kRO; ro++ )
                {
                    long long flip_ro = 2*halfKRO - ro;

                    flipY(flip_ro, flip_e1) = y[ro+e1*kRO];
                }
            }

            // perform the convolution
            #pragma omp parallel for default(none) private(n, ro, e1) shared(num, x, RO, E1, z, halfKRO, halfKE1, pKer)
            for ( n=0; n<(long long)num; n++ )
            {
                const T* pX = x + n*RO*E1;
                T* pZ = z + n*RO*E1;

                long long kro, ke1, dro, de1;

                for ( e1=0; e1<(long long)E1; e1++ )
                {
                    for ( ro=0; ro<(long long)RO; ro++ )
                    {
                        pZ[ro + e1*RO] = 0;
                        for ( ke1=-halfKE1; ke1<=halfKE1; ke1++ )
                        {
                            de1 = ke1 + e1;
                            if ( de1 < 0 )
                            {
                                de1 += E1;
                            }
                            else if ( de1 >= (long long)E1 )
                            {
                                de1 -= E1;
                            }

                            for ( kro=-halfKRO; kro<=halfKRO; kro++ )
                            {
                                dro = kro + ro;
                                if ( dro < 0 )
                                {
                                    dro += RO;
                                }
                                else if ( dro >= (long long)RO )
                                {
                                    dro -= RO;
                                }

                                pZ[ro + e1*RO] += pKer[ kro+halfKRO + (ke1+halfKE1) * (2*halfKRO+1) ] * pX[dro + de1*RO];
                            }
                        }
                    }
                }
            }
        }
        catch(...)
        {
            GADGET_THROW("Errors happened in conv2(size_t RO, size_t E1, size_t num, const T* x, size_t kRO, size_t kE1, const T* y, T* z) ... ");
        }
    }

    template EXPORTCPUCOREMATH void conv2(size_t RO, size_t E1, size_t num, const float* x, size_t kRO, size_t kE1, const float* y, float* z);
    template EXPORTCPUCOREMATH void conv2(size_t RO, size_t E1, size_t num, const double* x, size_t kRO, size_t kE1, const double* y, double* z);
    template EXPORTCPUCOREMATH void conv2(size_t RO, size_t E1, size_t num, const GT_Complex8* x, size_t kRO, size_t kE1, const GT_Complex8* y, GT_Complex8* z);
    template EXPORTCPUCOREMATH void conv2(size_t RO, size_t E1, size_t num, const GT_Complex16* x, size_t kRO, size_t kE1, const GT_Complex16* y, GT_Complex16* z);

    /// --------------------------------------------------------------------

    template<typename T> 
    void conv3(size_t RO, size_t E1, size_t E2, size_t num, const T* x, size_t kRO, size_t kE1, size_t kE2, const T* y, T* z)
    {
        try
        {
            long long halfKRO = (long long)(kRO/2);
            long long halfKE1 = (long long)(kE1/2);
            long long halfKE2 = (long long)(kE2/2);

            hoNDArray<T> flipY(2*halfKRO+1, 2*halfKE1+1, 2*halfKE2+1);
            T* pKer = flipY.begin();

            long long n, e2;
            long long ro, e1;

            // flip the kernel
            for ( e2=0; e2<(long long)kE2; e2++ )
            {
                long long flip_e2 = 2*halfKE2 - e2;

                for ( e1=0; e1<(long long)kE1; e1++ )
                {
                    long long flip_e1 = 2*halfKE1 - e1;

                    for ( ro=0; ro<(long long)kRO; ro++ )
                    {
                        long long flip_ro = 2*halfKRO - ro;

                        flipY(flip_ro, flip_e1, flip_e2) = y[ro+e1*kRO+e2*kRO*kE1];
                    }
                }
            }

            // perform the convolution
            #pragma omp parallel for default(none) private(n) shared(num, x, RO, E1, E2, z, halfKRO, halfKE1, halfKE2, pKer) if ( num > 8 )
            for ( n=0; n<(long long)num; n++ )
            {
                const T* pX = x + n*RO*E1*E2;
                T* pZ = z + n*RO*E1*E2;

                long long kro, ke1, ke2, dro, de1, de2;

                #pragma omp parallel for default(none) private(ro, e1, e2, kro, ke1, ke2, dro, de1, de2) shared(pX, RO, E1, E2, pZ, halfKRO, halfKE1, halfKE2, pKer)
                for ( e2=0; e2<(long long)E2; e2++ )
                {
                    for ( e1=0; e1<(long long)E1; e1++ )
                    {
                        for ( ro=0; ro<(long long)RO; ro++ )
                        {
                            pZ[ro + e1*RO + e2*RO*E1] = 0;
                            for ( ke2=-halfKE2; ke2<=halfKE2; ke2++ )
                            {
                                de2 = ke2 + e2;
                                if ( de2 < 0 )
                                {
                                    de2 += E2;
                                }
                                else if ( de2 >= (long long)E2 )
                                {
                                    de2 -= E2;
                                }

                                for ( ke1=-halfKE1; ke1<=halfKE1; ke1++ )
                                {
                                    de1 = ke1 + e1;
                                    if ( de1 < 0 )
                                    {
                                        de1 += E1;
                                    }
                                    else if ( de1 >= (long long)E1 )
                                    {
                                        de1 -= E1;
                                    }

                                    for ( kro=-halfKRO; kro<=halfKRO; kro++ )
                                    {
                                        dro = kro + ro;
                                        if ( dro < 0 )
                                        {
                                            dro += RO;
                                        }
                                        else if ( dro >= (long long)RO )
                                        {
                                            dro -= RO;
                                        }

                                        pZ[ro + e1*RO + e2*RO*E1] += pKer[ kro+halfKRO + (ke1+halfKE1)*(2*halfKRO+1) + (ke2+halfKE2)*(2*halfKRO+1)*(2*halfKE1+1) ] * pX[dro + de1*RO + de2*RO*E1];
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        catch(...)
        {
            GADGET_THROW("Errors happened in conv3(size_t RO, size_t E1, size_t E2, size_t num, const T* x, size_t kRO, size_t kE1, size_t kE2, const T* y, T* z) ... ");
        }
    }

    template EXPORTCPUCOREMATH void conv3(size_t RO, size_t E1, size_t E2, size_t num, const float* x, size_t kRO, size_t kE1, size_t kE2, const float* y, float* z);
    template EXPORTCPUCOREMATH void conv3(size_t RO, size_t E1, size_t E2, size_t num, const double* x, size_t kRO, size_t kE1, size_t kE2, const double* y, double* z);
    template EXPORTCPUCOREMATH void conv3(size_t RO, size_t E1, size_t E2, size_t num, const GT_Complex8* x, size_t kRO, size_t kE1, size_t kE2, const GT_Complex8* y, GT_Complex8* z);
    template EXPORTCPUCOREMATH void conv3(size_t RO, size_t E1, size_t E2, size_t num, const GT_Complex16* x, size_t kRO, size_t kE1, size_t kE2, const GT_Complex16* y, GT_Complex16* z);

    /// --------------------------------------------------------------------

    template <typename T> 
    struct hoCompAscending
    {
        bool operator() (T a, T b) { return (a>=b); }
    };

    template <typename T> 
    struct hoCompDescending
    {
        bool operator() (T a, T b) { return (a<b); }
    };

    template <typename T> 
    void sort(size_t N, const T* x, T* r, bool isascending)
    {
        if ( r != x )
        {
            memcpy(r, x, sizeof(T)*N);
        }

        if ( isascending )
        {
            hoCompAscending<T> obj;
            std::sort(r, r+N, obj);
        }
        else
        {
            hoCompDescending<T> obj;
            std::sort(r, r+N, obj);
        }
    }

    template EXPORTCPUCOREMATH void sort(size_t N, const float* x, float* r, bool isascending);
    template EXPORTCPUCOREMATH void sort(size_t N, const double* x, double* r, bool isascending);

    /// --------------------------------------------------------------------

    template<typename T> 
    void fill( size_t N, T* x, T v)
    {
        long long n;
        #pragma omp parallel for default(none) private(n) shared(N, x, v)
        for ( n=0; n<(long long)N; n++ )
        {
            x[n] = v;
        }
    }

    template EXPORTCPUCOREMATH void fill(size_t N, float* x, float v);
    template EXPORTCPUCOREMATH void fill(size_t N, double* x, double v);
    template EXPORTCPUCOREMATH void fill(size_t N, GT_Complex8* x, GT_Complex8 v);
    template EXPORTCPUCOREMATH void fill(size_t N, GT_Complex16* x, GT_Complex16 v);
}}
