/** \file  hoNDMath_util.cpp
\brief math function utility
*/

#include "hoNDMath_util.h"
#include <stdexcept>
#include <iostream>
#include "GadgetronCommon.h"

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

    #define isamin_ isamin
    #define idamin_ idamin
    #define icamin_ icamin
    #define izamin_ izamin

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

        /// Finds the index of the element with the smallest absolute value.
        lapack_int isamin_(lapack_int* N, float* x, lapack_int* incx);
        lapack_int idamin_(lapack_int* N, double* x, lapack_int* incx);
        lapack_int icamin_(lapack_int* N, lapack_complex_float* x, lapack_int* incx);
        lapack_int izamin_(lapack_int* N, lapack_complex_double* x, lapack_int* incx);

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

            #ifdef WIN32
                r[n] = T( re2 + ar*re1 - ai*im1, im2 + ar*im1 + ai*re1 );
            #else
                r[n].real() = (re2 + ar*re1 - ai*im1);
                r[n].imag() = (im2 + ar*im1 + ai*re1);
            #endif // WIN32
        }
    }

    template <> EXPORTCPUCOREMATH void axpy(float a, size_t N, const float* x, const float* y, float* r)
    {
        try
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
        catch(...)
        {
            GADGET_THROW("Error happened in axpy(float a, size_t N, const float* x, const float* y, float* r) ... ");
        }
    }

    template <> EXPORTCPUCOREMATH void axpy(double a, size_t N, const double* x, const double* y, double* r)
    {
        try
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
        catch(...)
        {
            GADGET_THROW("Error happened in axpy(double a, size_t N, const double* x, const double* y, double* r) ... ");
        }
    }

    template <> EXPORTCPUCOREMATH void axpy(GT_Complex8 a, size_t N, const GT_Complex8* x, const GT_Complex8* y, GT_Complex8* r)
    {
        try
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
        catch(...)
        {
            GADGET_THROW("Error happened in axpy(GT_Complex8 a, size_t N, const GT_Complex8* x, const GT_Complex8* y, GT_Complex8* r) ... ");
        }
    }

    template <> EXPORTCPUCOREMATH void axpy(GT_Complex16 a, size_t N, const GT_Complex16* x, const GT_Complex16* y, GT_Complex16* r)
    {
        try
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
        catch(...)
        {
            GADGET_THROW("Error happened in axpy(GT_Complex16 a, size_t N, const GT_Complex16* x, const GT_Complex16* y, GT_Complex16* r) ... ");
        }
    }

    /// --------------------------------------------------------------------

#ifdef USE_MKL
    template <> EXPORTCPUCOREMATH void add(size_t N, const float* x, const float* y, float* r)
    {
        vsAdd(N, x, y, r);
        GADGET_CHECK_THROW(vmlGetErrStatus()==0);
    }

    template <> EXPORTCPUCOREMATH void add(size_t N, const double* x, const double* y, double* r)
    {
        vdAdd(N, x, y, r);
        GADGET_CHECK_THROW(vmlGetErrStatus()==0);
    }

    template <> EXPORTCPUCOREMATH void add(size_t N, const GT_Complex8* x, const GT_Complex8* y, GT_Complex8* r)
    {
        vcAdd(N, (MKL_Complex8*)(x), (MKL_Complex8*)(y), (MKL_Complex8*)(r));
        GADGET_CHECK_THROW(vmlGetErrStatus()==0);
    }

    template <> EXPORTCPUCOREMATH void add(size_t N, const GT_Complex16* x, const GT_Complex16* y, GT_Complex16* r)
    {
        vzAdd(N, (MKL_Complex16*)(x), (MKL_Complex16*)(y), (MKL_Complex16*)(r));
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
        vsSub(N, x, y, r);
        GADGET_CHECK_THROW(vmlGetErrStatus()==0);
    }

    template <> EXPORTCPUCOREMATH void subtract(size_t N, const double* x, const double* y, double* r)
    {
        vdSub(N, x, y, r);
        GADGET_CHECK_THROW(vmlGetErrStatus()==0);
    }

    template <> EXPORTCPUCOREMATH void subtract(size_t N, const GT_Complex8* x, const GT_Complex8* y, GT_Complex8* r)
    {
        vcSub(N, (MKL_Complex8*)(x), (MKL_Complex8*)(y), (MKL_Complex8*)(r));
        GADGET_CHECK_THROW(vmlGetErrStatus()==0);
    }

    template <> EXPORTCPUCOREMATH void subtract(size_t N, const GT_Complex16* x, const GT_Complex16* y, GT_Complex16* r)
    {
        vzSub(N, (MKL_Complex16*)(x), (MKL_Complex16*)(y), (MKL_Complex16*)(r));
        GADGET_CHECK_THROW(vmlGetErrStatus()==0);
    }

#else

    template <typename T> void subtract(size_t N, const T* x, const T* y, T* r)
    {
        try
        {
            long long n;
    #pragma omp parallel for default(none) private(n) shared(N, x, y, r) if (N>NumElementsUseThreading)
            for ( n=0; n<(long long)N; n++ )
            {
                r[n] = x[n] - y[n];
            }
        }
        catch(...)
        {
            GADGET_THROW("Error happened in subtract(size_t N, const T* x, const T* y, T* r) ... ");
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
        vsMul(N, x, y, r);
        GADGET_CHECK_THROW(vmlGetErrStatus()==0);
    }

    template <> EXPORTCPUCOREMATH void multiply(size_t N, const double* x, const double* y, double* r)
    {
        vdMul(N, x, y, r);
        GADGET_CHECK_THROW(vmlGetErrStatus()==0);
    }

    template <> EXPORTCPUCOREMATH void multiply(size_t N, const GT_Complex8* x, const GT_Complex8* y, GT_Complex8* r)
    {
        vcMul(N, (MKL_Complex8*)(x), (MKL_Complex8*)(y), (MKL_Complex8*)(r));
        GADGET_CHECK_THROW(vmlGetErrStatus()==0);
    }

    template <> EXPORTCPUCOREMATH void multiply(size_t N, const GT_Complex16* x, const GT_Complex16* y, GT_Complex16* r)
    {
        vzMul(N, (MKL_Complex16*)(x), (MKL_Complex16*)(y), (MKL_Complex16*)(r));
        GADGET_CHECK_THROW(vmlGetErrStatus()==0);
    }

#else

    template <typename T> 
    void multiply(size_t N, const T* x, const T* y, T* r)
    {
        try
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
        catch(...)
        {
            GADGET_THROW("Error happened in multiply(size_t N, const T* x, const T* y, T* r) ... ");
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

            #ifdef WIN32
                r[i] = std::complex<float>(a*c-b*d, a*d+b*c);
            #else
                r[i].real() = (a*c-b*d);
                r[i].imag() = (a*d+b*c);
            #endif // WIN32
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

            #ifdef WIN32
                r[i] = std::complex<double>(a*c-b*d, a*d+b*c);
            #else
                r[i].real() = (a*c-b*d);
                r[i].imag() = (a*d+b*c);
            #endif // WIN32
        }
    }

#endif // USE_MKL

    /// --------------------------------------------------------------------

#ifdef USE_MKL
    template <> EXPORTCPUCOREMATH void divide(size_t N, const float* x, const float* y, float* r)
    {
        vsDiv(N, x, y, r);
        GADGET_CHECK_THROW(vmlGetErrStatus()==0);
    }

    template <> EXPORTCPUCOREMATH void divide(size_t N, const double* x, const double* y, double* r)
    {
        vdDiv(N, x, y, r);
        GADGET_CHECK_THROW(vmlGetErrStatus()==0);
    }

    template <> EXPORTCPUCOREMATH void divide(size_t N, const GT_Complex8* x, const GT_Complex8* y, GT_Complex8* r)
    {
        vcDiv(N, (MKL_Complex8*)(x), (MKL_Complex8*)(y), (MKL_Complex8*)(r));
        GADGET_CHECK_THROW(vmlGetErrStatus()==0);
    }

    template <> EXPORTCPUCOREMATH void divide(size_t N, const GT_Complex16* x, const GT_Complex16* y, GT_Complex16* r)
    {
        vzDiv(N, (MKL_Complex16*)(x), (MKL_Complex16*)(y), (MKL_Complex16*)(r));
        GADGET_CHECK_THROW(vmlGetErrStatus()==0);
    }

#else

    template <typename T> 
    void divide(size_t N, const T* x, const T* y, T* r)
    {
        try
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
        catch(...)
        {
            GADGET_THROW("Error happened in divide(size_t N, const T* x, const T* y, T* r) ... ");
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

            #ifdef WIN32
                r[i] = std::complex<float>( (a*c+b*d)*m, (b*c-a*d)*m );
            #else
                r[i].real() = ((a*c+b*d)*m);
                r[i].imag() = ((b*c-a*d)*m);
            #endif // WIN32
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

            #ifdef WIN32
                r[i] = std::complex<double>( (a*c+b*d)*m, (b*c-a*d)*m );
            #else
                r[i].real() = ((a*c+b*d)*m);
                r[i].imag() = ((b*c-a*d)*m);
            #endif // WIN32
        }
    }
#endif // USE_MKL

    /// --------------------------------------------------------------------

#ifdef USE_MKL
    template <> EXPORTCPUCOREMATH void sqrt(size_t N, const float* x, float* r)
    {
        vsSqrt(N, x, r);
        GADGET_CHECK_THROW(vmlGetErrStatus()==0);
    }

    template <> EXPORTCPUCOREMATH void sqrt(size_t N, const double* x, double* r)
    {
        vdSqrt(N, x, r);
        GADGET_CHECK_THROW(vmlGetErrStatus()==0);
    }

    template <> EXPORTCPUCOREMATH void sqrt(size_t N, const GT_Complex8* x, GT_Complex8* r)
    {
        vcSqrt(N, (MKL_Complex8*)(x), (MKL_Complex8*)(r));
        GADGET_CHECK_THROW(vmlGetErrStatus()==0);
    }

    template <> EXPORTCPUCOREMATH void sqrt(size_t N, const GT_Complex16* x, GT_Complex16* r)
    {
        vzSqrt(N, (MKL_Complex16*)(x), (MKL_Complex16*)(r));
        GADGET_CHECK_THROW(vmlGetErrStatus()==0);
    }

#else

    template <typename T> 
    void sqrt(size_t N, const T* x, T* r)
    {
        try
        {
            long long n;
#pragma omp parallel for default(none) private(n) shared(N, x, r) if (N>NumElementsUseThreading)
            for ( n=0; n<(long long)N; n++ )
            {
                r[n] = std::sqrt(x[n]);
            }
        }
        catch(...)
        {
            GADGET_THROW("Error happened in sqrt(size_t N, const T* x, T* r) ... ");
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
                if ( v2 < v )
                {
                    v = v2;
                    ind = n;
                }
            }

            r = x[ind];
        }
        catch(...)
        {
            GADGET_THROW("Error happened in minAbsolute(size_t N, const T* x, T& r, size_t& ind) ... ");
        }
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
        vcMulByConj(N, (MKL_Complex8*)(x), (MKL_Complex8*)(y), (MKL_Complex8*)(r));
        GADGET_CHECK_THROW(vmlGetErrStatus()==0);
    }

    template <> EXPORTCPUCOREMATH void multiplyConj(size_t N, const GT_Complex16* x, const GT_Complex16* y, GT_Complex16* r)
    {
        vzMulByConj(N, (MKL_Complex16*)(x), (MKL_Complex16*)(y), (MKL_Complex16*)(r));
        GADGET_CHECK_THROW(vmlGetErrStatus()==0);
    }

#else

    template <typename T> 
    void multiplyConj(size_t N, const T* x, const T* y, T* r)
    {
        try
        {
            long long n;

#pragma omp parallel for private(n) shared(N, x, y, r) if (N>NumElementsUseThreading)
            for ( n=0; n<(long long)N; n++ )
            {
                const typename realType<T>::Type a = x[n].real();
                const typename realType<T>::Type b = x[n].imag();
                const typename realType<T>::Type c = y[n].real();
                const typename realType<T>::Type d = y[n].imag();

                #ifdef WIN32
                    r[n] = T(a*c + b*d, c*b - a*d);
                #else
                    r[n].real() = (a*c + b*d);
                    r[n].imag() = (c*b - a*d);
                #endif // WIN32
            }
        }
        catch(...)
        {
            GADGET_THROW("Error happened in multiplyConj(const hoNDArray<T>& x, const hoNDArray<T>& y, hoNDArray<T>& r) ... ");
        }
    }

    template EXPORTCPUCOREMATH void multiplyConj(size_t N, const GT_Complex8* x, const GT_Complex8* y, GT_Complex8* r);
    template EXPORTCPUCOREMATH void multiplyConj(size_t N, const GT_Complex16* x, const GT_Complex16* y, GT_Complex16* r);
#endif // USE_MKL
    /// --------------------------------------------------------------------

#ifdef USE_MKL

    template <> EXPORTCPUCOREMATH void conjugate(size_t N, const GT_Complex8* x, GT_Complex8* r)
    {
        vcConj(N, (MKL_Complex8*)(x), (MKL_Complex8*)(r));
        GADGET_CHECK_THROW(vmlGetErrStatus()==0);
    }

    template <> EXPORTCPUCOREMATH void conjugate(size_t N, const GT_Complex16* x, GT_Complex16* r)
    {
        vzConj(N, (MKL_Complex16*)(x), (MKL_Complex16*)(r));
        GADGET_CHECK_THROW(vmlGetErrStatus()==0);
    }

#else

    template <typename T> 
    void conjugate(size_t N, const T* x, T* r)
    {
        try
        {
            long long n;

#pragma omp parallel for default(none) private(n) shared(N, x, r) if (N>NumElementsUseThreading)
            for ( n=0; n<(long long)N; n++ )
            {
                #ifdef WIN32
                    r[n] = std::conj(x[n]);
                #else
                    r[n].real() = x[n].real();
                    r[n].imag() = -x[n].imag();
                #endif // WIN32
            }
        }
        catch(...)
        {
            GADGET_THROW("Error happened in conjugate(size_t N, const T* x, T* r) ... ");
        }
    }

    template EXPORTCPUCOREMATH void conjugate(size_t N, const GT_Complex8* x, GT_Complex8* r);
    template EXPORTCPUCOREMATH void conjugate(size_t N, const GT_Complex16* x, GT_Complex16* r);
#endif // USE_MKL

    /// --------------------------------------------------------------------

    template <typename T> 
    void addEpsilon(size_t N, T* x)
    {
        try
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
        catch(...)
        {
            GADGET_THROW("Error happened in addEpsilon(size_t N, T* x) ... ");
        }
    }

    template EXPORTCPUCOREMATH void addEpsilon(size_t N, float* x);
    template EXPORTCPUCOREMATH void addEpsilon(size_t N, double* x);
    template EXPORTCPUCOREMATH void addEpsilon(size_t N, GT_Complex8* x);
    template EXPORTCPUCOREMATH void addEpsilon(size_t N, GT_Complex16* x);

    /// --------------------------------------------------------------------

    template <> EXPORTCPUCOREMATH void norm2(size_t N, const float* x, float& r)
    {
        try
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
        catch(...)
        {
            GADGET_THROW("Error happened in norm2(size_t N, const float* x, float& r) ... ");
        }
    }

    template <> EXPORTCPUCOREMATH void norm2(size_t N, const double* x, double& r)
    {
        try
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
        catch(...)
        {
            GADGET_THROW("Error happened in norm2(size_t N, const double* x, double& r) ... ");
        }
    }

    template <> EXPORTCPUCOREMATH void norm2(size_t N, const GT_Complex8* x, float& r)
    {
        try
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
        catch(...)
        {
            GADGET_THROW("Error happened in norm2(size_t N, const GT_Complex8* x, float& r) ... ");
        }
    }

    template <> EXPORTCPUCOREMATH void norm2(size_t N, const GT_Complex16* x, double& r)
    {
        try
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
        catch(...)
        {
            GADGET_THROW("Error happened in norm2(size_t N, const GT_Complex16* x, double& r) ... ");
        }
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
        try
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
        catch(...)
        {
            GADGET_THROW("Error happened in norm1(size_t N, const T* x, typename realType<T>::Type& r) ... ");
        }
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
        try
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
        catch(...)
        {
            GADGET_THROW("Error happened in dotc(size_t N, const GT_Complex8* x, const GT_Complex8* y, GT_Complex8& r) ... ");
        }
    }

    template <> EXPORTCPUCOREMATH void dotc(size_t N, const GT_Complex16* x, const GT_Complex16* y, GT_Complex16& r)
    {
        try
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
        catch(...)
        {
            GADGET_THROW("Error happened in dotc(size_t N, const GT_Complex16* x, const GT_Complex16* y, GT_Complex16& r) ... ");
        }
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
        try
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
        catch(...)
        {
            GADGET_THROW("Error happened in dotu(size_t N, const float* x, const float* y, float& r) ... ");
        }
    }

    template <> EXPORTCPUCOREMATH void dotu(size_t N, const double* x, const double* y, double& r)
    {
        try
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
        catch(...)
        {
            GADGET_THROW("Error happened in dotu(size_t N, const double* x, const double* y, double& r) ... ");
        }
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
        try
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
        catch(...)
        {
            GADGET_THROW("Error happened in dotu(size_t N, const GT_Complex8* x, const GT_Complex8* y, GT_Complex8& r) ... ");
        }
    }

    template <> EXPORTCPUCOREMATH void dotu(size_t N, const GT_Complex16* x, const GT_Complex16* y, GT_Complex16& r)
    {
        try
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
        catch(...)
        {
            GADGET_THROW("Error happened in dotu(size_t N, const GT_Complex16* x, const GT_Complex16* y, GT_Complex16& r) ... ");
        }
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
        try
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
        catch(...)
        {
            GADGET_THROW("Error happened in asum(size_t N, const float* x, float& r) ... ");
        }
    }

    template <> EXPORTCPUCOREMATH void asum(size_t N, const double* x, double& r)
    {
        try
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
        catch(...)
        {
            GADGET_THROW("Error happened in asum(size_t N, const double* x, double& r) ... ");
        }
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
        try
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
        catch(...)
        {
            GADGET_THROW("Error happened in asum(size_t N, const GT_Complex8* x, float& r) ... ");
        }
    }

    template <> EXPORTCPUCOREMATH void asum(size_t N, const GT_Complex16* x, double& r)
    {
        try
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
        catch(...)
        {
            GADGET_THROW("Error happened in asum(size_t N, const GT_Complex16* x, double& r) ... ");
        }
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

    template <> EXPORTCPUCOREMATH size_t amin(size_t N, const float* x)
    {
        try
        {
            lapack_int num = (lapack_int)(N);
            lapack_int incx = 1;

            return isamin_(&num, (float*)(x), &incx);
        }
        catch(...)
        {
            GADGET_THROW("Error happened in size_t amin(size_t N, const float* x) ... ");
        }
    }

    template <> EXPORTCPUCOREMATH size_t amin(size_t N, const double* x)
    {
        try
        {
            lapack_int num = (lapack_int)(N);
            lapack_int incx = 1;

            return idamin_(&num, (double*)(x), &incx);
        }
        catch(...)
        {
            GADGET_THROW("Error happened in size_t amin(size_t N, const double* x) ... ");
        }
    }

    template <> EXPORTCPUCOREMATH size_t amin(size_t N, const GT_Complex8* x)
    {
        try
        {
            lapack_int num = (lapack_int)(N);
            lapack_int incx = 1;

            return icamin_(&num, (lapack_complex_float*)(x), &incx);
        }
        catch(...)
        {
            GADGET_THROW("Error happened in size_t amin(size_t N, const GT_Complex8* x) ... ");
        }
    }

    template <> EXPORTCPUCOREMATH size_t amin(size_t N, const GT_Complex16* x)
    {
        try
        {
            lapack_int num = (lapack_int)(N);
            lapack_int incx = 1;

            return izamin_(&num, (lapack_complex_double*)(x), &incx);
        }
        catch(...)
        {
            GADGET_THROW("Error happened in size_t amin(size_t N, const GT_Complex16* x) ... ");
        }
    }

    /// --------------------------------------------------------------------

    template <> EXPORTCPUCOREMATH size_t amax(size_t N, const float* x)
    {
        try
        {
            lapack_int num = (lapack_int)(N);
            lapack_int incx = 1;

            return isamax_(&num, (float*)(x), &incx);
        }
        catch(...)
        {
            GADGET_THROW("Error happened in size_t amax(size_t N, const float* x) ... ");
        }
    }

    template <> EXPORTCPUCOREMATH size_t amax(size_t N, const double* x)
    {
        try
        {
            lapack_int num = (lapack_int)(N);
            lapack_int incx = 1;

            return idamax_(&num, (double*)(x), &incx);
        }
        catch(...)
        {
            GADGET_THROW("Error happened in size_t amax(size_t N, const double* x) ... ");
        }
    }

    template <> EXPORTCPUCOREMATH size_t amax(size_t N, const GT_Complex8* x)
    {
        try
        {
            lapack_int num = (lapack_int)(N);
            lapack_int incx = 1;

            return icamax_(&num, (lapack_complex_float*)(x), &incx);
        }
        catch(...)
        {
            GADGET_THROW("Error happened in size_t amax(size_t N, const GT_Complex8* x) ... ");
        }
    }

    template <> EXPORTCPUCOREMATH size_t amax(size_t N, const GT_Complex16* x)
    {
        try
        {
            lapack_int num = (lapack_int)(N);
            lapack_int incx = 1;

            return izamax_(&num, (lapack_complex_double*)(x), &incx);
        }
        catch(...)
        {
            GADGET_THROW("Error happened in size_t amax(size_t N, const GT_Complex16* x) ... ");
        }
    }

    /// --------------------------------------------------------------------

#ifdef USE_MKL
    template <> EXPORTCPUCOREMATH void absolute(size_t N, const float* x, float* r)
    {
        vsAbs(N, x, r);
        GADGET_CHECK_THROW(vmlGetErrStatus()==0);
    }

    template <> EXPORTCPUCOREMATH void absolute(size_t N, const double* x, double* r)
    {
        vdAbs(N, x, r);
        GADGET_CHECK_THROW(vmlGetErrStatus()==0);
    }

    template <> EXPORTCPUCOREMATH void absolute(size_t N, const GT_Complex8* x, float* r)
    {
        vcAbs(N, (MKL_Complex8*)(x), r);
        GADGET_CHECK_THROW(vmlGetErrStatus()==0);
    }

    template <> EXPORTCPUCOREMATH void absolute(size_t N, const GT_Complex16* x, double* r)
    {
        vzAbs(N, (MKL_Complex16*)(x), r);
        GADGET_CHECK_THROW(vmlGetErrStatus()==0);
    }

#else

    template <typename T> 
    void absolute(size_t N, const T* x, typename realType<T>::Type* r)
    {
        try
        {
            long long n;

            #pragma omp parallel for default(none) private(n) shared(N, x, r) if (N>NumElementsUseThreading)
            for ( n=0; n<(long long)N; n++ )
            {
                r[n]= GT_ABS(x[n]);
            }
        }
        catch(...)
        {
            GADGET_THROW("Error happened in absolute(size_t N, const T* x, typename realType<T>::Type* r) ... ");
        }
    }

    template EXPORTCPUCOREMATH void absolute(size_t N, const float* x, float* r);
    template EXPORTCPUCOREMATH void absolute(size_t N, const double* x, double* r);

    template <> EXPORTCPUCOREMATH void absolute(size_t N, const GT_Complex8* x, float* r)
    {
        try
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
        catch(...)
        {
            GADGET_THROW("Error happened in absolute(size_t N, const GT_Complex8* x, float* r) ... ");
        }
    }

    template <> EXPORTCPUCOREMATH void absolute(size_t N, const GT_Complex16* x, double* r)
    {
        try
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
        catch(...)
        {
            GADGET_THROW("Error happened in absolute(size_t N, const GT_Complex16* x, double* r) ... ");
        }
    }
#endif // USE_MKL

    template <typename T> 
    void absolute(size_t N, const std::complex<T>* x, std::complex<T>* r)
    {
        try
        {
            long long n;

            #pragma omp parallel for default(none) private(n) shared(N, x, r) if (N>NumElementsUseThreading)
            for ( n=0; n<(long long)N; n++ )
            {
                const std::complex<T>& c = x[n];
                const T re = c.real();
                const T im = c.imag();

                r[n] = std::complex<T>( std::sqrt( (re*re) + (im * im) ), 0);
            }
        }
        catch(...)
        {
            GADGET_THROW("Error happened in absolute(size_t N, const std::complex<T>* x, std::complex<T>* r) ... ");
        }
    }

    template EXPORTCPUCOREMATH void absolute(size_t N, const GT_Complex8* x, GT_Complex8* r);
    template EXPORTCPUCOREMATH void absolute(size_t N, const GT_Complex16* x, GT_Complex16* r);

    /// --------------------------------------------------------------------

#ifdef USE_MKL

    template <> EXPORTCPUCOREMATH void argument(size_t N, const GT_Complex8* x, float* r)
    {
        vcArg(N, (MKL_Complex8*)(x), r);
        GADGET_CHECK_THROW(vmlGetErrStatus()==0);
    }

    template <> EXPORTCPUCOREMATH void argument(size_t N, const GT_Complex16* x, double* r)
    {
        vzArg(N, (MKL_Complex16*)(x), r);
        GADGET_CHECK_THROW(vmlGetErrStatus()==0);
    }

#else

    template <typename T> 
    void argument(size_t N, const T* x, typename realType<T>::Type* r)
    {
        try
        {
            long long n;

            #pragma omp parallel for default(none) private(n) shared(N, x, r) if (N>NumElementsUseThreading)
            for ( n=0; n<(long long)N; n++ )
            {
                r[n] = std::arg( x[n] );
            }
        }
        catch(...)
        {
            GADGET_THROW("Error happened in argument(size_t N, const T* x, typename realType<T>::Type>* r) ... ");
        }
    }

    template EXPORTCPUCOREMATH void argument(size_t N, const GT_Complex8* x, float* r);
    template EXPORTCPUCOREMATH void argument(size_t N, const GT_Complex16* x, double* r);
#endif // USE_MKL

    /// --------------------------------------------------------------------

    template <typename T> 
    void inv(size_t N, const T* x, T* r)
    {
        try
        {
            T v(1.0);
            long long n;

            #pragma omp parallel for default(none) private(n) shared(N, x, r, v) if (N>NumElementsUseThreading)
            for ( n=0; n<(long long)N; n++ )
            {
                r[n] = v/x[n];
            }
        }
        catch(...)
        {
             GADGET_THROW("Error happened in inv(size_t N, const T* x, T* r) ... ");
        }
    }

#ifdef USE_MKL

    template <> EXPORTCPUCOREMATH void inv(size_t N, const float* x, float* r)
    {
        vsInv(N, x, r);
        GADGET_CHECK_THROW(vmlGetErrStatus()==0);
    }

    template <> EXPORTCPUCOREMATH void inv(size_t N, const double* x, double* r)
    {
        vdInv(N, x, r);
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
        try
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
        catch(...)
        {
            GADGET_THROW("Errors happened in sort(size_t N, const T* x, T* r, bool isascending) ... ");
        }
    }

    template EXPORTCPUCOREMATH void sort(size_t N, const float* x, float* r, bool isascending);
    template EXPORTCPUCOREMATH void sort(size_t N, const double* x, double* r, bool isascending);

    /// --------------------------------------------------------------------

    template<typename T> 
    void fill( size_t N, T* x, T v)
    {
        try
        {
            long long n;
            #pragma omp parallel for default(none) private(n) shared(N, x, v)
            for ( n=0; n<(long long)N; n++ )
            {
                x[n] = v;
            }
        }
        catch(...)
        {
            GADGET_THROW("Errors happened in fill( size_t N, T* x, T v) ... ");
        }
    }

    template EXPORTCPUCOREMATH void fill(size_t N, float* x, float v);
    template EXPORTCPUCOREMATH void fill(size_t N, double* x, double v);
    template EXPORTCPUCOREMATH void fill(size_t N, GT_Complex8* x, GT_Complex8 v);
    template EXPORTCPUCOREMATH void fill(size_t N, GT_Complex16* x, GT_Complex16 v);

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

            #ifdef WIN32
                x[n] = T(re*ar-im*ai, re*ai+im*ar);
            #else
                x[n].real() = (re*ar-im*ai);
                x[n].imag() = (re*ai+im*ar);
            #endif // WIN32
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

            #ifdef WIN32
                x[n] = T(re*a, im*a);
            #else
                x[n].real() = (re*a);
                x[n].imag() = (im*a);
            #endif // WIN32
        }
    }

    // -----------------------------------

    template <> EXPORTCPUCOREMATH void scal(size_t N, float a, float* x)
    {
        try
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
        catch(...)
        {
            GADGET_THROW("Error happened in scal(size_t N, float a, float* x) ... ");
        }
    }

    // -----------------------------------

    template <> EXPORTCPUCOREMATH void scal(size_t N, double a, double* x)
    {
        try
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
        catch(...)
        {
            GADGET_THROW("Error happened in scal(size_t N, double a, double* x) ... ");
        }
    }

    // -----------------------------------

    template <> EXPORTCPUCOREMATH void scal(size_t N, GT_Complex8 a, GT_Complex8* x)
    {
        try
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
        catch(...)
        {
            GADGET_THROW("Error happened in scal(size_t N, GT_Complex8 a, GT_Complex8* x) ... ");
        }
    }

    // -----------------------------------

    template <> EXPORTCPUCOREMATH void scal(size_t N, GT_Complex16 a, GT_Complex16* x)
    {
        try
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
        catch(...)
        {
            GADGET_THROW("Error happened in scal(size_t N, GT_Complex16 a, GT_Complex16* x) ... ");
        }
    }

    // -----------------------------------

    template <> EXPORTCPUCOREMATH void scal(size_t N, float a, GT_Complex8* x)
    {
        try
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
        catch(...)
        {
            GADGET_THROW("Error happened in scal(size_t N, float a, GT_Complex8* x) ... ");
        }
    }

    // -----------------------------------

    template <> EXPORTCPUCOREMATH void scal(size_t N, double a, GT_Complex16* x)
    {
        try
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
        catch(...)
        {
            GADGET_THROW("Error happened in scal(size_t N, double a, GT_Complex16* x) ... ");
        }
    }
}}
