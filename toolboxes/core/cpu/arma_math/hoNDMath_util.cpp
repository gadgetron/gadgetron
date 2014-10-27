/** \file  hoNDMath_util.cpp
\brief math function utility
*/

#include "hoNDMath_util.h"
#include <stdexcept>
#include <iostream>
#include "GadgetronCommon.h"

#ifndef lapack_int
    #define lapack_int int
#endif // lapack_int

#ifndef lapack_complex_float
    #define lapack_complex_float GT_Complex8
#endif // lapack_complex_float

#ifndef lapack_complex_double
    #define lapack_complex_double GT_Complex16
#endif // #ifndef lapack_complex_double

//Declaration of BLAS and LAPACK routines
extern "C"
{
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

#define NumElementsUseThreading 128*128*6

const int TwoGBLimit = (2147483647);

namespace Gadgetron { namespace math {

    static int get_num_threads(float computing_cost, size_t N)
    {

    }

    /// --------------------------------------------------------------------

    template <> EXPORTCPUCOREMATH void scal(size_t N, float a, float* x)
    {
        long long n;
#pragma omp parallel for private(n) if (N>NumElementsUseThreading)
        for (n = 0; n < (long long)N; n++)
        {
            x[n] *= a;
        }
    }

    // -----------------------------------

    template <> EXPORTCPUCOREMATH void scal(size_t N, double a, double* x)
    {
        long long n;
#pragma omp parallel for private(n) if (N>NumElementsUseThreading)
        for (n = 0; n < (long long)N; n++)
        {
            x[n] *= a;
        }
    }

    // -----------------------------------

    template <> EXPORTCPUCOREMATH void scal(size_t N, GT_Complex8 a, GT_Complex8* x)
    {
        long long n;

#pragma omp parallel for private(n) if (N>NumElementsUseThreading)
        for (n = 0; n < (long long)N; n++)
        {
            const GT_Complex8& c = x[n];
            const float re = c.real();
            const float im = c.imag();

            const float ar = a.real();
            const float ai = a.imag();

            reinterpret_cast<float(&)[2]>(x[n])[0] = re*ar-im*ai;
            reinterpret_cast<float(&)[2]>(x[n])[1] = re*ai+im*ar;
        }
    }

    // -----------------------------------

    template <> EXPORTCPUCOREMATH void scal(size_t N, GT_Complex16 a, GT_Complex16* x)
    {
        long long n;

#pragma omp parallel for private(n) if (N>NumElementsUseThreading)
        for (n = 0; n < (long long)N; n++)
        {
            const GT_Complex16& c = x[n];
            const double re = c.real();
            const double im = c.imag();

            const double ar = a.real();
            const double ai = a.imag();

            reinterpret_cast<double(&)[2]>(x[n])[0] = re*ar-im*ai;
            reinterpret_cast<double(&)[2]>(x[n])[1] = re*ai+im*ar;
        }
    }

    // -----------------------------------

    template <> EXPORTCPUCOREMATH void scal(size_t N, float a, GT_Complex8* x)
    {
        long long n;

#pragma omp parallel for private(n) if (N>NumElementsUseThreading)
        for (n = 0; n < (long long)N; n++)
        {
            const GT_Complex8& c = x[n];
            const float re = c.real();
            const float im = c.imag();

            reinterpret_cast<float(&)[2]>(x[n])[0] = re*a;
            reinterpret_cast<float(&)[2]>(x[n])[1] = im*a;
        }
    }

    // -----------------------------------

    template <> EXPORTCPUCOREMATH void scal(size_t N, double a, GT_Complex16* x)
    {
        long long n;

#pragma omp parallel for private(n) if (N>NumElementsUseThreading)
        for (n = 0; n < (long long)N; n++)
        {
            const GT_Complex16& c = x[n];
            const double re = c.real();
            const double im = c.imag();

            reinterpret_cast<double(&)[2]>(x[n])[0] = re*a;
            reinterpret_cast<double(&)[2]>(x[n])[1] = im*a;
        }
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
        long long n;

        #pragma omp parallel for default(none) private(n) shared(N, r, a , x, y) if(N>NumElementsUseThreading)
        for ( n=0; n<(long long)N; ++n)
        {
            r[n] = a*x[n] + y[n];
        }
    }

    template <> EXPORTCPUCOREMATH void axpy(double a, size_t N, const double* x, const double* y, double* r)
    {
        long long n;

        #pragma omp parallel for default(none) private(n) shared(N, r, a , x, y) if(N>NumElementsUseThreading)
        for ( n=0; n<(long long)N; ++n)
        {
            r[n] = a*x[n] + y[n];
        }
    }

    template <> EXPORTCPUCOREMATH void axpy(GT_Complex8 a, size_t N, const GT_Complex8* x, const GT_Complex8* y, GT_Complex8* r)
    {
        axpy_64bit_mode(a, N, x, y, r);
    }

    template <> EXPORTCPUCOREMATH void axpy(GT_Complex16 a, size_t N, const GT_Complex16* x, const GT_Complex16* y, GT_Complex16* r)
    {
        axpy_64bit_mode(a, N, x, y, r);
    }

    /// --------------------------------------------------------------------

    template <typename T> void add(size_t N, const T* x, const T* y, T* r)
    {
        long long n;
#pragma omp parallel for default(none) private(n) shared(N, x, y, r) if (N>NumElementsUseThreading)
        for ( n=0; n<(long long)N; n++ )
        {
            r[n] = x[n] + y[n];
        }
    }

    template EXPORTCPUCOREMATH void add(size_t N, const float* x, const float* y, float* r);
    template EXPORTCPUCOREMATH void add(size_t N, const double* x, const double* y, double* r);
    template EXPORTCPUCOREMATH void add(size_t N, const GT_Complex8* x, const GT_Complex8* y, GT_Complex8* r);
    template EXPORTCPUCOREMATH void add(size_t N, const GT_Complex16* x, const GT_Complex16* y, GT_Complex16* r);

    /// --------------------------------------------------------------------

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

    /// --------------------------------------------------------------------

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

    /// --------------------------------------------------------------------

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

    /// --------------------------------------------------------------------

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

    /// --------------------------------------------------------------------

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

#pragma omp parallel for private(n) if (N>NumElementsUseThreading)
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
        const double eps = std::numeric_limits<double>::epsilon();

        long long n;

#pragma omp parallel for private(n) if (N>NumElementsUseThreading)
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

        long long numOfChunk = N/TwoGBLimit;
        lapack_int len = TwoGBLimit;

        float res = 0;
        float* pX = (float*)(x);

        long long ii;
        for ( ii=0; ii<numOfChunk; ii++ )
        {
            r = snrm2_(&len, (float*)(pX), &incx);
            pX += len;
            res += r*r;
        }

        if ( (size_t)(numOfChunk*TwoGBLimit) < N )
        {
            len = (lapack_int)((long long)N - numOfChunk*(long long)TwoGBLimit);
            r = snrm2_(&len, (float*)(pX), &incx);
            res += r*r;
        }

        r = std::sqrt(r);
    }

    template <> EXPORTCPUCOREMATH void norm2(size_t N, const double* x, double& r)
    {
        lapack_int num = (lapack_int)N;
        lapack_int incx = 1;

        long long numOfChunk = N/TwoGBLimit;
        lapack_int len = TwoGBLimit;

        double res = 0;
        double* pX = (double*)(x);

        long long ii;
        for ( ii=0; ii<numOfChunk; ii++ )
        {
            r = dnrm2_(&len, pX, &incx);
            pX += len;
            res += r*r;
        }

        if ( (size_t)(numOfChunk*TwoGBLimit) < N )
        {
            len = (lapack_int)((long long)N - numOfChunk*(long long)TwoGBLimit);
            r = dnrm2_(&len, pX, &incx);
            res += r*r;
        }

        r = std::sqrt(res);
    }

    template <> EXPORTCPUCOREMATH void norm2(size_t N, const GT_Complex8* x, float& r)
    {
        lapack_int num = (lapack_int)N;
        lapack_int incx = 1;

        long long numOfChunk = N/TwoGBLimit;
        lapack_int len = TwoGBLimit;

        float res = 0;
        lapack_complex_float* pX = (lapack_complex_float*)(x);

        long long ii;
        for ( ii=0; ii<numOfChunk; ii++ )
        {
            r = scnrm2_(&len, pX, &incx);
            pX += len;
            res += r*r;
        }

        if ( (size_t)(numOfChunk*TwoGBLimit) < N )
        {
            len = (lapack_int)((long long)N - numOfChunk*(long long)TwoGBLimit);
            r = scnrm2_(&len, pX, &incx);
            res += r*r;
        }

        r = std::sqrt(res);
    }

    template <> EXPORTCPUCOREMATH void norm2(size_t N, const GT_Complex16* x, double& r)
    {
        lapack_int num = (lapack_int)N;
        lapack_int incx = 1;

        long long numOfChunk = N/TwoGBLimit;
        lapack_int len = TwoGBLimit;

        double res = 0;
        lapack_complex_double* pX = (lapack_complex_double*)(x);

        long long ii;
        for ( ii=0; ii<numOfChunk; ii++ )
        {
            r = dzasum_(&len, pX, &incx);
            pX += len;
            res += r*r;
        }

        if ( (size_t)(numOfChunk*TwoGBLimit) < N )
        {
            len = (lapack_int)((long long)N - numOfChunk*(long long)TwoGBLimit);
            r = dzasum_(&len, pX, &incx);
            res += r*r;
        }

        r = std::sqrt(res);
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
        lapack_int incx=1, incy=1;

        long long numOfChunk = N/TwoGBLimit;
        lapack_int len = TwoGBLimit;

        lapack_complex_float res(0, 0);
        r = res;

        lapack_complex_float* pX = (lapack_complex_float*)(x);
        lapack_complex_float* pY = (lapack_complex_float*)(y);

        long long ii;
        for ( ii=0; ii<numOfChunk; ii++ )
        {
            cdotc_(&res, &len, (lapack_complex_float*)(pX), &incx, (lapack_complex_float*)(pY), &incy);
            pX += len; pY += len;
            r += res;
        }

        if ( (size_t)(numOfChunk*TwoGBLimit) < N )
        {
            len = (lapack_int)((long long)N - numOfChunk*(long long)TwoGBLimit);
            cdotc_(&res, &len, (lapack_complex_float*)(pX), &incx, (lapack_complex_float*)(pY), &incy);
            r += res;
        }
    }

    template <> EXPORTCPUCOREMATH void dotc(size_t N, const GT_Complex16* x, const GT_Complex16* y, GT_Complex16& r)
    {
        lapack_int incx=1, incy=1;

        long long numOfChunk = N/TwoGBLimit;
        lapack_int len = TwoGBLimit;

        lapack_complex_double res(0, 0);
        r = res;

        lapack_complex_double* pX = (lapack_complex_double*)(x);
        lapack_complex_double* pY = (lapack_complex_double*)(y);

        long long ii;
        for ( ii=0; ii<numOfChunk; ii++ )
        {
            zdotc_(&res, &len, (lapack_complex_double*)(pX), &incx, (lapack_complex_double*)(pY), &incy);
            pX += len; pY += len;
            r += res;
        }

        if ( (size_t)(numOfChunk*TwoGBLimit) < N )
        {
            len = (lapack_int)((long long)N - numOfChunk*(long long)TwoGBLimit);
            zdotc_(&res, &len, (lapack_complex_double*)(pX), &incx, (lapack_complex_double*)(pY), &incy);
            r += res;
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

    template <> EXPORTCPUCOREMATH void dotu(size_t N, const GT_Complex8* x, const GT_Complex8* y, GT_Complex8& r)
    {
        lapack_int incx=1, incy=1;

        long long numOfChunk = N/TwoGBLimit;
        lapack_int len = TwoGBLimit;

        lapack_complex_float res(0, 0);
        r = res;

        lapack_complex_float* pX = (lapack_complex_float*)(x);
        lapack_complex_float* pY = (lapack_complex_float*)(y);

        long long ii;
        for ( ii=0; ii<numOfChunk; ii++ )
        {
            cdotu_(&res, &len, (lapack_complex_float*)(pX), &incx, (lapack_complex_float*)(pY), &incy);
            pX += len; pY += len;
            r += res;
        }

        if ( (size_t)(numOfChunk*TwoGBLimit) < N )
        {
            len = (lapack_int)((long long)N - numOfChunk*(long long)TwoGBLimit);
            cdotu_(&res, &len, (lapack_complex_float*)(pX), &incx, (lapack_complex_float*)(pY), &incy);
            r += res;
        }
    }

    template <> EXPORTCPUCOREMATH void dotu(size_t N, const GT_Complex16* x, const GT_Complex16* y, GT_Complex16& r)
    {
        lapack_int incx=1, incy=1;

        long long numOfChunk = N/TwoGBLimit;
        lapack_int len = TwoGBLimit;

        lapack_complex_double res(0, 0);
        r = res;

        lapack_complex_double* pX = (lapack_complex_double*)(x);
        lapack_complex_double* pY = (lapack_complex_double*)(y);

        long long ii;
        for ( ii=0; ii<numOfChunk; ii++ )
        {
            zdotu_(&res, &len, (lapack_complex_double*)(pX), &incx, (lapack_complex_double*)(pY), &incy);
            pX += len; pY += len;
            r += res;
        }

        if ( (size_t)(numOfChunk*TwoGBLimit) < N )
        {
            len = (lapack_int)((long long)N - numOfChunk*(long long)TwoGBLimit);
            zdotu_(&res, &len, (lapack_complex_double*)(pX), &incx, (lapack_complex_double*)(pY), &incy);
            r += res;
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
        lapack_int incx=1;

        long long numOfChunk = N/TwoGBLimit;
        lapack_int len = TwoGBLimit;

        float res(0);
        r = res;

        float* pX = (float*)(x);

        long long ii;
        for ( ii=0; ii<numOfChunk; ii++ )
        {
            res = sasum_(&len, pX, &incx);
            pX += len;
            r += res;
        }

        if ( (size_t)(numOfChunk*TwoGBLimit) < N )
        {
            len = (lapack_int)((long long)N - numOfChunk*(long long)TwoGBLimit);
            res = sasum_(&len, pX, &incx);
            r += res;
        }
    }

    template <> EXPORTCPUCOREMATH void asum(size_t N, const double* x, double& r)
    {
        lapack_int incx=1;

        long long numOfChunk = N/TwoGBLimit;
        lapack_int len = TwoGBLimit;

        double res(0);
        r = res;

        double* pX = (double*)(x);

        long long ii;
        for ( ii=0; ii<numOfChunk; ii++ )
        {
            res = dasum_(&len, pX, &incx);
            pX += len;
            r += res;
        }

        if ( (size_t)(numOfChunk*TwoGBLimit) < N )
        {
            len = (lapack_int)((long long)N - numOfChunk*(long long)TwoGBLimit);
            res = dasum_(&len, pX, &incx);
            r += res;
        }
    }

    template <> EXPORTCPUCOREMATH void asum(size_t N, const GT_Complex8* x, float& r)
    {
        lapack_int incx=1;

        long long numOfChunk = N/TwoGBLimit;
        lapack_int len = TwoGBLimit;

        float res(0);
        r = res;

        lapack_complex_float* pX = (lapack_complex_float*)(x);

        long long ii;
        for ( ii=0; ii<numOfChunk; ii++ )
        {
            res = scasum_(&len, (lapack_complex_float*)(pX), &incx);
            pX += len;
            r += res;
        }

        if ( (size_t)(numOfChunk*TwoGBLimit) < N )
        {
            len = (lapack_int)((long long)N - numOfChunk*(long long)TwoGBLimit);
            res = scasum_(&len, (lapack_complex_float*)(pX), &incx);
            r += res;
        }
    }

    template <> EXPORTCPUCOREMATH void asum(size_t N, const GT_Complex16* x, double& r)
    {
        lapack_int incx=1;

        long long numOfChunk = N/TwoGBLimit;
        lapack_int len = TwoGBLimit;

        double res(0);
        r = res;

        lapack_complex_double* pX = (lapack_complex_double*)(x);

        long long ii;
        for ( ii=0; ii<numOfChunk; ii++ )
        {
            res = dzasum_(&len, (lapack_complex_double*)(pX), &incx);
            pX += len;
            r += res;
        }

        if ( (size_t)(numOfChunk*TwoGBLimit) < N )
        {
            len = (lapack_int)((long long)N - numOfChunk*(long long)TwoGBLimit);
            res = dzasum_(&len, (lapack_complex_double*)(pX), &incx);
            r += res;
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

    template EXPORTCPUCOREMATH void inv(size_t N, const float* x, float* r);
    template EXPORTCPUCOREMATH void inv(size_t N, const double* x, double* r);
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
