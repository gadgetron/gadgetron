#include "hoNDArray_math_util.h"

#ifndef lapack_int
    #define lapack_int int
#endif // lapack_int

#ifndef lapack_complex_float
    #define lapack_complex_float  std::complex<float> 
#endif // lapack_complex_float

#ifndef lapack_complex_double
    #define lapack_complex_double  std::complex<double> 
#endif // #ifndef lapack_complex_double

//Declaration of BLAS and LAPACK routines
extern "C"
{
    /// Finds the index of the element with the maximal absolute value.
    lapack_int isamax_(lapack_int* N, float* x, lapack_int* incx);
    lapack_int idamax_(lapack_int* N, double* x, lapack_int* incx);
    lapack_int icamax_(lapack_int* N, lapack_complex_float* x, lapack_int* incx);
    lapack_int izamax_(lapack_int* N, lapack_complex_double* x, lapack_int* incx);
}

#define NumElementsUseThreading 64*1024

namespace Gadgetron
{
    // --------------------------------------------------------------------------------

    inline void add(size_t N, const float* x, const float* y, float* r)
    {
        long long n;

        #pragma omp parallel for default(none) private(n) shared(N, r, x, y) if(N>NumElementsUseThreading)
        for ( n=0; n<(long long)N; ++n)
        {
            r[n] = x[n] + y[n];
        }
    }

    inline void add(size_t N, const double* x, const double* y, double* r)
    {
        long long n;

        #pragma omp parallel for default(none) private(n) shared(N, r, x, y) if(N>NumElementsUseThreading)
        for ( n=0; n<(long long)N; ++n)
        {
            r[n] = x[n] + y[n];
        }
    }

    inline void add(size_t N, const  std::complex<float> * x, const  std::complex<float> * y,  std::complex<float> * r)
    {
        long long n;

        #pragma omp parallel for default(none) private(n) shared(N, r, x, y) if(N>NumElementsUseThreading)
        for ( n=0; n<(long long)N; ++n)
        {
            const  std::complex<float> & vx = x[n];
            const float re1 = vx.real();
            const float im1 = vx.imag();

            const  std::complex<float> & vy = y[n];
            const float re2 = vy.real();
            const float im2 = vy.imag();

            reinterpret_cast<float(&)[2]>(r[n])[0] = re1 + re2;
            reinterpret_cast<float(&)[2]>(r[n])[1] = im1 + im2;
        }
    }

    inline void add(size_t N, const  std::complex<double> * x, const  std::complex<double> * y,  std::complex<double> * r)
    {
        long long n;

        #pragma omp parallel for default(none) private(n) shared(N, r, x, y) if(N>NumElementsUseThreading)
        for ( n=0; n<(long long)N; ++n)
        {
            const  std::complex<double> & vx = x[n];
            const double re1 = vx.real();
            const double im1 = vx.imag();

            const  std::complex<double> & vy = y[n];
            const double re2 = vy.real();
            const double im2 = vy.imag();

            reinterpret_cast<double(&)[2]>(r[n])[0] = re1 + re2;
            reinterpret_cast<double(&)[2]>(r[n])[1] = im1 + im2;
        }
    }

    template <typename T> 
    bool add(const hoNDArray<T>& x, const hoNDArray<T>& y, hoNDArray<T>& r)
    {
        GADGET_DEBUG_CHECK_RETURN_FALSE(x.get_number_of_elements()==y.get_number_of_elements());
        if ( r.get_number_of_elements()!=x.get_number_of_elements())
        {
            r = x;
        }

        add(x.get_number_of_elements(), x.begin(), y.begin(), r.begin());

        return true;
    }

    template EXPORTCPUCOREMATH bool add(const hoNDArray<float>& x, const hoNDArray<float>& y, hoNDArray<float>& r);
    template EXPORTCPUCOREMATH bool add(const hoNDArray<double>& x, const hoNDArray<double>& y, hoNDArray<double>& r);
    template EXPORTCPUCOREMATH bool add(const hoNDArray< std::complex<float> >& x, const hoNDArray< std::complex<float> >& y, hoNDArray< std::complex<float> >& r);
    template EXPORTCPUCOREMATH bool add(const hoNDArray< std::complex<double> >& x, const hoNDArray< std::complex<double> >& y, hoNDArray< std::complex<double> >& r);

    // --------------------------------------------------------------------------------

    inline void subtract(size_t N, const float* x, const float* y, float* r)
    {
        long long n;

        #pragma omp parallel for default(none) private(n) shared(N, r, x, y) if(N>NumElementsUseThreading)
        for ( n=0; n<(long long)N; ++n)
        {
            r[n] = x[n] - y[n];
        }
    }

    inline void subtract(size_t N, const double* x, const double* y, double* r)
    {
        long long n;

        #pragma omp parallel for default(none) private(n) shared(N, r, x, y) if(N>NumElementsUseThreading)
        for ( n=0; n<(long long)N; ++n)
        {
            r[n] = x[n] - y[n];
        }
    }

    inline void subtract(size_t N, const  std::complex<float> * x, const  std::complex<float> * y,  std::complex<float> * r)
    {
        long long n;

        #pragma omp parallel for default(none) private(n) shared(N, r, x, y) if(N>NumElementsUseThreading)
        for ( n=0; n<(long long)N; ++n)
        {
            const  std::complex<float> & vx = x[n];
            const float re1 = vx.real();
            const float im1 = vx.imag();

            const  std::complex<float> & vy = y[n];
            const float re2 = vy.real();
            const float im2 = vy.imag();

            reinterpret_cast<float(&)[2]>(r[n])[0] = re1 - re2;
            reinterpret_cast<float(&)[2]>(r[n])[1] = im1 - im2;
        }
    }

    inline void subtract(size_t N, const  std::complex<double> * x, const  std::complex<double> * y,  std::complex<double> * r)
    {
        long long n;

        #pragma omp parallel for default(none) private(n) shared(N, r, x, y) if(N>NumElementsUseThreading)
        for ( n=0; n<(long long)N; ++n)
        {
            const  std::complex<double> & vx = x[n];
            const double re1 = vx.real();
            const double im1 = vx.imag();

            const  std::complex<double> & vy = y[n];
            const double re2 = vy.real();
            const double im2 = vy.imag();

            reinterpret_cast<double(&)[2]>(r[n])[0] = re1 - re2;
            reinterpret_cast<double(&)[2]>(r[n])[1] = im1 - im2;
        }
    }

    template <typename T> 
    bool subtract(const hoNDArray<T>& x, const hoNDArray<T>& y, hoNDArray<T>& r)
    {
        GADGET_DEBUG_CHECK_RETURN_FALSE(x.get_number_of_elements()==y.get_number_of_elements());
        if ( r.get_number_of_elements()!=x.get_number_of_elements())
        {
            r = x;
        }

        subtract(x.get_number_of_elements(), x.begin(), y.begin(), r.begin());

        return true;
    }

    template EXPORTCPUCOREMATH bool subtract(const hoNDArray<float>& x, const hoNDArray<float>& y, hoNDArray<float>& r);
    template EXPORTCPUCOREMATH bool subtract(const hoNDArray<double>& x, const hoNDArray<double>& y, hoNDArray<double>& r);
    template EXPORTCPUCOREMATH bool subtract(const hoNDArray< std::complex<float> >& x, const hoNDArray< std::complex<float> >& y, hoNDArray< std::complex<float> >& r);
    template EXPORTCPUCOREMATH bool subtract(const hoNDArray< std::complex<double> >& x, const hoNDArray< std::complex<double> >& y, hoNDArray< std::complex<double> >& r);

    // --------------------------------------------------------------------------------

    template <typename T> 
    inline void multiply(size_t N, const T* x, const T* y, T* r)
    {
        long long n;
        #pragma omp parallel for default(none) private(n) shared(N, x, y, r) if (N>NumElementsUseThreading)
        for ( n=0; n<(long long)N; n++ )
        {
            const T& a = x[n];
            const T& b = y[n];
            r[n] = a*b;
        }
    }

    inline void multiply(size_t N, const  std::complex<float> * x, const  std::complex<float> * y,  std::complex<float> * r)
    {
        long long n;
        #pragma omp parallel for default(none) private(n) shared(N, x, y, r) if (N>NumElementsUseThreading)
        for (n = 0; n < (long long)N; n++)
        {
            const std::complex<float>& a1 = x[n];
            const std::complex<float>& b1 = y[n];
            const float a = a1.real();
            const float b = a1.imag();
            const float c = b1.real();
            const float d = b1.imag();

            reinterpret_cast<float(&)[2]>(r[n])[0] = a*c-b*d;
            reinterpret_cast<float(&)[2]>(r[n])[1] = a*d+b*c;
        }
    }

    inline void multiply(size_t N, const  std::complex<double> * x, const  std::complex<double> * y,  std::complex<double> * r)
    {
        long long n;
        #pragma omp parallel for default(none) private(n) shared(N, x, y, r) if (N>NumElementsUseThreading)
        for (n = 0; n < (long long)N; n++)
        {
            const std::complex<double>& a1 = x[n];
            const std::complex<double>& b1 = y[n];
            const double a = a1.real();
            const double b = a1.imag();
            const double c = b1.real();
            const double d = b1.imag();

            reinterpret_cast<double(&)[2]>(r[n])[0] = a*c-b*d;
            reinterpret_cast<double(&)[2]>(r[n])[1] = a*d+b*c;
        }
    }

    template <typename T> 
    bool multiply(const hoNDArray<T>& x, const hoNDArray<T>& y, hoNDArray<T>& r)
    {
        GADGET_DEBUG_CHECK_RETURN_FALSE(x.get_number_of_elements()==y.get_number_of_elements());
        if ( r.get_number_of_elements()!=x.get_number_of_elements())
        {
            r = x;
        }

        multiply(x.get_number_of_elements(), x.begin(), y.begin(), r.begin());

        return true;
    }

    template EXPORTCPUCOREMATH bool multiply(const hoNDArray<float>& x, const hoNDArray<float>& y, hoNDArray<float>& r);
    template EXPORTCPUCOREMATH bool multiply(const hoNDArray<double>& x, const hoNDArray<double>& y, hoNDArray<double>& r);
    template EXPORTCPUCOREMATH bool multiply(const hoNDArray< std::complex<float> >& x, const hoNDArray< std::complex<float> >& y, hoNDArray< std::complex<float> >& r);
    template EXPORTCPUCOREMATH bool multiply(const hoNDArray< std::complex<double> >& x, const hoNDArray< std::complex<double> >& y, hoNDArray< std::complex<double> >& r);

    // --------------------------------------------------------------------------------

    template <typename T> 
    inline void divide(size_t N, const T* x, const T* y, T* r)
    {
        long long n;
        #pragma omp parallel for default(none) private(n) shared(N, x, y, r) if (N>NumElementsUseThreading)
        for ( n=0; n<(long long)N; n++ )
        {
            const T& a = x[n];
            const T& b = y[n];
            r[n] = a/b;
        }
    }

    inline void divide(size_t N, const  std::complex<float> * x, const  std::complex<float> * y,  std::complex<float> * r)
    {
        long long n;
        #pragma omp parallel for default(none) private(n) shared(N, x, y, r) if (N>NumElementsUseThreading)
        for (n = 0; n < (long long)N; n++)
        {
            const std::complex<float>& a1 = x[n];
            const std::complex<float>& b1 = y[n];
            const float a = a1.real();
            const float b = a1.imag();
            const float c = b1.real();
            const float d = b1.imag();

            const float m = 1/(c*c+d*d);

            reinterpret_cast<float(&)[2]>(r[n])[0] = (a*c+b*d)*m;
            reinterpret_cast<float(&)[2]>(r[n])[1] = (b*c-a*d)*m;
        }
    }

    inline void divide(size_t N, const  std::complex<double> * x, const  std::complex<double> * y,  std::complex<double> * r)
    {
        long long n;
        #pragma omp parallel for default(none) private(n) shared(N, x, y, r) if (N>NumElementsUseThreading)
        for (n = 0; n < (long long)N; n++)
        {
            const std::complex<double>& a1 = x[n];
            const std::complex<double>& b1 = y[n];
            const double a = a1.real();
            const double b = a1.imag();
            const double c = b1.real();
            const double d = b1.imag();

            const double m = 1/(c*c+d*d);

            reinterpret_cast<double(&)[2]>(r[n])[0] = (a*c+b*d)*m;
            reinterpret_cast<double(&)[2]>(r[n])[1] = (b*c-a*d)*m;
        }
    }

    template <typename T> 
    bool divide(const hoNDArray<T>& x, const hoNDArray<T>& y, hoNDArray<T>& r)
    {
        GADGET_DEBUG_CHECK_RETURN_FALSE(x.get_number_of_elements()==y.get_number_of_elements());
        if ( r.get_number_of_elements()!=x.get_number_of_elements())
        {
            r = x;
        }

        divide(x.get_number_of_elements(), x.begin(), y.begin(), r.begin());

        return true;
    }

    template EXPORTCPUCOREMATH bool divide(const hoNDArray<float>& x, const hoNDArray<float>& y, hoNDArray<float>& r);
    template EXPORTCPUCOREMATH bool divide(const hoNDArray<double>& x, const hoNDArray<double>& y, hoNDArray<double>& r);
    template EXPORTCPUCOREMATH bool divide(const hoNDArray< std::complex<float> >& x, const hoNDArray< std::complex<float> >& y, hoNDArray< std::complex<float> >& r);
    template EXPORTCPUCOREMATH bool divide(const hoNDArray< std::complex<double> >& x, const hoNDArray< std::complex<double> >& y, hoNDArray< std::complex<double> >& r);

    // --------------------------------------------------------------------------------

    template <typename T> 
    bool sqrt(const hoNDArray<T>& x, hoNDArray<T>& r)
    {
        if ( r.get_number_of_elements()!=x.get_number_of_elements())
        {
            r.create(x.get_dimensions());
        }

        size_t N = x.get_number_of_elements();
        const T* pX = x.begin();
        T* pR = r.begin();

        long long n;
        #pragma omp parallel for default(none) private(n) shared(N, pX, pR) if (N>NumElementsUseThreading)
        for ( n=0; n<(long long)N; n++ )
        {
            pR[n] = std::sqrt(pX[n]);
        }

        return true;
    }

    template EXPORTCPUCOREMATH bool sqrt(const hoNDArray<float>& x, hoNDArray<float>& r);
    template EXPORTCPUCOREMATH bool sqrt(const hoNDArray<double>& x, hoNDArray<double>& r);
    template EXPORTCPUCOREMATH bool sqrt(const hoNDArray< std::complex<float> >& x, hoNDArray< std::complex<float> >& r);
    template EXPORTCPUCOREMATH bool sqrt(const hoNDArray< std::complex<double> >& x, hoNDArray< std::complex<double> >& r);

    // --------------------------------------------------------------------------------

    template <typename T> 
    bool minAbsolute(const hoNDArray<T>& x, T& r, size_t& ind)
    {
        size_t N = x.get_number_of_elements();
        const T* pX = x.begin();

        ind = 0;
        if ( N == 0 ) return true;

        long long n;

        typename realType<T>::Type v = abs(pX[0]);
        typename realType<T>::Type v2;

        ind = 0;
        for ( n=1; n<(long long)N; n++ )
        {
            v2 = std::abs(pX[n]);
            if ( v2 < v )
            {
                v = v2;
                ind = n;
            }
        }

        r = pX[ind];

        return true;
    }

    template EXPORTCPUCOREMATH bool minAbsolute(const hoNDArray<float>& x, float& r, size_t& ind);
    template EXPORTCPUCOREMATH bool minAbsolute(const hoNDArray<double>& x, double& r, size_t& ind);
    template EXPORTCPUCOREMATH bool minAbsolute(const hoNDArray< std::complex<float> >& x,  std::complex<float> & r, size_t& ind);
    template EXPORTCPUCOREMATH bool minAbsolute(const hoNDArray< std::complex<double> >& x,  std::complex<double> & r, size_t& ind);

    // --------------------------------------------------------------------------------

    template <typename T> 
    bool maxAbsolute(const hoNDArray<T>& x, T& r, size_t& ind)
    {
        size_t N = x.get_number_of_elements();
        const T* pX = x.begin();

        ind = 0;
        if ( N == 0 ) return true;

        long long n;

        typename realType<T>::Type v = abs(pX[0]);
        typename realType<T>::Type v2;

        ind = 0;
        for ( n=1; n<(long long)N; n++ )
        {
            v2 = std::abs(pX[n]);
            if ( v2 > v )
            {
                v = v2;
                ind = n;
            }
        }

        r = pX[ind];

        return true;
    }

    template EXPORTCPUCOREMATH bool maxAbsolute(const hoNDArray<float>& x, float& r, size_t& ind);
    template EXPORTCPUCOREMATH bool maxAbsolute(const hoNDArray<double>& x, double& r, size_t& ind);
    template EXPORTCPUCOREMATH bool maxAbsolute(const hoNDArray< std::complex<float> >& x,  std::complex<float> & r, size_t& ind);
    template EXPORTCPUCOREMATH bool maxAbsolute(const hoNDArray< std::complex<double> >& x,  std::complex<double> & r, size_t& ind);

    // --------------------------------------------------------------------------------

    inline void multiplyConj(size_t N, const  std::complex<float> * x, const  std::complex<float> * y,  std::complex<float> * r)
    {
        long long n;

        #pragma omp parallel for default(none) private(n) shared(N, x, y, r) if (N>NumElementsUseThreading)
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

    inline void multiplyConj(size_t N, const  std::complex<double> * x, const  std::complex<double> * y,  std::complex<double> * r)
    {
        long long n;

        #pragma omp parallel for default(none) private(n) shared(N, x, y, r) if (N>NumElementsUseThreading)
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

    template <typename T> 
    bool multiplyConj(const hoNDArray<T>& x, const hoNDArray<T>& y, hoNDArray<T>& r)
    {
        GADGET_DEBUG_CHECK_RETURN_FALSE(x.get_number_of_elements()==y.get_number_of_elements());
        if ( r.get_number_of_elements()!=x.get_number_of_elements())
        {
            r = x;
        }

        multiplyConj(x.get_number_of_elements(), x.begin(), y.begin(), r.begin());

        return true;
    }

    template EXPORTCPUCOREMATH bool multiplyConj(const hoNDArray< std::complex<float> >& x, const hoNDArray< std::complex<float> >& y, hoNDArray< std::complex<float> >& r);
    template EXPORTCPUCOREMATH bool multiplyConj(const hoNDArray< std::complex<double> >& x, const hoNDArray< std::complex<double> >& y, hoNDArray< std::complex<double> >& r);

    // --------------------------------------------------------------------------------

    inline void conjugate(size_t N, const  std::complex<float> * x,  std::complex<float> * r)
    {
        long long n;

        #pragma omp parallel for default(none) private(n) shared(N, x, r) if (N>NumElementsUseThreading)
        for ( n=0; n<(long long)N; n++ )
        {
            reinterpret_cast<float(&)[2]>(r[n])[0] = reinterpret_cast< const float(&)[2]>(x[n])[0];
            reinterpret_cast<float(&)[2]>(r[n])[1] = -(reinterpret_cast< const float(&)[2]>(x[n])[1]);
        }
    }

    inline void conjugate(size_t N, const  std::complex<double> * x,  std::complex<double> * r)
    {
        long long n;

        #pragma omp parallel for default(none) private(n) shared(N, x, r) if (N>NumElementsUseThreading)
        for ( n=0; n<(long long)N; n++ )
        {
            reinterpret_cast<double(&)[2]>(r[n])[0] = reinterpret_cast< const double(&)[2]>(x[n])[0];
            reinterpret_cast<double(&)[2]>(r[n])[1] = -(reinterpret_cast<const double(&)[2]>(x[n])[1]);
        }
    }

    template <typename T> 
    bool conjugate(const hoNDArray<T>& x, hoNDArray<T>& r)
    {
        if ( r.get_number_of_elements()!=x.get_number_of_elements())
        {
            r.create(x.get_dimensions());
        }

        conjugate(x.get_number_of_elements(), x.begin(), r.begin());

        return true;
    }

    template EXPORTCPUCOREMATH bool conjugate(const hoNDArray< std::complex<float> >& x, hoNDArray< std::complex<float> >& r);
    template EXPORTCPUCOREMATH bool conjugate(const hoNDArray< std::complex<double> >& x, hoNDArray< std::complex<double> >& r);

    // --------------------------------------------------------------------------------

    template <typename T> 
    inline void addEpsilon(size_t N, T* x)
    {
        typename realType<T>::Type eps = std::numeric_limits<typename realType<T>::Type>::epsilon();

        long long n;

        #pragma omp parallel for default(none) private(n) shared(N, x, eps) if (N>NumElementsUseThreading)
        for (n=0; n<(long long)N; n++ )
        {
            if ( std::abs(x[n]) < eps )
            {
                x[n] += eps;
            }
        }
    }

    inline void addEpsilon(size_t N,  std::complex<float> * x)
    {
        const float eps = std::numeric_limits<float>::epsilon();

        long long n;

        #pragma omp parallel for private(n) if (N>NumElementsUseThreading)
        for (n=0; n<(long long)N; n++ )
        {
            if ( std::abs(x[n]) < eps )
            {
                reinterpret_cast<float(&)[2]>(x[n])[0] += eps;
            }
        }
    }

    inline void addEpsilon(size_t N,  std::complex<double> * x)
    {
        const double eps = std::numeric_limits<double>::epsilon();

        long long n;

        #pragma omp parallel for private(n) if (N>NumElementsUseThreading)
        for (n=0; n<(long long)N; n++ )
        {
            if ( std::abs(x[n]) < eps )
            {
                reinterpret_cast<double(&)[2]>(x[n])[0] += eps;
            }
        }
    }

    template <typename T> 
    bool addEpsilon(hoNDArray<T>& x)
    {
        addEpsilon(x.get_number_of_elements(), x.begin());
        return true;
    }

    template EXPORTCPUCOREMATH bool addEpsilon(hoNDArray<float>& x);
    template EXPORTCPUCOREMATH bool addEpsilon(hoNDArray<double>& x);
    template EXPORTCPUCOREMATH bool addEpsilon(hoNDArray< std::complex<float> >& x);
    template EXPORTCPUCOREMATH bool addEpsilon(hoNDArray< std::complex<double> >& x);

    // --------------------------------------------------------------------------------

    inline void norm2(size_t N, const float* x, float& r)
    {
        long long i;

        float sum(0);

        #pragma omp parallel for private(i) reduction(+:sum) if (N>NumElementsUseThreading)
        for (i = 0; i < (long long)N; i++)
        {
            const float& re = x[i];
            sum += ( re*re );
        }

        r = std::sqrt(sum);
    }

    inline void norm2(size_t N, const double* x, double& r)
    {
        long long i;

        double sum(0);

        #pragma omp parallel for private(i) reduction(+:sum) if (N>NumElementsUseThreading)
        for (i = 0; i < (long long)N; i++)
        {
            const double& re = x[i];
            sum += ( re*re );
        }

        r = std::sqrt(sum);
    }

    inline void norm2(size_t N, const  std::complex<float> * x, float& r)
    {
        long long i;

        float sum(0);

        #pragma omp parallel for private(i) reduction(+:sum) if (N>NumElementsUseThreading)
        for (i = 0; i < (long long)N; i++)
        {
            const std::complex<float>& c = x[i];
            const float re = c.real();
            const float im = c.imag();
            sum += ( (re*re) + (im * im) );
        }

        r = std::sqrt(sum);
    }

    inline void norm2(size_t N, const  std::complex<double> * x, double& r)
    {
        long long i;

        double sum(0);

        #pragma omp parallel for private(i) reduction(+:sum) if (N>NumElementsUseThreading)
        for (i = 0; i < (long long)N; i++)
        {
            const std::complex<double>& c = x[i];
            const double re = c.real();
            const double im = c.imag();
            sum += ( (re*re) + (im * im) );
        }

        r = std::sqrt(sum);
    }

    template <typename T> 
    bool norm2(const hoNDArray<T>& x, typename realType<T>::Type& r)
    {
        norm2(x.get_number_of_elements(), x.begin(), r);
        return true;
    }

    template EXPORTCPUCOREMATH bool norm2(const hoNDArray<float>& x, float& r);
    template EXPORTCPUCOREMATH bool norm2(const hoNDArray<double>& x, double& r);
    template EXPORTCPUCOREMATH bool norm2(const hoNDArray< std::complex<float> >& x, float& r);
    template EXPORTCPUCOREMATH bool norm2(const hoNDArray< std::complex<double> >& x, double& r);

    template <typename T> inline 
    typename realType<T>::Type norm2(const hoNDArray<T>& x)
    {
        typename realType<T>::Type r;
        norm2(x, r);
        return r;
    }

    template EXPORTCPUCOREMATH float norm2(const hoNDArray<float>& x);
    template EXPORTCPUCOREMATH double norm2(const hoNDArray<double>& x);
    template EXPORTCPUCOREMATH float norm2(const hoNDArray< std::complex<float> >& x);
    template EXPORTCPUCOREMATH double norm2(const hoNDArray< std::complex<double> >& x);

    // --------------------------------------------------------------------------------

    template <typename T> inline 
    void norm1(size_t N, const T* x, typename realType<T>::Type& r)
    {
        long long n;

        typename realType<T>::Type norm1Sum(0);

        #pragma omp parallel for private(n) reduction(+:norm1Sum) if (N>NumElementsUseThreading)
        for (n=0; n<(long long)N; n++)
        {
            const T& c = x[n];
            norm1Sum += GT_ABS(c);
        }

        r = norm1Sum;
    }

    inline void norm1(size_t N, const  std::complex<float> * x, float& r)
    {
        long long i;
        float sum = 0.0f;
        #pragma omp parallel for private(i) reduction(+:sum) if (N>NumElementsUseThreading)
        for (i = 0; i < (long long)N; i++)
        {
            const std::complex<float>& c = x[i];
            const float re = c.real();
            const float im = c.imag();
            sum += std::sqrt( (re*re) + (im * im) );
        }

        r = sum;
    }

    inline void norm1(size_t N, const  std::complex<double> * x, double& r)
    {
        long long i;
        double sum = 0.0;
        #pragma omp parallel for private(i) reduction(+:sum) if (N>NumElementsUseThreading)
        for (i = 0; i < (long long)N; i++)
        {
            const std::complex<double>& c = x[i];
            const double re = c.real();
            const double im = c.imag();
            sum += std::sqrt( (re*re) + (im * im) );
        }

        r = sum;
    }

    template <typename T> 
    bool norm1(const hoNDArray<T>& x, typename realType<T>::Type& r)
    {
        norm1(x.get_number_of_elements(), x.begin(), r);
        return true;
    }

    template EXPORTCPUCOREMATH bool norm1(const hoNDArray<float>& x, float& r);
    template EXPORTCPUCOREMATH bool norm1(const hoNDArray<double>& x, double& r);
    template EXPORTCPUCOREMATH bool norm1(const hoNDArray< std::complex<float> >& x, float& r);
    template EXPORTCPUCOREMATH bool norm1(const hoNDArray< std::complex<double> >& x, double& r);

    template <typename T> inline 
    typename realType<T>::Type norm1(const hoNDArray<T>& x)
    {
        typename realType<T>::Type r;
        norm1(x, r);
        return r;
    }

    template EXPORTCPUCOREMATH float norm1(const hoNDArray<float>& x);
    template EXPORTCPUCOREMATH double norm1(const hoNDArray<double>& x);
    template EXPORTCPUCOREMATH float norm1(const hoNDArray< std::complex<float> >& x);
    template EXPORTCPUCOREMATH double norm1(const hoNDArray< std::complex<double> >& x);

    // --------------------------------------------------------------------------------

    inline void dotc(size_t N, const  std::complex<float> * x, const  std::complex<float> * y,  std::complex<float> & r)
    {
        long long n;

        float sum(0);

        float sa(0), sb(0);

        #pragma omp parallel for private(n) reduction(+:sa) if (N>NumElementsUseThreading)
        for (n = 0; n < (long long)N; n++)
        {
            const float a = x[n].real();
            const float b = x[n].imag();
            const float c = y[n].real();
            const float d = y[n].imag();

            sa += (a*c + b*d);
            sb += (c*b - a*d);
        }

        reinterpret_cast<float(&)[2]>(r)[0] = sa;
        reinterpret_cast<float(&)[2]>(r)[1] = sb;
    }

    inline void dotc(size_t N, const  std::complex<double> * x, const  std::complex<double> * y,  std::complex<double> & r)
    {
        long long n;

        double sum(0);

        double sa(0), sb(0);

        #pragma omp parallel for private(n) reduction(+:sa) if (N>NumElementsUseThreading)
        for (n = 0; n < (long long)N; n++)
        {
            const double a = x[n].real();
            const double b = x[n].imag();
            const double c = y[n].real();
            const double d = y[n].imag();

            sa += (a*c + b*d);
            sb += (c*b - a*d);
        }

        reinterpret_cast<double(&)[2]>(r)[0] = sa;
        reinterpret_cast<double(&)[2]>(r)[1] = sb;
    }

    template <typename T> 
    bool dotc(const hoNDArray<T>& x, const hoNDArray<T>& y, T& r)
    {
        GADGET_DEBUG_CHECK_RETURN_FALSE(x.get_number_of_elements()==y.get_number_of_elements());
        dotc(x.get_number_of_elements(), x.begin(), y.begin(), r);
        return true;
    }

    template EXPORTCPUCOREMATH bool dotc(const hoNDArray< std::complex<float> >& x, const hoNDArray< std::complex<float> >& y,  std::complex<float> & r);
    template EXPORTCPUCOREMATH bool dotc(const hoNDArray< std::complex<double> >& x, const hoNDArray< std::complex<double> >& y,  std::complex<double> & r);

    template <typename T> 
    T dotc(const hoNDArray<T>& x, const hoNDArray<T>& y)
    {
        T r;
        dotc(x, y, r);
        return r;
    }

    template EXPORTCPUCOREMATH std::complex<float> dotc(const hoNDArray< std::complex<float> >& x, const hoNDArray< std::complex<float> >& y);
    template EXPORTCPUCOREMATH std::complex<double> dotc(const hoNDArray< std::complex<double> >& x, const hoNDArray< std::complex<double> >& y);

    // --------------------------------------------------------------------------------

    inline void dotu(size_t N, const float* x, const float* y, float& r)
    {
        long long n;

        float res(0);

        #pragma omp parallel for private(n) reduction(+:res) if (N>NumElementsUseThreading)
        for (n=0; n<(long long)N; n++)
        {
            res += x[n]*y[n];
        }

        r = res;
    }

    inline void dotu(size_t N, const double* x, const double* y, double& r)
    {
        long long n;

        double res(0);

        #pragma omp parallel for private(n) reduction(+:res) if (N>NumElementsUseThreading)
        for (n=0; n<(long long)N; n++)
        {
            res += x[n]*y[n];
        }

        r = res;
    }

    inline void dotu(size_t N, const  std::complex<float> * x, const  std::complex<float> * y,  std::complex<float> & r)
    {
        long long n;

         std::complex<float>  sum(0);

        float sa(0), sb(0);
        #pragma omp parallel for private(n) reduction(+:sa) reduction(+:sb) if (N>NumElementsUseThreading)
        for (n = 0; n < (long long)N; n++)
        {
            const float a = x[n].real();
            const float b = x[n].imag();
            const float c = y[n].real();
            const float d = y[n].imag();

            sa += (a*c - b*d);
            sb += (c*b + a*d);
        }

        reinterpret_cast<float(&)[2]>(r)[0] = sa;
        reinterpret_cast<float(&)[2]>(r)[1] = sb;
    }

    inline void dotu(size_t N, const  std::complex<double> * x, const  std::complex<double> * y,  std::complex<double> & r)
    {
        long long n;

         std::complex<double>  sum(0);

        double sa(0), sb(0);
        #pragma omp parallel for private(n) reduction(+:sa) reduction(+:sb) if (N>NumElementsUseThreading)
        for (n = 0; n < (long long)N; n++)
        {
            const double a = x[n].real();
            const double b = x[n].imag();
            const double c = y[n].real();
            const double d = y[n].imag();

            sa += (a*c - b*d);
            sb += (c*b + a*d);
        }

        reinterpret_cast<double(&)[2]>(r)[0] = sa;
        reinterpret_cast<double(&)[2]>(r)[1] = sb;
    }

    template <typename T> 
    bool dotu(const hoNDArray<T>& x, const hoNDArray<T>& y, T& r)
    {
        GADGET_DEBUG_CHECK_RETURN_FALSE(x.get_number_of_elements()==y.get_number_of_elements());
        dotu(x.get_number_of_elements(), x.begin(), y.begin(), r);

        return true;
    }

    template EXPORTCPUCOREMATH bool dotu(const hoNDArray<float>& x, const hoNDArray<float>& y, float& r);
    template EXPORTCPUCOREMATH bool dotu(const hoNDArray<double>& x, const hoNDArray<double>& y, double& r);
    template EXPORTCPUCOREMATH bool dotu(const hoNDArray< std::complex<float> >& x, const hoNDArray< std::complex<float> >& y, std::complex<float>& r);
    template EXPORTCPUCOREMATH bool dotu(const hoNDArray< std::complex<double> >& x, const hoNDArray< std::complex<double> >& y, std::complex<double>& r);

    template <typename T> 
    T dotu(const hoNDArray<T>& x, const hoNDArray<T>& y)
    {
        T r = 0;
        dotu(x, y, r);
        return r;
    }

    template EXPORTCPUCOREMATH float dotu(const hoNDArray<float>& x, const hoNDArray<float>& y);
    template EXPORTCPUCOREMATH double dotu(const hoNDArray<double>& x, const hoNDArray<double>& y);
    template EXPORTCPUCOREMATH  std::complex<float>  dotu(const hoNDArray< std::complex<float> >& x, const hoNDArray< std::complex<float> >& y);
    template EXPORTCPUCOREMATH  std::complex<double>  dotu(const hoNDArray< std::complex<double> >& x, const hoNDArray< std::complex<double> >& y);

    // --------------------------------------------------------------------------------

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

    inline void absolute(size_t N, const  std::complex<float> * x, float* r)
    {
        long long n;

        #pragma omp parallel for default(none) private(n) shared(N, x, r) if (N>NumElementsUseThreading)
        for ( n=0; n<(long long)N; n++ )
        {
            const  std::complex<float> & c = x[n];
            const float re = c.real();
            const float im = c.imag();
            r[n]= std::sqrt( (re*re) + (im * im) );
        }
    }

    inline void absolute(size_t N, const  std::complex<double> * x, double* r)
    {
        long long n;

        #pragma omp parallel for default(none) private(n) shared(N, x, r) if (N>NumElementsUseThreading)
        for ( n=0; n<(long long)N; n++ )
        {
            const  std::complex<double> & c = x[n];
            const double re = c.real();
            const double im = c.imag();
            r[n]= std::sqrt( (re*re) + (im * im) );
        }
    }

    template <typename T> 
    bool absolute(const hoNDArray<T>& x, hoNDArray<typename realType<T>::Type>& r)
    {
        if ( r.get_number_of_elements()!=x.get_number_of_elements())
        {
            r.create(x.get_dimensions());
        }

        absolute(x.get_number_of_elements(), x.begin(), r.begin());

        return true;
    }

    template EXPORTCPUCOREMATH bool absolute(const hoNDArray<float>& x, hoNDArray<float>& r);
    template EXPORTCPUCOREMATH bool absolute(const hoNDArray<double>& x, hoNDArray<double>& r);
    template EXPORTCPUCOREMATH bool absolute(const hoNDArray< std::complex<float> >& x, hoNDArray<float>& r);
    template EXPORTCPUCOREMATH bool absolute(const hoNDArray< std::complex<double> >& x, hoNDArray<double>& r);

    // --------------------------------------------------------------------------------

    inline void absolute(size_t N, const std::complex<float>* x, std::complex<float>* r)
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

    inline void absolute(size_t N, const std::complex<double>* x, std::complex<double>* r)
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

    template <typename T> 
    bool absolute(const hoNDArray< std::complex<T> >& x, hoNDArray< std::complex<T> >& r)
    {
        if ( r.get_number_of_elements()!=x.get_number_of_elements())
        {
            r.create(x.get_dimensions());
        }

        absolute(x.get_number_of_elements(), x.begin(), r.begin());

        return true;
    }

    template EXPORTCPUCOREMATH bool absolute(const hoNDArray< std::complex<float> >& x, hoNDArray< std::complex<float> >& r);
    template EXPORTCPUCOREMATH bool absolute(const hoNDArray< std::complex<double> >& x, hoNDArray< std::complex<double> >& r);

    // --------------------------------------------------------------------------------

    template <typename T> 
    bool argument(const hoNDArray<T>& x, hoNDArray<typename realType<T>::Type>& r)
    {
        if ( r.get_number_of_elements()!=x.get_number_of_elements())
        {
            r.create(x.get_dimensions());
        }

        size_t N = x.get_number_of_elements();
        const T* pX = x.begin();
        typename realType<T>::Type* pR = r.begin();

        long long n;

        #pragma omp parallel for default(none) private(n) shared(N, pX, pR) if (N>NumElementsUseThreading)
        for ( n=0; n<(long long)N; n++ )
        {
            pR[n] = std::arg( pX[n] );
        }

        return true;
    }

    template EXPORTCPUCOREMATH bool argument(const hoNDArray< std::complex<float> >& x, hoNDArray<float>& r);
    template EXPORTCPUCOREMATH bool argument(const hoNDArray< std::complex<double> >& x, hoNDArray<double>& r);

    // --------------------------------------------------------------------------------

    template <typename T> 
    bool inv(const hoNDArray<T>& x, hoNDArray<T>& r)
    {
        if ( !r.dimensions_equal(&x) )
        {
            r = x;
        }

        size_t N = x.get_number_of_elements();
        const T* pX = x.begin();
        T* pR = r.begin();

        T v(1.0);
        long long n;

        #pragma omp parallel for default(none) private(n) shared(N, pX, pR, v) if (N>NumElementsUseThreading)
        for ( n=0; n<(long long)N; n++ )
        {
            pR[n] = v/pX[n];
        }

        return true;
    }

    template EXPORTCPUCOREMATH bool inv(const hoNDArray<float>& x, hoNDArray<float>& r);
    template EXPORTCPUCOREMATH bool inv(const hoNDArray<double>& x, hoNDArray<double>& r);
    template EXPORTCPUCOREMATH bool inv(const hoNDArray< std::complex<float> >& x, hoNDArray< std::complex<float> >& r);
    template EXPORTCPUCOREMATH bool inv(const hoNDArray< std::complex<double> >& x, hoNDArray< std::complex<double> >& r);

    // --------------------------------------------------------------------------------

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

    template<typename T> 
    bool conv2(const hoNDArray<T>& x, const hoNDArray<T>& y, hoNDArray<T>& z)
    {
        try
        {
            if ( !z.dimensions_equal(&x) )
            {
                z = x;
            }

            long long RO = (long long) x.get_size(0);
            long long E1 = (long long) x.get_size(1);
            long long num = ((long long) x.get_number_of_elements()) / (RO*E1);

            long long kRO = (long long) y.get_size(0);
            long long kE1 = (long long) y.get_size(1);

            conv2(RO, E1, num, x.begin(), kRO, kE1, y.begin(), z.begin());
        }
        catch(...)
        {
            GADGET_ERROR_MSG("Errors happened in conv2(const hoNDArray<T>& x, const hoNDArray<T>& y, hoNDArray<T>& z) ... ");
            return false;
        }

        return true;
    }

    template EXPORTCPUCOREMATH bool conv2(const hoNDArray<float>& x, const hoNDArray<float>& y, hoNDArray<float>& z);
    template EXPORTCPUCOREMATH bool conv2(const hoNDArray<double>& x, const hoNDArray<double>& y, hoNDArray<double>& z);
    template EXPORTCPUCOREMATH bool conv2(const hoNDArray< std::complex<float> >& x, const hoNDArray< std::complex<float> >& y, hoNDArray< std::complex<float> >& z);
    template EXPORTCPUCOREMATH bool conv2(const hoNDArray< std::complex<double> >& x, const hoNDArray< std::complex<double> >& y, hoNDArray< std::complex<double> >& z);

    // --------------------------------------------------------------------------------

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

    template<typename T> 
    bool conv3(const hoNDArray<T>& x, const hoNDArray<T>& y, hoNDArray<T>& z)
    {
        try
        {
            if ( !z.dimensions_equal(&x) )
            {
                z = x;
            }

            long long RO = (long long) x.get_size(0);
            long long E1 = (long long) x.get_size(1);
            long long E2 = (long long) x.get_size(2);
            long long num = ((long long)x.get_number_of_elements()) / (RO*E1*E2);

            long long kRO = (long long) y.get_size(0);
            long long kE1 = (long long) y.get_size(1);
            long long kE2 = (long long) y.get_size(2);

            conv3(RO, E1, E2, num, x.begin(), kRO, kE1, kE2, y.begin(), z.begin());
        }
        catch(...)
        {
            GADGET_ERROR_MSG("Errors happened in conv3(const hoNDArray<T>& x, const hoNDArray<T>& y, hoNDArray<T>& z) ... ");
            return false;
        }

        return true;
    }

    template EXPORTCPUCOREMATH bool conv3(const hoNDArray<float>& x, const hoNDArray<float>& y, hoNDArray<float>& z);
    template EXPORTCPUCOREMATH bool conv3(const hoNDArray<double>& x, const hoNDArray<double>& y, hoNDArray<double>& z);
    template EXPORTCPUCOREMATH bool conv3(const hoNDArray< std::complex<float> >& x, const hoNDArray< std::complex<float> >& y, hoNDArray< std::complex<float> >& z);
    template EXPORTCPUCOREMATH bool conv3(const hoNDArray< std::complex<double> >& x, const hoNDArray< std::complex<double> >& y, hoNDArray< std::complex<double> >& z);

    // --------------------------------------------------------------------------------

    inline void axpy(float a, size_t N, const float* x, const float* y, float* r)
    {
        long long n;

        #pragma omp parallel for default(none) private(n) shared(N, r, a , x, y) if(N>NumElementsUseThreading)
        for ( n=0; n<(long long)N; ++n)
        {
            r[n] = a*x[n] + y[n];
        }
    }

    inline void axpy(double a, size_t N, const double* x, const double* y, double* r)
    {
        long long n;

        #pragma omp parallel for default(none) private(n) shared(N, r, a , x, y) if(N>NumElementsUseThreading)
        for ( n=0; n<(long long)N; ++n)
        {
            r[n] = a*x[n] + y[n];
        }
    }

    inline void axpy( std::complex<float>  a, size_t N, const  std::complex<float> * x, const  std::complex<float> * y,  std::complex<float> * r)
    {
        long long n;

        #pragma omp parallel for default(none) private(n) shared(N, r, a, x, y) if(N>NumElementsUseThreading)
        for ( n=0; n<(long long)N; ++n)
        {
            const  std::complex<float> & vx = x[n];
            const float re1 = vx.real();
            const float im1 = vx.imag();

            const  std::complex<float> & vy = y[n];
            const float re2 = vy.real();
            const float im2 = vy.imag();

            const float ar = a.real();
            const float ai = a.imag();

            reinterpret_cast<float(&)[2]>(r[n])[0] = re2 + ar*re1 - ai*im1;
            reinterpret_cast<float(&)[2]>(r[n])[1] = im2 + ar*im1 + ai*re1;
        }
    }

    inline void axpy( std::complex<double>  a, size_t N, const  std::complex<double> * x, const  std::complex<double> * y,  std::complex<double> * r)
    {
        long long n;

        #pragma omp parallel for default(none) private(n) shared(N, r, a, x, y) if(N>NumElementsUseThreading)
        for ( n=0; n<(long long)N; ++n)
        {
            const  std::complex<double> & vx = x[n];
            const double re1 = vx.real();
            const double im1 = vx.imag();

            const  std::complex<double> & vy = y[n];
            const double re2 = vy.real();
            const double im2 = vy.imag();

            const double ar = a.real();
            const double ai = a.imag();

            reinterpret_cast<double(&)[2]>(r[n])[0] = re2 + ar*re1 - ai*im1;
            reinterpret_cast<double(&)[2]>(r[n])[1] = im2 + ar*im1 + ai*re1;
        }
    }

    template <typename T> 
    bool axpy(T a, const hoNDArray<T>& x, const hoNDArray<T>& y, hoNDArray<T>& r)
    {
        GADGET_CHECK_RETURN_FALSE(x.get_number_of_elements()==y.get_number_of_elements());

        if ( r.get_number_of_elements() != x.get_number_of_elements() )
        {
            r = y;
        }
        else
        {
            if ( &r != &y )
            {
                memcpy(r.begin(), y.begin(), r.get_number_of_bytes());
            }
        }

        axpy(a, x.get_number_of_elements(), x.begin(), y.begin(), r.begin());

        return true;
    }

    template EXPORTCPUCOREMATH bool axpy(float a, const hoNDArray<float>& x, const hoNDArray<float>& y, hoNDArray<float>& r);
    template EXPORTCPUCOREMATH bool axpy(double a, const hoNDArray<double>& x, const hoNDArray<double>& y, hoNDArray<double>& r);
    template EXPORTCPUCOREMATH bool axpy( std::complex<float>  a, const hoNDArray< std::complex<float> >& x, const hoNDArray< std::complex<float> >& y, hoNDArray< std::complex<float> >& r);
    template EXPORTCPUCOREMATH bool axpy( std::complex<double>  a, const hoNDArray< std::complex<double> >& x, const hoNDArray< std::complex<double> >& y, hoNDArray< std::complex<double> >& r);

    // --------------------------------------------------------------------------------

    inline void scal(size_t N, float a, float* x)
    {
        long long n;
        #pragma omp parallel for default(none) private(n) shared(N, x, a) if (N>NumElementsUseThreading)
        for (n = 0; n < (long long)N; n++)
        {
            x[n] *= a;
        }
    }

    inline void scal(size_t N, double a, double* x)
    {
        long long n;
        #pragma omp parallel for default(none) private(n) shared(N, x, a) if (N>NumElementsUseThreading)
        for (n = 0; n < (long long)N; n++)
        {
            x[n] *= a;
        }
    }

    inline void scal(size_t N,  std::complex<float>  a,  std::complex<float> * x)
    {
        long long n;

        #pragma omp parallel for default(none) private(n) shared(N, x, a) if (N>NumElementsUseThreading)
        for (n = 0; n < (long long)N; n++)
        {
            const  std::complex<float> & c = x[n];
            const float re = c.real();
            const float im = c.imag();

            const float ar = a.real();
            const float ai = a.imag();

            reinterpret_cast<float(&)[2]>(x[n])[0] = re*ar-im*ai;
            reinterpret_cast<float(&)[2]>(x[n])[1] = re*ai+im*ar;
        }
    }

    inline void scal(size_t N,  std::complex<double>  a,  std::complex<double> * x)
    {
        long long n;

        #pragma omp parallel for default(none) private(n) shared(N, x, a) if (N>NumElementsUseThreading)
        for (n = 0; n < (long long)N; n++)
        {
            const  std::complex<double> & c = x[n];
            const double re = c.real();
            const double im = c.imag();

            const double ar = a.real();
            const double ai = a.imag();

            reinterpret_cast<double(&)[2]>(x[n])[0] = re*ar-im*ai;
            reinterpret_cast<double(&)[2]>(x[n])[1] = re*ai+im*ar;
        }
    }

    template <typename T> 
    bool scal(T a, hoNDArray<T>& x)
    {
        scal(x.get_number_of_elements(), a, x.begin());
        return true;
    }

    template EXPORTCPUCOREMATH bool scal(float a, hoNDArray<float>& x);
    template EXPORTCPUCOREMATH bool scal(double a, hoNDArray<double>& x);
    template EXPORTCPUCOREMATH bool scal( std::complex<float>  a, hoNDArray< std::complex<float> >& x);
    template EXPORTCPUCOREMATH bool scal( std::complex<double>  a, hoNDArray< std::complex<double> >& x);

    // --------------------------------------------------------------------------------

    inline void scal(size_t N, float a,  std::complex<float> * x)
    {
        long long n;

        #pragma omp parallel for default(none) private(n) shared(N, x, a) if (N>NumElementsUseThreading)
        for (n = 0; n < (long long)N; n++)
        {
            const  std::complex<float> & c = x[n];
            const float re = c.real();
            const float im = c.imag();

            reinterpret_cast<float(&)[2]>(x[n])[0] = re*a;
            reinterpret_cast<float(&)[2]>(x[n])[1] = im*a;
        }
    }

    inline void scal(size_t N, double a,  std::complex<double> * x)
    {
        long long n;

        #pragma omp parallel for default(none) private(n) shared(N, x, a) if (N>NumElementsUseThreading)
        for (n = 0; n < (long long)N; n++)
        {
            const  std::complex<double> & c = x[n];
            const double re = c.real();
            const double im = c.imag();

            reinterpret_cast<double(&)[2]>(x[n])[0] = re*a;
            reinterpret_cast<double(&)[2]>(x[n])[1] = im*a;
        }
    }

    template <typename T> 
    bool scal(T a, hoNDArray< std::complex<T> >& x)
    {
        scal(x.get_number_of_elements(), a, x.begin());
        return true;
    }

    template EXPORTCPUCOREMATH bool scal(float a, hoNDArray< std::complex<float> >& x);
    template EXPORTCPUCOREMATH bool scal(double a, hoNDArray< std::complex<double> >& x);

    // --------------------------------------------------------------------------------

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

    template <typename T> 
    bool sort(const hoNDArray<T>& x, hoNDArray<T>& r, bool isascending)
    {
        if ( &r != &x )
        {
            if ( r.get_number_of_elements()!=x.get_number_of_elements())
            {
                r = x;
            }
            else
            {
                memcpy(r.begin(), x.begin(), x.get_number_of_bytes());
            }
        }

        sort(x.get_number_of_elements(), x.begin(), r.begin(), isascending);

        return true;
    }

    template EXPORTCPUCOREMATH bool sort(const hoNDArray<float>& x, hoNDArray<float>& r, bool isascending);
    template EXPORTCPUCOREMATH bool sort(const hoNDArray<double>& x, hoNDArray<double>& r, bool isascending);

// --------------------------------------------------------------------------------

    template<typename T> void fill( hoNDArray<T>* x, T val)
    {
        size_t N = x->get_number_of_elements();
        T* pX = x->begin();

        long long n;
        #pragma omp parallel for default(none) private(n) shared(N, pX, val) if (N>NumElementsUseThreading)
        for ( n=0; n<(long long)N; n++ )
        {
            pX[n] = val;
        }
    }

    template EXPORTCPUCOREMATH void fill( hoNDArray<float>* x, float val);
    template EXPORTCPUCOREMATH void fill( hoNDArray<double>* x, double val);
    template EXPORTCPUCOREMATH void fill( hoNDArray< std::complex<float> >* x,  std::complex<float>  val);
    template EXPORTCPUCOREMATH void fill( hoNDArray< std::complex<double> >* x,  std::complex<double>  val);

    // --------------------------------------------------------------------------------

    template<typename T> void fill( hoNDArray<T>& x, T val )
    {
        Gadgetron::fill( &x, val);
    }

    template EXPORTCPUCOREMATH void fill( hoNDArray<float>& x, float val);
    template EXPORTCPUCOREMATH void fill( hoNDArray<double>& x, double val);
    template EXPORTCPUCOREMATH void fill( hoNDArray< std::complex<float> >& x,  std::complex<float>  val);
    template EXPORTCPUCOREMATH void fill( hoNDArray< std::complex<double> >& x,  std::complex<double>  val);

    // --------------------------------------------------------------------------------

    inline void asum(size_t N, const float* x, float& r)
    {
        long long i;
        float sum(0);
        #pragma omp parallel for private(i) reduction(+:sum) if (N>NumElementsUseThreading)
        for (i = 0; i < (long long)N; i++)
        {
            sum += GT_ABS(x[i]);
        }

        r = sum;
    }

    inline void asum(size_t N, const double* x, double& r)
    {
        long long i;
        double sum(0);
        #pragma omp parallel for private(i) reduction(+:sum) if (N>NumElementsUseThreading)
        for (i = 0; i < (long long)N; i++)
        {
            sum += GT_ABS(x[i]);
        }

        r = sum;
    }

    inline void asum(size_t N, const  std::complex<float> * x, float& r)
    {
        long long i;
        float sum(0);
        #pragma omp parallel for private(i) reduction(+:sum) if (N>NumElementsUseThreading)
        for (i = 0; i < (long long)N; i++)
        {
            const  std::complex<float> & c = x[i];
            const float re = c.real();
            const float im = c.imag();
            sum += ( GT_ABS(re) + GT_ABS(im) );
        }

        r = sum;
    }

    inline void asum(size_t N, const  std::complex<double> * x, double& r)
    {
        long long i;
        double sum(0);
        #pragma omp parallel for private(i) reduction(+:sum) if (N>NumElementsUseThreading)
        for (i = 0; i < (long long)N; i++)
        {
            const  std::complex<double> & c = x[i];
            const double re = c.real();
            const double im = c.imag();
            sum += ( GT_ABS(re) + GT_ABS(im) );
        }

        r = sum;
    }

    template<class T> void asum(const hoNDArray<T>& x, typename realType<T>::Type& r)
    {
        asum(x.get_number_of_elements(), x.begin(), r);
    }

    template EXPORTCPUCOREMATH void asum( const hoNDArray<float>& x, float& r);
    template EXPORTCPUCOREMATH void asum( const hoNDArray<double>& x, double& r);
    template EXPORTCPUCOREMATH void asum( const hoNDArray< std::complex<float> >& x, float& r);
    template EXPORTCPUCOREMATH void asum( const hoNDArray< std::complex<double> >& x, double& r);

    template<class T> typename realType<T>::Type asum(const hoNDArray<T>& x)
    {
        typename realType<T>::Type r;
        asum(x, r);
        return r;
    }

    template EXPORTCPUCOREMATH float asum( const hoNDArray<float>& x);
    template EXPORTCPUCOREMATH double asum( const hoNDArray<double>& x);
    template EXPORTCPUCOREMATH float asum( const hoNDArray< std::complex<float> >& x);
    template EXPORTCPUCOREMATH double asum( const hoNDArray< std::complex<double> >& x);

    // --------------------------------------------------------------------------------

    inline size_t amax(size_t N, const float* x)
    {
        lapack_int num = (lapack_int)(N);
        lapack_int incx = 1;

        return isamax_(&num, (float*)(x), &incx);
    }

    inline size_t amax(size_t N, const double* x)
    {
        lapack_int num = (lapack_int)(N);
        lapack_int incx = 1;

        return idamax_(&num, (double*)(x), &incx);
    }

    inline size_t amax(size_t N, const  std::complex<float> * x)
    {
        lapack_int num = (lapack_int)(N);
        lapack_int incx = 1;

        return icamax_(&num, (lapack_complex_float*)(x), &incx);
    }

    inline size_t amax(size_t N, const  std::complex<double> * x)
    {
        lapack_int num = (lapack_int)(N);
        lapack_int incx = 1;

        return izamax_(&num, (lapack_complex_double*)(x), &incx);
    }

    template<class T> size_t amax(const hoNDArray<T>& x)
    {
        return amax(x.get_number_of_elements(), x.begin());
    }

    template EXPORTCPUCOREMATH size_t amax( const hoNDArray<float>& x);
    template EXPORTCPUCOREMATH size_t amax( const hoNDArray<double>& x);
    template EXPORTCPUCOREMATH size_t amax( const hoNDArray< std::complex<float> >& x);
    template EXPORTCPUCOREMATH size_t amax( const hoNDArray< std::complex<double> >& x);

    // --------------------------------------------------------------------------------

    template<class T> 
    bool real_imag_to_complex(const hoNDArray<typename realType<T>::Type>& real, const hoNDArray<typename realType<T>::Type>& imag, hoNDArray<T>& cplx)
    {
        try
        {
            GADGET_CHECK_RETURN_FALSE(real.dimensions_equal(&imag));

            if ( !cplx.dimensions_equal(&real) )
            {
                cplx.create(real.get_dimensions());
            }

            T* pRes = cplx.begin();
            const typename realType<T>::Type* pReal = real.begin();
            const typename realType<T>::Type* pImag = imag.begin();

            size_t N = real.get_number_of_elements();

            long long n;
            #pragma omp parallel for private(n) shared(N, pRes, pReal, pImag)
            for ( n=0; n<(long long)N; n++ )
            {
                pRes[n] = T(pReal[n], pImag[n]);
            }
        }
        catch(...)
        {
            GADGET_ERROR_MSG("Errors in real_imag_to_complex(...) ... ");
            return false;
        }

        return true;
    }

    template EXPORTCPUCOREMATH bool real_imag_to_complex(const hoNDArray<float>& real, const hoNDArray<float>& imag, hoNDArray< std::complex<float> >& cplx);
    template EXPORTCPUCOREMATH bool real_imag_to_complex(const hoNDArray<double>& real, const hoNDArray<double>& imag, hoNDArray< std::complex<double> >& cplx);

    // --------------------------------------------------------------------------------

    template<class T> 
    bool complex_to_real_imag(const hoNDArray<T>& cplx, hoNDArray<typename realType<T>::Type>& real, hoNDArray<typename realType<T>::Type>& imag)
    {
        try
        {
            if ( !real.dimensions_equal(&cplx) )
            {
                real.create(cplx.get_dimensions());
            }

            if ( !imag.dimensions_equal(&cplx) )
            {
                imag.create(cplx.get_dimensions());
            }

            const T* pRes = cplx.begin();
            typename realType<T>::Type* pReal = real.begin();
            typename realType<T>::Type* pImag = imag.begin();

            size_t N = real.get_number_of_elements();

            long long n;
            #pragma omp parallel for default(none) private(n) shared(N, pRes, pReal, pImag)
            for ( n=0; n<(long long)N; n++ )
            {
                pReal[n] = pRes[n].real();
                pImag[n] = pRes[n].imag();
            }
        }
        catch(...)
        {
            GADGET_ERROR_MSG("Errors in complex_to_real_imag(...) ... ");
            return false;
        }

        return true;
    }

    template EXPORTCPUCOREMATH bool complex_to_real_imag(const hoNDArray< std::complex<float> >& cplx, hoNDArray<float>& real, hoNDArray<float>& imag);
    template EXPORTCPUCOREMATH bool complex_to_real_imag(const hoNDArray< std::complex<double> >& cplx, hoNDArray<double>& real, hoNDArray<double>& imag);

    template <> EXPORTCPUCOREMATH
    bool complex_to_real_imag(const hoNDArray<float>& cplx, hoNDArray<float>& real, hoNDArray<float>& imag)
    {
        try
        {
            if ( !real.dimensions_equal(&cplx) )
            {
                real.create(cplx.get_dimensions());
            }

            if ( !imag.dimensions_equal(&cplx) )
            {
                imag.create(cplx.get_dimensions());
            }

            const float* pRes = cplx.begin();
            float* pReal = real.begin();
            float* pImag = imag.begin();

            size_t N = real.get_number_of_elements();

            long long n;
            #pragma omp parallel for default(none) private(n) shared(N, pRes, pReal, pImag)
            for ( n=0; n<(long long)N; n++ )
            {
                pReal[n] = pRes[n];
                pImag[n] = 0;
            }
        }
        catch(...)
        {
            GADGET_ERROR_MSG("Errors in complex_to_real_imag(...) ... ");
            return false;
        }

        return true;
    }

    template<> EXPORTCPUCOREMATH 
    bool complex_to_real_imag(const hoNDArray<double>& cplx, hoNDArray<double>& real, hoNDArray<double>& imag)
    {
        try
        {
            if ( !real.dimensions_equal(&cplx) )
            {
                real.create(cplx.get_dimensions());
            }

            if ( !imag.dimensions_equal(&cplx) )
            {
                imag.create(cplx.get_dimensions());
            }

            const double* pRes = cplx.begin();
            double* pReal = real.begin();
            double* pImag = imag.begin();

            size_t N = real.get_number_of_elements();

            long long n;
            #pragma omp parallel for default(none) private(n) shared(N, pRes, pReal, pImag)
            for ( n=0; n<(long long)N; n++ )
            {
                pReal[n] = pRes[n];
                pImag[n] = 0;
            }
        }
        catch(...)
        {
            GADGET_ERROR_MSG("Errors in complex_to_real_imag(...) ... ");
            return false;
        }

        return true;
    }

    // --------------------------------------------------------------------------------

    template<class T> 
    bool complex_to_real(const hoNDArray<T>& cplx, hoNDArray<typename realType<T>::Type>& real)
    {
        try
        {
            if ( !real.dimensions_equal(&cplx) )
            {
                real.create(cplx.get_dimensions());
            }

            const T* pRes = cplx.begin();
            typename realType<T>::Type* pReal = real.begin();

            size_t N = real.get_number_of_elements();

            long long n;
            #pragma omp parallel for default(none) private(n) shared(N, pRes, pReal)
            for ( n=0; n<(long long)N; n++ )
            {
                pReal[n] = pRes[n].real();
            }
        }
        catch(...)
        {
            GADGET_ERROR_MSG("Errors in complex_to_real(...) ... ");
            return false;
        }

        return true;
    }

    template EXPORTCPUCOREMATH bool complex_to_real(const hoNDArray< std::complex<float> >& cplx, hoNDArray<float>& real);
    template EXPORTCPUCOREMATH bool complex_to_real(const hoNDArray< std::complex<double> >& cplx, hoNDArray<double>& real);

    template<class T> 
    bool complex_to_real(const hoNDArray<T>& cplx, hoNDArray<T>& real)
    {
        try
        {
            if ( !real.dimensions_equal(&cplx) )
            {
                real.create(cplx.get_dimensions());
            }

            const T* pRes = cplx.begin();
            T* pReal = real.begin();

            size_t N = real.get_number_of_elements();

            long long n;
            #pragma omp parallel for private(n) shared(N, pRes, pReal)
            for ( n=0; n<(long long)N; n++ )
            {
                pReal[n] = T(pRes[n].real(), 0);
            }
        }
        catch(...)
        {
            GADGET_ERROR_MSG("Errors in complex_to_real(...) ... ");
            return false;
        }

        return true;
    }

    template EXPORTCPUCOREMATH bool complex_to_real(const hoNDArray< std::complex<float> >& cplx, hoNDArray< std::complex<float> >& real);
    template EXPORTCPUCOREMATH bool complex_to_real(const hoNDArray< std::complex<double> >& cplx, hoNDArray< std::complex<double> >& real);

    template<class T> 
    bool complex_to_real(hoNDArray<T>& cplx)
    {
        try
        {
            T* pRes = cplx.begin();

            size_t N = cplx.get_number_of_elements();

            long long n;
            #pragma omp parallel for private(n) shared(N, pRes)
            for ( n=0; n<(long long)N; n++ )
            {
                pRes[n] = T(pRes[n].real(), 0);
            }
        }
        catch(...)
        {
            GADGET_ERROR_MSG("Errors in complex_to_real(...) ... ");
            return false;
        }

        return true;
    }

    template EXPORTCPUCOREMATH bool complex_to_real(hoNDArray< std::complex<float> >& cplx);
    template EXPORTCPUCOREMATH bool complex_to_real(hoNDArray< std::complex<double> >& cplx);

    // --------------------------------------------------------------------------------

    template<class T> 
    bool complex_to_imag(const hoNDArray<T>& cplx, hoNDArray<typename realType<T>::Type>& imag)
    {
        try
        {
            if ( !imag.dimensions_equal(&cplx) )
            {
                imag.create(cplx.get_dimensions());
            }

            const T* pRes = cplx.begin();
            typename realType<T>::Type* pImag = imag.begin();

            size_t N = imag.get_number_of_elements();

            long long n;
            #pragma omp parallel for default(none) private(n) shared(N, pRes, pImag)
            for ( n=0; n<(long long)N; n++ )
            {
                pImag[n] = pRes[n].imag();
            }
        }
        catch(...)
        {
            GADGET_ERROR_MSG("Errors in complex_to_imag(...) ... ");
            return false;
        }

        return true;
    }

    template EXPORTCPUCOREMATH bool complex_to_imag(const hoNDArray< std::complex<float> >& cplx, hoNDArray<float>& imag);
    template EXPORTCPUCOREMATH bool complex_to_imag(const hoNDArray< std::complex<double> >& cplx, hoNDArray<double>& imag);

    template<class T> 
    bool complex_to_imag(const hoNDArray<T>& cplx, hoNDArray<T>& imag)
    {
        try
        {
            if ( !imag.dimensions_equal(&cplx) )
            {
                imag.create(cplx.get_dimensions());
            }

            const T* pRes = cplx.begin();
            T* pImag = imag.begin();

            size_t N = imag.get_number_of_elements();

            long long n;
            #pragma omp parallel for private(n) shared(N, pRes, pImag)
            for ( n=0; n<(long long)N; n++ )
            {
                pImag[n] = T(0, pRes[n].imag());
            }
        }
        catch(...)
        {
            GADGET_ERROR_MSG("Errors in complex_to_imag(...) ... ");
            return false;
        }

        return true;
    }

    template EXPORTCPUCOREMATH bool complex_to_imag(const hoNDArray< std::complex<float> >& cplx, hoNDArray< std::complex<float> >& imag);
    template EXPORTCPUCOREMATH bool complex_to_imag(const hoNDArray< std::complex<double> >& cplx, hoNDArray< std::complex<double> >& imag);

    template<class T> 
    bool complex_to_imag(hoNDArray<T>& cplx)
    {
        try
        {
            T* pRes = cplx.begin();

            size_t N = cplx.get_number_of_elements();

            long long n;
            #pragma omp parallel for private(n) shared(N, pRes)
            for ( n=0; n<(long long)N; n++ )
            {
                pRes[n] = T( pRes[n].real(), 0);
            }
        }
        catch(...)
        {
            GADGET_ERROR_MSG("Errors in complex_to_imag(...) ... ");
            return false;
        }

        return true;
    }

    template EXPORTCPUCOREMATH bool complex_to_imag(hoNDArray< std::complex<float> >& cplx);
    template EXPORTCPUCOREMATH bool complex_to_imag(hoNDArray< std::complex<double> >& cplx);

    // --------------------------------------------------------------------------------

    template<class T> 
    bool real_to_complex(const hoNDArray<typename realType<T>::Type>& real, hoNDArray<T>& cplx)
    {
        try
        {
            if ( !cplx.dimensions_equal(&real) )
            {
                cplx.create(real.get_dimensions());
            }

            const typename realType<T>::Type* pReal = real.begin();
            T* pRes = cplx.begin();

            size_t N = real.get_number_of_elements();

            long long n;
            #pragma omp parallel for private(n) shared(N, pRes, pReal)
            for ( n=0; n<(long long)N; n++ )
            {
                pRes[n] = T(pReal[n], 0);
            }
        }
        catch(...)
        {
            GADGET_ERROR_MSG("Errors in real_to_complex(...) ... ");
            return false;
        }

        return true;
    }

    template EXPORTCPUCOREMATH bool real_to_complex(const hoNDArray< float >& real, hoNDArray< std::complex<float> >& cplx);
    template EXPORTCPUCOREMATH bool real_to_complex(const hoNDArray< double >& real, hoNDArray< std::complex<double> >& cplx);

    // --------------------------------------------------------------------------------

    template <class T>
    bool minValue(const hoNDArray<T>& a, T& v)
    {
        typedef T ValueType;

        try
        {
            const ValueType* pA = a.begin();
            size_t n = a.get_number_of_elements();
            v = pA[0];

            size_t ii;
            for (ii=1; ii<n; ii++)
            {
                if (pA[ii]<v) v = pA[ii];
            }
        }
        catch(...)
        {
            GADGET_ERROR_MSG("Errors in minValue(const hoNDArray<T>& a, T& v) ... ");
            return false;
        }

        return true;
    }

    template EXPORTCPUCOREMATH bool minValue(const hoNDArray<float>& a, float& v);
    template EXPORTCPUCOREMATH bool minValue(const hoNDArray<double>& a, double& v);

    template <class T>
    bool maxValue(const hoNDArray<T>& a, T& v)
    {
        typedef T ValueType;

        try
        {
            const ValueType* pA = a.begin();
            size_t n = a.get_number_of_elements();
            v = pA[0];

            size_t ii;
            for (ii=1; ii<n; ii++)
            {
                if (pA[ii]>v) v = pA[ii];
            }
        }
        catch(...)
        {
            GADGET_ERROR_MSG("Errors in maxValue(const hoNDArray<T>& a, T& v) ... ");
            return false;
        }

        return true;
    }

    template EXPORTCPUCOREMATH bool maxValue(const hoNDArray<float>& a, float& v);
    template EXPORTCPUCOREMATH bool maxValue(const hoNDArray<double>& a, double& v);

    // --------------------------------------------------------------------------------
}
