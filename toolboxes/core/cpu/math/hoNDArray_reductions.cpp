#include "hoNDArray_reductions.h"
#include "hoArmadillo.h"

#ifndef lapack_int
    #define lapack_int int
#endif // lapack_int

#ifndef lapack_complex_float
    #define lapack_complex_float  std::complex<float> 
#endif // lapack_complex_float

#ifndef lapack_complex_double
    #define lapack_complex_double  std::complex<double> 
#endif // #ifndef lapack_complex_double

#define NumElementsUseThreading 64*1024

//Declaration of BLAS and LAPACK routines
extern "C"
{
    /// Finds the index of the element with the maximal absolute value.
    lapack_int isamax_(lapack_int* N, float* x, lapack_int* incx);
    lapack_int idamax_(lapack_int* N, double* x, lapack_int* incx);
    lapack_int icamax_(lapack_int* N, lapack_complex_float* x, lapack_int* incx);
    lapack_int izamax_(lapack_int* N, lapack_complex_double* x, lapack_int* incx);
}

namespace Gadgetron{

    // --------------------------------------------------------------------------------

    template<class REAL> REAL max(hoNDArray<REAL>* data){
        return as_arma_col(data).max();
    }

    // --------------------------------------------------------------------------------

    template<class REAL> REAL min(hoNDArray<REAL>* data){
        return as_arma_col(data).min();
    }

    // --------------------------------------------------------------------------------

    template<class T> T mean(hoNDArray<T>* data){
        return (typename stdType<T>::Type) arma::mean(as_arma_col(data));
    }

    // --------------------------------------------------------------------------------

    template<class T> T sum(hoNDArray<T>* data){
        return (typename stdType<T>::Type) arma::sum(as_arma_col(data));
    }

    // --------------------------------------------------------------------------------

    template<class T> T stddev(hoNDArray<T>* data){
        return (typename stdType<T>::Type) arma::stddev(as_arma_col(data));
    }

    // --------------------------------------------------------------------------------

    template<class T> T var(hoNDArray<T>* data) {
        return (typename stdType<T>::Type) arma::var(as_arma_col(data));
    }

    // --------------------------------------------------------------------------------

     template<class T> T median(hoNDArray<T>* data){
        return (typename stdType<T>::Type) arma::median(as_arma_col(data));
    }

    // --------------------------------------------------------------------------------

    template<class T> T dot( hoNDArray<T> *x, hoNDArray<T> *y, bool cc )
    {
        if( x == 0x0 || y == 0x0 )
            throw std::runtime_error("Gadgetron::dot(): Invalid input array");

        if( x->get_number_of_elements() != y->get_number_of_elements() )
            throw std::runtime_error("Gadgetron::dot(): Array sizes mismatch");

        arma::Col<typename stdType<T>::Type> xM = as_arma_col(x);
        arma::Col<typename stdType<T>::Type> yM = as_arma_col(y);
        typename stdType<T>::Type res = (cc) ? arma::cdot(xM,yM) : arma::dot(xM,yM);
        return *((T*)(&res));
    }

    // --------------------------------------------------------------------------------

    inline void asum(size_t N, const float* x, float& r)
    {
        long long i;
        float sum(0);
        #pragma omp parallel for private(i) reduction(+:sum) if (N>NumElementsUseThreading)
        for (i = 0; i < (long long)N; i++)
        {
            sum += std::abs(x[i]);
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
            sum += std::abs(x[i]);
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
            sum += ( std::abs(re) + std::abs(im) );
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
            sum += ( std::abs(re) + std::abs(im) );
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

    template<class T> typename realType<T>::Type asum( hoNDArray<T> *x )
    {
        if( x == 0x0 )
            throw std::runtime_error("Gadgetron::asum(): Invalid input array");

        typedef typename realType<T>::Type realT;
        arma::Col<typename stdType<T>::Type> xM = as_arma_col(x);
        return realT(arma::norm(xM,1));
    }

    template<class T> T asum( hoNDArray< std::complex<T> > *x )
    {
        if( x == 0x0 )
            throw std::runtime_error("Gadgetron::asum(): Invalid input array");

        return arma::norm(arma::abs(real(as_arma_col(x)))+arma::abs(imag(as_arma_col(x))),1);
    }

    template<class T> T asum( hoNDArray< complext<T> > *x )
    {
        if( x == 0x0 )
            throw std::runtime_error("Gadgetron::asum(): Invalid input array");

        return arma::norm(arma::abs(real(as_arma_col(x)))+arma::abs(imag(as_arma_col(x))),1);
    }

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
            norm1Sum += std::abs(c);
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

    inline void norm1(size_t N, const  complext<float> * x, float& r)
    {
        norm1(N, (std::complex<float> *)x, r);
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

    inline void norm1(size_t N, const  complext<double> * x, double& r)
    {
        norm1(N, (std::complex<double> *)x, r);
    }

    template <typename T> 
    void norm1(const hoNDArray<T>& x, typename realType<T>::Type& r)
    {
        norm1(x.get_number_of_elements(), x.begin(), r);
    }

    template EXPORTCPUCOREMATH void norm1(const hoNDArray<float>& x, float& r);
    template EXPORTCPUCOREMATH void norm1(const hoNDArray<double>& x, double& r);
    template EXPORTCPUCOREMATH void norm1(const hoNDArray< std::complex<float> >& x, float& r);
    template EXPORTCPUCOREMATH void norm1(const hoNDArray< complext<float> >& x, float& r);
    template EXPORTCPUCOREMATH void norm1(const hoNDArray< std::complex<double> >& x, double& r);
    template EXPORTCPUCOREMATH void norm1(const hoNDArray< complext<double> >& x, double& r);

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
    template EXPORTCPUCOREMATH float norm1(const hoNDArray< complext<float> >& x);
    template EXPORTCPUCOREMATH double norm1(const hoNDArray< std::complex<double> >& x);
    template EXPORTCPUCOREMATH double norm1(const hoNDArray< complext<double> >& x);

    template<class T> typename realType<T>::Type nrm1( hoNDArray<T> *x )
    {
        if( x == 0x0 )
            throw std::runtime_error("Gadgetron::nrm2(): Invalid input array");

        /*typedef typename realType<T>::Type realT;
        arma::Col<typename stdType<T>::Type> xM = as_arma_col(x);
        return realT(arma::norm(xM,1));*/

        return norm1(*x);
    }

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

    inline void norm2(size_t N, const  complext<float> * x, float& r)
    {
        norm2(N, (std::complex<float> *)x, r);
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

    inline void norm2(size_t N, const  complext<double> * x, double& r)
    {
        norm2(N, (std::complex<double> *)x, r);
    }

    template <typename T> 
    void norm2(const hoNDArray<T>& x, typename realType<T>::Type& r)
    {
        norm2(x.get_number_of_elements(), x.begin(), r);
    }

    template EXPORTCPUCOREMATH void norm2(const hoNDArray<float>& x, float& r);
    template EXPORTCPUCOREMATH void norm2(const hoNDArray<double>& x, double& r);
    template EXPORTCPUCOREMATH void norm2(const hoNDArray< std::complex<float> >& x, float& r);
    template EXPORTCPUCOREMATH void norm2(const hoNDArray< complext<float> >& x, float& r);
    template EXPORTCPUCOREMATH void norm2(const hoNDArray< std::complex<double> >& x, double& r);
    template EXPORTCPUCOREMATH void norm2(const hoNDArray< complext<double> >& x, double& r);

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
    template EXPORTCPUCOREMATH float norm2(const hoNDArray< complext<float> >& x);
    template EXPORTCPUCOREMATH double norm2(const hoNDArray< std::complex<double> >& x);
    template EXPORTCPUCOREMATH double norm2(const hoNDArray< complext<double> >& x);

    template<class T> typename realType<T>::Type nrm2( hoNDArray<T> *x )
    {
        if( x == 0x0 )
            throw std::runtime_error("Gadgetron::nrm2(): Invalid input array");

        /*typedef typename realType<T>::Type realT;
        arma::Col<typename stdType<T>::Type> xM = as_arma_col(x);
        return realT(arma::norm(xM,2));*/

        return norm2(*x);
    }

    // --------------------------------------------------------------------------------

    template <typename T> 
    void minAbsolute(const hoNDArray<T>& x, T& r, size_t& ind)
    {
        size_t N = x.get_number_of_elements();
        const T* pX = x.begin();

        ind = 0;
        if ( N == 0 ) return;

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
    }

    template EXPORTCPUCOREMATH void minAbsolute(const hoNDArray<float>& x, float& r, size_t& ind);
    template EXPORTCPUCOREMATH void minAbsolute(const hoNDArray<double>& x, double& r, size_t& ind);
    template EXPORTCPUCOREMATH void minAbsolute(const hoNDArray< std::complex<float> >& x,  std::complex<float> & r, size_t& ind);
    template EXPORTCPUCOREMATH void minAbsolute(const hoNDArray< std::complex<double> >& x,  std::complex<double> & r, size_t& ind);

    template<class T> size_t amin( hoNDArray<T> *x )
    {
        if( x == 0x0 )
            throw std::runtime_error("Gadgetron::amin(): Invalid input array");

        typedef typename realType<T>::Type realT;
        arma::Col<realT> xM = arma::abs(as_arma_col(x));
        arma::uword idx;
        realT min = xM.min(idx);
        return idx;
    }

    template<class T> size_t amin( hoNDArray< std::complex<T> > *x )
    {
        if( x == 0x0 )
            throw std::runtime_error("Gadgetron::amin(): Invalid input array");

        arma::Col<T> xM = arma::abs(real(as_arma_col(x)))+arma::abs(imag(as_arma_col(x)));
        arma::uword idx;
        T min = xM.min(idx);
        return idx;
    }

    template<class T> size_t amin( hoNDArray< complext<T> > *x )
    {
        if( x == 0x0 )
            throw std::runtime_error("Gadgetron::amin(): Invalid input array");

        arma::Col<T> xM = arma::abs(real(as_arma_col(x)))+arma::abs(imag(as_arma_col(x)));
        arma::uword idx;
        T min = xM.min(idx);
        return idx;
    }

    // --------------------------------------------------------------------------------

    template <typename T> 
    void maxAbsolute(const hoNDArray<T>& x, T& r, size_t& ind)
    {
        size_t N = x.get_number_of_elements();
        const T* pX = x.begin();

        ind = 0;
        if ( N == 0 ) return;

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
    }

    template EXPORTCPUCOREMATH void maxAbsolute(const hoNDArray<float>& x, float& r, size_t& ind);
    template EXPORTCPUCOREMATH void maxAbsolute(const hoNDArray<double>& x, double& r, size_t& ind);
    template EXPORTCPUCOREMATH void maxAbsolute(const hoNDArray< std::complex<float> >& x,  std::complex<float> & r, size_t& ind);
    template EXPORTCPUCOREMATH void maxAbsolute(const hoNDArray< std::complex<double> >& x,  std::complex<double> & r, size_t& ind);

    // --------------------------------------------------------------------------------

    inline size_t amax(size_t N, const float* x)
    {
        lapack_int num = (lapack_int)(N);
        lapack_int incx = 1;

        return isamax_(&num, (float*)(x), &incx) - size_t(1);
    }

    inline size_t amax(size_t N, const double* x)
    {
        lapack_int num = (lapack_int)(N);
        lapack_int incx = 1;

        return idamax_(&num, (double*)(x), &incx) - size_t(1);
    }

    inline size_t amax(size_t N, const  std::complex<float> * x)
    {
        lapack_int num = (lapack_int)(N);
        lapack_int incx = 1;

        return icamax_(&num, (lapack_complex_float*)(x), &incx) - size_t(1);
    }

    inline size_t amax(size_t N, const  std::complex<double> * x)
    {
        lapack_int num = (lapack_int)(N);
        lapack_int incx = 1;

        return izamax_(&num, (lapack_complex_double*)(x), &incx) - size_t(1);
    }

    template<class T> size_t amax(const hoNDArray<T>& x)
    {
        return amax(x.get_number_of_elements(), x.begin());
    }

    template EXPORTCPUCOREMATH size_t amax( const hoNDArray<float>& x);
    template EXPORTCPUCOREMATH size_t amax( const hoNDArray<double>& x);
    template EXPORTCPUCOREMATH size_t amax( const hoNDArray< std::complex<float> >& x);
    template EXPORTCPUCOREMATH size_t amax( const hoNDArray< std::complex<double> >& x);

    template<class T> size_t amax( hoNDArray<T> *x )
    {
        if( x == 0x0 )
            throw std::runtime_error("Gadgetron::amax(): Invalid input array");

        typedef typename realType<T>::Type realT;
        arma::Col<realT> xM = arma::abs(as_arma_col(x));
        arma::uword idx;
        realT max = xM.max(idx);
        return idx;
    }

    template<class T> size_t amax( hoNDArray< std::complex<T> > *x )
    {
        if( x == 0x0 )
            throw std::runtime_error("Gadgetron::amax(): Invalid input array");

        arma::Col<T> xM = arma::abs(real(as_arma_col(x)))+arma::abs(imag(as_arma_col(x)));
        arma::uword idx;
        T max = xM.max(idx);
        return idx;
    }

    template<class T> size_t amax( hoNDArray< complext<T> > *x )
    {
        if( x == 0x0 )
            throw std::runtime_error("Gadgetron::amax(): Invalid input array");

        arma::Col<T> xM = arma::abs(real(as_arma_col(x)))+arma::abs(imag(as_arma_col(x)));
        arma::uword idx;
        T max = xM.max(idx);
        return idx;
    }

    // --------------------------------------------------------------------------------

    template EXPORTCPUCOREMATH float max(hoNDArray<float>*);
    template EXPORTCPUCOREMATH float min(hoNDArray<float>*);
    template EXPORTCPUCOREMATH float mean(hoNDArray<float>*);
    template EXPORTCPUCOREMATH float median(hoNDArray<float>*);
    template EXPORTCPUCOREMATH float sum(hoNDArray<float>*);
    template EXPORTCPUCOREMATH float stddev(hoNDArray<float>*);
    template EXPORTCPUCOREMATH float var(hoNDArray<float>*);

    template EXPORTCPUCOREMATH double max(hoNDArray<double>*);
    template EXPORTCPUCOREMATH double min(hoNDArray<double>*);
    template EXPORTCPUCOREMATH double mean(hoNDArray<double>*);
    template EXPORTCPUCOREMATH double median(hoNDArray<double>*);
    template EXPORTCPUCOREMATH double sum(hoNDArray<double>*);
    template EXPORTCPUCOREMATH double stddev(hoNDArray<double>*);
    template EXPORTCPUCOREMATH double var(hoNDArray<double>*);

    template EXPORTCPUCOREMATH complext<double> mean(hoNDArray<complext<double> >*);
    template EXPORTCPUCOREMATH complext<double> median(hoNDArray<complext<double> >*);
    template EXPORTCPUCOREMATH complext<double> sum(hoNDArray<complext<double> >*);
    template EXPORTCPUCOREMATH complext<double> stddev(hoNDArray<complext<double> >*);
    template EXPORTCPUCOREMATH complext<double> var(hoNDArray<complext<double> >*);

    template EXPORTCPUCOREMATH complext<float> mean(hoNDArray<complext<float> >*);
    template EXPORTCPUCOREMATH complext<float> median(hoNDArray<complext<float> >*);
    template EXPORTCPUCOREMATH complext<float> sum(hoNDArray<complext<float> >*);
    template EXPORTCPUCOREMATH complext<float> stddev(hoNDArray<complext<float> >*);
    template EXPORTCPUCOREMATH complext<float> var(hoNDArray<complext<float> >*);

    template EXPORTCPUCOREMATH std::complex<double> mean(hoNDArray<std::complex<double> >*);
    template EXPORTCPUCOREMATH std::complex<double> sum(hoNDArray<std::complex<double> >*);
    template EXPORTCPUCOREMATH std::complex<double> stddev(hoNDArray<std::complex<double> >*);
    template EXPORTCPUCOREMATH std::complex<double> var(hoNDArray<std::complex<double> >*);

    template EXPORTCPUCOREMATH std::complex<float> mean(hoNDArray<std::complex<float> >*);
    template EXPORTCPUCOREMATH std::complex<float> sum(hoNDArray<std::complex<float> >*);
    template EXPORTCPUCOREMATH std::complex<float> stddev(hoNDArray<std::complex<float> >*);
    template EXPORTCPUCOREMATH std::complex<float> var(hoNDArray<std::complex<float> >*);

    template EXPORTCPUCOREMATH float dot<float>( hoNDArray<float>*, hoNDArray<float>*, bool );
    template EXPORTCPUCOREMATH float asum<float>( hoNDArray<float>* );
    template EXPORTCPUCOREMATH float nrm2<float>( hoNDArray<float>* );

    template EXPORTCPUCOREMATH size_t amin<float>( hoNDArray<float>* );
    template EXPORTCPUCOREMATH size_t amax<float>( hoNDArray<float>* );

    template EXPORTCPUCOREMATH double dot<double>( hoNDArray<double>*, hoNDArray<double>*, bool );
    template EXPORTCPUCOREMATH double asum<double>( hoNDArray<double>* );
    template EXPORTCPUCOREMATH double nrm2<double>( hoNDArray<double>* );

    template EXPORTCPUCOREMATH size_t amin<double>( hoNDArray<double>* );
    template EXPORTCPUCOREMATH size_t amax<double>( hoNDArray<double>* );

    template EXPORTCPUCOREMATH std::complex<float> dot< std::complex<float> >( hoNDArray< std::complex<float> >*, hoNDArray< std::complex<float> >*, bool );
    template EXPORTCPUCOREMATH float asum<float>( hoNDArray< std::complex<float> >* );
    template EXPORTCPUCOREMATH float nrm2< std::complex<float> >( hoNDArray< std::complex<float> >* );
    template EXPORTCPUCOREMATH float nrm1< std::complex<float> >( hoNDArray< std::complex<float> >* );
    template EXPORTCPUCOREMATH size_t amin<float>( hoNDArray< std::complex<float> >* );
    template EXPORTCPUCOREMATH size_t amax<float>( hoNDArray< std::complex<float> >* );

    template EXPORTCPUCOREMATH std::complex<double> dot< std::complex<double> >( hoNDArray< std::complex<double> >*, hoNDArray< std::complex<double> >*, bool );
    template EXPORTCPUCOREMATH double asum<double>( hoNDArray< std::complex<double> >* );
    template EXPORTCPUCOREMATH double nrm2< std::complex<double> >( hoNDArray< std::complex<double> >* );
    template EXPORTCPUCOREMATH double nrm1< std::complex<double> >( hoNDArray< std::complex<double> >* );
    template EXPORTCPUCOREMATH size_t amin<double>( hoNDArray< std::complex<double> >* );
    template EXPORTCPUCOREMATH size_t amax<double>( hoNDArray< std::complex<double> >* );

    template EXPORTCPUCOREMATH complext<float> dot< complext<float> >( hoNDArray< complext<float> >*, hoNDArray< complext<float> >*, bool );
    template EXPORTCPUCOREMATH float asum<float>( hoNDArray< complext<float> >* );
    template EXPORTCPUCOREMATH float nrm2< complext<float> >( hoNDArray< complext<float> >* );
    template EXPORTCPUCOREMATH float nrm1< complext<float> >( hoNDArray< complext<float> >* );
    template EXPORTCPUCOREMATH size_t amin<float>( hoNDArray< complext<float> >* );
    template EXPORTCPUCOREMATH size_t amax<float>( hoNDArray< complext<float> >* );

    template EXPORTCPUCOREMATH complext<double> dot< complext<double> >( hoNDArray< complext<double> >*, hoNDArray< complext<double> >*, bool );
    template EXPORTCPUCOREMATH double asum<double>( hoNDArray< complext<double> >* );
    template EXPORTCPUCOREMATH double nrm2< complext<double> >( hoNDArray< complext<double> >* );
    template EXPORTCPUCOREMATH size_t amin<double>( hoNDArray< complext<double> >* );
    template EXPORTCPUCOREMATH size_t amax<double>( hoNDArray< complext<double> >* );

    // --------------------------------------------------------------------------------

    inline void dotc(size_t N, const  std::complex<float> * x, const  std::complex<float> * y,  std::complex<float> & r)
    {
        long long n;

        float sa(0), sb(0);

        #pragma omp parallel for private(n) reduction(+:sa) reduction(+:sb) if (N>NumElementsUseThreading)
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

        double sa(0), sb(0);

        #pragma omp parallel for private(n) reduction(+:sa) reduction(+:sb) if (N>NumElementsUseThreading)
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

    inline void dotc(size_t N, const  complext<float> * x, const  complext<float> * y,  complext<float> & r)
    {
        long long n;

        float sa(0), sb(0);

        #pragma omp parallel for private(n) reduction(+:sa) reduction(+:sb) if (N>NumElementsUseThreading)
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

    inline void dotc(size_t N, const  complext<double> * x, const  complext<double> * y,  complext<double> & r)
    {
        long long n;

        double sa(0), sb(0);

        #pragma omp parallel for private(n) reduction(+:sa) reduction(+:sb) if (N>NumElementsUseThreading)
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
    void dotc(const hoNDArray<T>& x, const hoNDArray<T>& y, T& r)
    {
        GADGET_DEBUG_CHECK_THROW(x.get_number_of_elements()==y.get_number_of_elements());
        dotc(x.get_number_of_elements(), x.begin(), y.begin(), r);
    }

    template EXPORTCPUCOREMATH void dotc(const hoNDArray< std::complex<float> >& x, const hoNDArray< std::complex<float> >& y,  std::complex<float> & r);
    template EXPORTCPUCOREMATH void dotc(const hoNDArray< std::complex<double> >& x, const hoNDArray< std::complex<double> >& y,  std::complex<double> & r);
    template EXPORTCPUCOREMATH void dotc(const hoNDArray< complext<float> >& x, const hoNDArray< complext<float> >& y,  complext<float> & r);
    template EXPORTCPUCOREMATH void dotc(const hoNDArray< complext<double> >& x, const hoNDArray< complext<double> >& y,  complext<double> & r);

    template <typename T> 
    T dotc(const hoNDArray<T>& x, const hoNDArray<T>& y)
    {
        T r;
        dotc(x, y, r);
        return r;
    }

    template EXPORTCPUCOREMATH std::complex<float> dotc(const hoNDArray< std::complex<float> >& x, const hoNDArray< std::complex<float> >& y);
    template EXPORTCPUCOREMATH std::complex<double> dotc(const hoNDArray< std::complex<double> >& x, const hoNDArray< std::complex<double> >& y);
    template EXPORTCPUCOREMATH complext<float> dotc(const hoNDArray< complext<float> >& x, const hoNDArray< complext<float> >& y);
    template EXPORTCPUCOREMATH complext<double> dotc(const hoNDArray< complext<double> >& x, const hoNDArray< complext<double> >& y);

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

    inline void dotu(size_t N, const  complext<float> * x, const  complext<float> * y,  complext<float> & r)
    {
        long long n;

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

    inline void dotu(size_t N, const  complext<double> * x, const  complext<double> * y,  complext<double> & r)
    {
        long long n;

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
    void dotu(const hoNDArray<T>& x, const hoNDArray<T>& y, T& r)
    {
        GADGET_DEBUG_CHECK_THROW(x.get_number_of_elements()==y.get_number_of_elements());
        dotu(x.get_number_of_elements(), x.begin(), y.begin(), r);
    }

    template EXPORTCPUCOREMATH void dotu(const hoNDArray<float>& x, const hoNDArray<float>& y, float& r);
    template EXPORTCPUCOREMATH void dotu(const hoNDArray<double>& x, const hoNDArray<double>& y, double& r);
    template EXPORTCPUCOREMATH void dotu(const hoNDArray< std::complex<float> >& x, const hoNDArray< std::complex<float> >& y, std::complex<float>& r);
    template EXPORTCPUCOREMATH void dotu(const hoNDArray< std::complex<double> >& x, const hoNDArray< std::complex<double> >& y, std::complex<double>& r);
    template EXPORTCPUCOREMATH void dotu(const hoNDArray< complext<float> >& x, const hoNDArray< complext<float> >& y, complext<float>& r);
    template EXPORTCPUCOREMATH void dotu(const hoNDArray< complext<double> >& x, const hoNDArray< complext<double> >& y, complext<double>& r);

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
    template EXPORTCPUCOREMATH  complext<float>  dotu(const hoNDArray< complext<float> >& x, const hoNDArray< complext<float> >& y);
    template EXPORTCPUCOREMATH  complext<double>  dotu(const hoNDArray< complext<double> >& x, const hoNDArray< complext<double> >& y);

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
    void sort(const hoNDArray<T>& x, hoNDArray<T>& r, bool isascending)
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
    }

    template EXPORTCPUCOREMATH void sort(const hoNDArray<float>& x, hoNDArray<float>& r, bool isascending);
    template EXPORTCPUCOREMATH void sort(const hoNDArray<double>& x, hoNDArray<double>& r, bool isascending);

    // --------------------------------------------------------------------------------

    template <typename T>
    struct hoCompAscendingIndex
    {
        typedef std::pair<size_t, T> PairType;
        bool operator() (const PairType& a, const PairType& b) { return (a.second < b.second); }
    };

    template <typename T>
    struct hoCompDescendingIndex
    {
        typedef std::pair<size_t, T> PairType; 
        bool operator() (const PairType& a, const PairType& b) { return (a.second >= b.second); }
    };

    template <typename T>
    void sort(size_t N, const T* x, T* r, std::vector<size_t>& ind, bool isascending)
    {
        if (r != x)
        {
            memcpy(r, x, sizeof(T)*N);
        }

        ind.resize(N, 0);

        std::vector< std::pair<size_t, T> > x_v(N);

        size_t n;
        for (n = 0; n < N; n++)
        {
            x_v[n].first = n;
            x_v[n].second = x[n];
        }

        if (isascending)
        {
            hoCompAscendingIndex<T> obj;
            std::sort(x_v.begin(), x_v.end(), obj);
        }
        else
        {
            hoCompDescendingIndex<T> obj;
            std::sort(x_v.begin(), x_v.end(), obj);
        }

        for (n = 0; n < N; n++)
        {
            ind[n] = x_v[n].first;
            r[n] = x_v[n].second;
        }
    }

    template <typename T>
    void sort(const hoNDArray<T>& x, hoNDArray<T>& r, std::vector<size_t>& ind, bool isascending)
    {
        if (&r != &x)
        {
            if (r.get_number_of_elements() != x.get_number_of_elements())
            {
                r = x;
            }
            else
            {
                memcpy(r.begin(), x.begin(), x.get_number_of_bytes());
            }
        }

        sort(x.get_number_of_elements(), x.begin(), r.begin(), ind, isascending);
    }

    template EXPORTCPUCOREMATH void sort(const hoNDArray<float>& x, hoNDArray<float>& r, std::vector<size_t>& ind, bool isascending);
    template EXPORTCPUCOREMATH void sort(const hoNDArray<double>& x, hoNDArray<double>& r, std::vector<size_t>& ind, bool isascending);

    // --------------------------------------------------------------------------------

    template <class T>
    void minValue(const hoNDArray<T>& a, T& v)
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
            GADGET_THROW("Errors in minValue(const hoNDArray<T>& a, T& v) ... ");
        }
    }

    template EXPORTCPUCOREMATH void minValue(const hoNDArray<float>& a, float& v);
    template EXPORTCPUCOREMATH void minValue(const hoNDArray<double>& a, double& v);

    template <class T>
    void maxValue(const hoNDArray<T>& a, T& v)
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
            GADGET_THROW("Errors in maxValue(const hoNDArray<T>& a, T& v) ... ");
        }
    }

    template EXPORTCPUCOREMATH void maxValue(const hoNDArray<float>& a, float& v);
    template EXPORTCPUCOREMATH void maxValue(const hoNDArray<double>& a, double& v);

    // --------------------------------------------------------------------------------
}
