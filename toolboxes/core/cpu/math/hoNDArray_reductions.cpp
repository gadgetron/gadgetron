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



namespace Gadgetron{

    // --------------------------------------------------------------------------------

    template<class REAL> REAL max(const hoNDArray<REAL>* data){
        return as_arma_col(data).max();
    }

    // --------------------------------------------------------------------------------

    template<class REAL> REAL min(const hoNDArray<REAL>* data){
        return as_arma_col(data).min();
    }

    // --------------------------------------------------------------------------------

    template<class T> T mean(const hoNDArray<T>* data){
        return (typename stdType<T>::Type) arma::mean(as_arma_col(data));
    }

    // --------------------------------------------------------------------------------

    template<class T> T sum(const hoNDArray<T>* data){
        return (typename stdType<T>::Type) arma::sum(as_arma_col(data));
    }

    // --------------------------------------------------------------------------------

    template<class T> T stddev(const hoNDArray<T>* data){
        return (typename stdType<T>::Type) arma::stddev(as_arma_col(data));
    }

    // --------------------------------------------------------------------------------

    template<class T> T var(const hoNDArray<T>* data) {
        return (typename stdType<T>::Type) arma::var(as_arma_col(data));
    }

    // --------------------------------------------------------------------------------

     template<class T> T median(const hoNDArray<T>* data){
        return (typename stdType<T>::Type) arma::median(as_arma_col(data));
    }

    // --------------------------------------------------------------------------------

    // --------------------------------------------------------------------------------


    // --------------------------------------------------------------------------------


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

    template<class T> size_t amin(const  hoNDArray<T> *x )
    {
        if( x == 0x0 )
            throw std::runtime_error("Gadgetron::amin(): Invalid input array");

        typedef typename realType<T>::Type realT;
        arma::Col<realT> xM = arma::abs(as_arma_col(x));
        arma::uword idx;
        realT min = xM.min(idx);
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

    template EXPORTCPUCOREMATH float max(const hoNDArray<float>*);
    template EXPORTCPUCOREMATH float min(const hoNDArray<float>*);
    template EXPORTCPUCOREMATH float mean(const hoNDArray<float>*);
    template EXPORTCPUCOREMATH float median(const hoNDArray<float>*);
    template EXPORTCPUCOREMATH float sum(const hoNDArray<float>*);
    template EXPORTCPUCOREMATH float stddev(const hoNDArray<float>*);
    template EXPORTCPUCOREMATH float var(const hoNDArray<float>*);

    template EXPORTCPUCOREMATH double max(const hoNDArray<double>*);
    template EXPORTCPUCOREMATH double min(const hoNDArray<double>*);
    template EXPORTCPUCOREMATH double mean(const hoNDArray<double>*);
    template EXPORTCPUCOREMATH double median(const hoNDArray<double>*);
    template EXPORTCPUCOREMATH double sum(const hoNDArray<double>*);
    template EXPORTCPUCOREMATH double stddev(const hoNDArray<double>*);
    template EXPORTCPUCOREMATH double var(const hoNDArray<double>*);

    template EXPORTCPUCOREMATH complext<double> mean(const hoNDArray<complext<double> >*);
    template EXPORTCPUCOREMATH complext<double> median(const hoNDArray<complext<double> >*);
    template EXPORTCPUCOREMATH complext<double> sum(const hoNDArray<complext<double> >*);
    template EXPORTCPUCOREMATH complext<double> stddev(const hoNDArray<complext<double> >*);
    template EXPORTCPUCOREMATH complext<double> var(const hoNDArray<complext<double> >*);

    template EXPORTCPUCOREMATH complext<float> mean(const hoNDArray<complext<float> >*);
    template EXPORTCPUCOREMATH complext<float> median(const hoNDArray<complext<float> >*);
    template EXPORTCPUCOREMATH complext<float> sum(const hoNDArray<complext<float> >*);
    template EXPORTCPUCOREMATH complext<float> stddev(const hoNDArray<complext<float> >*);
    template EXPORTCPUCOREMATH complext<float> var(const hoNDArray<complext<float> >*);

    template EXPORTCPUCOREMATH std::complex<double> mean(const hoNDArray<std::complex<double> >*);
    template EXPORTCPUCOREMATH std::complex<double> sum(const hoNDArray<std::complex<double> >*);
    template EXPORTCPUCOREMATH std::complex<double> stddev(const hoNDArray<std::complex<double> >*);
    template EXPORTCPUCOREMATH std::complex<double> var(const hoNDArray<std::complex<double> >*);

    template EXPORTCPUCOREMATH std::complex<float> mean(const hoNDArray<std::complex<float> >*);
    template EXPORTCPUCOREMATH std::complex<float> sum(const hoNDArray<std::complex<float> >*);
    template EXPORTCPUCOREMATH std::complex<float> stddev(const hoNDArray<std::complex<float> >*);
    template EXPORTCPUCOREMATH std::complex<float> var(const hoNDArray<std::complex<float> >*);


    template EXPORTCPUCOREMATH size_t amin<std::complex<float>>( const hoNDArray< std::complex<float> >* );
    template EXPORTCPUCOREMATH size_t amin<std::complex<double>>( const hoNDArray< std::complex<double> >* );

    template EXPORTCPUCOREMATH size_t amin<complext<float>>( const hoNDArray< complext<float> >* );
    template EXPORTCPUCOREMATH size_t amin<complext<double>>( const hoNDArray< complext<double> >* );

    template EXPORTCPUCOREMATH size_t amin<float>( const hoNDArray< float >* );
    template EXPORTCPUCOREMATH size_t amin<double>( const hoNDArray< double >* );




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



    // --------------------------------------------------------------------------------
}
