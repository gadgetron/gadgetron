#include "hoNDArray_elemwise.h"
#include "hoNDArray_operators.h"
#include "hoNDArray_blas.h"
#include "complext.h"
#include "hoArmadillo.h"

#ifdef USE_OMP
#include <omp.h>
#endif

namespace Gadgetron{

    template<class T> boost::shared_ptr< hoNDArray<typename realType<T>::Type> > abs( hoNDArray<T> *x )
    { 
        if( x == 0x0 )
            throw std::runtime_error("Gadgetron::abs(): Invalid input array");

        boost::shared_ptr< hoNDArray<typename realType<T>::Type> > result(new hoNDArray<typename realType<T>::Type>());
        result->create(x->get_dimensions());
        arma::Col<typename realType<T>::Type> aRes = as_arma_col(result.get());
        aRes = arma::abs(as_arma_col(x));
        return result;
    }

    template<class T> void abs_inplace( hoNDArray<T> *x )
    { 
        if( x == 0x0 )
            throw std::runtime_error("Gadgetron::abs_inplace(): Invalid input array");

        arma::Col<typename realType<T>::Type> aRes = as_arma_col(x);
        aRes = arma::abs(aRes);
    }  

    template<class T> boost::shared_ptr< hoNDArray<typename realType<T>::Type> > abs_square( hoNDArray<T> *x )
    { 
        if( x == 0x0 )
            throw std::runtime_error("Gadgetron::abs_square(): Invalid input array");

        boost::shared_ptr< hoNDArray<typename realType<T>::Type> > result(new hoNDArray<typename realType<T>::Type>());
        result->create(x->get_dimensions());
        arma::Col<typename realType<T>::Type> aRes = as_arma_col(result.get());
        aRes = arma::square(abs(as_arma_col(x)));
        return result;
    }

    template<class T> boost::shared_ptr< hoNDArray<T> > sqrt( hoNDArray<T> *x )
    { 
        if( x == 0x0 )
            throw std::runtime_error("Gadgetron::sqrt(): Invalid input array");

        boost::shared_ptr< hoNDArray<T> > result(new hoNDArray<T>());
        result->create(x->get_dimensions());
        arma::Col<typename stdType<T>::Type> aRes = as_arma_col(result.get());
        aRes = arma::sqrt(as_arma_col(x));
        return result;
    }

    template<class T> void sqrt_inplace( hoNDArray<T> *x )
    { 
        if( x == 0x0 )
            throw std::runtime_error("Gadgetron::sqrt_inplace(): Invalid input array");

        arma::Col<typename stdType<T>::Type> aRes = as_arma_col(x);
        aRes = arma::sqrt(aRes);
    }

    template<class T> boost::shared_ptr< hoNDArray<T> > square( hoNDArray<T> *x )
    { 
        if( x == 0x0 )
            throw std::runtime_error("Gadgetron::square(): Invalid input array");

        boost::shared_ptr< hoNDArray<T> > result(new hoNDArray<T>());
        result->create(x->get_dimensions());
        arma::Col<typename stdType<T>::Type> aRes = as_arma_col(result.get());
        aRes = arma::square(as_arma_col(x));
        return result;
    }

    template<class T> void square_inplace( hoNDArray<T> *x )
    { 
        if( x == 0x0 )
            throw std::runtime_error("Gadgetron::square_inplace(): Invalid input array");

        arma::Col<typename stdType<T>::Type> aRes = as_arma_col(x);
        aRes = arma::square(aRes);
    }  

    template<class T> boost::shared_ptr< hoNDArray<T> > reciprocal( hoNDArray<T> *x )
    { 
        if( x == 0x0 )
            throw std::runtime_error("Gadgetron::reciprocal(): Invalid input array");

        arma::Col<typename stdType<T>::Type> ones(x->get_number_of_elements());
        ones.ones();
        boost::shared_ptr< hoNDArray<T> > result(new hoNDArray<T>());
        result->create(x->get_dimensions());
        arma::Col<typename stdType<T>::Type> aRes = as_arma_col(result.get());
        aRes = ones/as_arma_col(x);
        return result;
    }

    template<class T> void reciprocal_inplace( hoNDArray<T> *x )
    { 
        if( x == 0x0 )
            throw std::runtime_error("Gadgetron::reciprocal_inplace(): Invalid input array");

        arma::Col<typename stdType<T>::Type> aRes = as_arma_col(x);
        arma::Col<typename stdType<T>::Type> ones(x->get_number_of_elements());
        ones.ones();
        aRes = ones/aRes;
    }

    template<class T> boost::shared_ptr< hoNDArray<T> > reciprocal_sqrt( hoNDArray<T> *x )
    { 
        if( x == 0x0 )
            throw std::runtime_error("Gadgetron::reciprocal_sqrt(): Invalid input array");

        arma::Col<typename stdType<T>::Type> ones(x->get_number_of_elements());
        ones.ones();   
        boost::shared_ptr< hoNDArray<T> > result(new hoNDArray<T>());
        result->create(x->get_dimensions());
        arma::Col<typename stdType<T>::Type> aRes = as_arma_col(result.get());
        aRes = ones/arma::sqrt(as_arma_col(x));
        return result;
    }

    template<class T> void reciprocal_sqrt_inplace( hoNDArray<T> *x )
    { 
        if( x == 0x0 )
            throw std::runtime_error("Gadgetron::reciprocal_sqrt_inplace(): Invalid input array");

        arma::Col<typename stdType<T>::Type> ones(x->get_number_of_elements());
        ones.ones();
        arma::Col<typename stdType<T>::Type> aRes = as_arma_col(x);
        aRes = ones/arma::sqrt(aRes);
    }

    template<class T> boost::shared_ptr< hoNDArray<T> > sgn( hoNDArray<T> *x )
    {
        if( x == 0x0 )
            throw std::runtime_error("Gadgetron::sgn(): Invalid input array");

        boost::shared_ptr< hoNDArray<T> > res( new hoNDArray<T>() );
        res->create(x->get_dimensions());   
#ifdef USE_OMP
#pragma omp parallel for
#endif
        for( long long i = 0; i < res->get_number_of_elements(); i++ ){
            res->get_data_ptr()[i] = sgn(x->get_data_ptr()[i]);
        }
        return res;
    }

    template<class T> void sgn_inplace( hoNDArray<T> *x )
    {
        if( x == 0x0 )
            throw std::runtime_error("Gadgetron::sgn_inplace(): Invalid input array");

#ifdef USE_OMP
#pragma omp parallel for
#endif
        for( long long i = 0; i < x->get_number_of_elements(); i++ ) 
            x->get_data_ptr()[i] = sgn(x->get_data_ptr()[i]);
    }

    template<class T> boost::shared_ptr< hoNDArray<typename realType<T>::Type> > real( hoNDArray<T> *x )
    {
        if( x == 0x0 )
            throw std::runtime_error("Gadgetron::real(): Invalid input array");

        boost::shared_ptr< hoNDArray<typename realType<T>::Type> > result(new hoNDArray<typename realType<T>::Type>());
        result->create(x->get_dimensions());
        arma::Col<typename realType<T>::Type> aRes = as_arma_col(result.get());
        aRes = arma::real(as_arma_col(x));
        return result;
    }

    template<class T> boost::shared_ptr< hoNDArray<typename realType<T>::Type> > imag( hoNDArray<T> *x )
    {
        if( x == 0x0 )
            throw std::runtime_error("Gadgetron::imag(): Invalid input array");

        boost::shared_ptr< hoNDArray<typename realType<T>::Type> > result(new hoNDArray<typename realType<T>::Type>());
        result->create(x->get_dimensions());
        arma::Col<typename realType<T>::Type> aRes = as_arma_col(result.get());
        aRes = arma::imag(as_arma_col(x));
        return result;
    }

    template<class T> boost::shared_ptr< hoNDArray<T> > conj( hoNDArray<T> *x )
    {
        if( x == 0x0 )
            throw std::runtime_error("Gadgetron::conj(): Invalid input array");

        boost::shared_ptr< hoNDArray<T> > result(new hoNDArray<T>());
        result->create(x->get_dimensions());
        arma::Col<typename stdType<T>::Type> aRes = as_arma_col(result.get());
        aRes = arma::conj(as_arma_col(x));
        return result;
    }

    template<class T> boost::shared_ptr< hoNDArray<T> > real_to_complex( hoNDArray<typename realType<T>::Type> *x )
    {
        if( x == 0x0 )
            BOOST_THROW_EXCEPTION(runtime_error("Gadgetron::real_to_complex(): Invalid input array"));

        boost::shared_ptr< hoNDArray<T> > result(new hoNDArray<T>());
        result->create(x->get_dimensions());
        arma::Col<typename stdType<T>::Type> aRes = as_arma_col(result.get());
        aRes = arma::Col<typename stdType<T>::Type>(as_arma_col(x), arma::Col<typename realType<T>::Type>(x->get_number_of_elements()).zeros());
        return result;
    }

    template<class T> boost::shared_ptr< hoNDArray<T> > real_imag_to_complex( hoNDArray<typename realType<T>::Type>* real, hoNDArray<typename realType<T>::Type>* imag )
    {
        if( real==0x0 || imag==0x0 )
            BOOST_THROW_EXCEPTION(runtime_error("Gadgetron::real_imag_to_complex(): Invalid input array"));

        if( real->get_number_of_elements() != imag->get_number_of_elements() )
            BOOST_THROW_EXCEPTION(runtime_error("Gadgetron::real_imag_to_complex(): Invalid input array"));

        boost::shared_ptr< hoNDArray<T> > result(new hoNDArray<T>());
        result->create(real->get_dimensions());

        T* pRes = result->begin();

        size_t N = real->get_number_of_elements();
        for ( size_t n=0; n<N; n++ )
        {
            pRes[n] = T(real->at(n), imag->at(n));
        }

        return result;
    }

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
            for ( n=0; n<N; n++ )
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
            for ( n=0; n<N; n++ )
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

    template<> 
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
            for ( n=0; n<N; n++ )
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

    template<> 
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
            for ( n=0; n<N; n++ )
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
            for ( n=0; n<N; n++ )
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

    template<class T> inline void clear( hoNDArray<T> *x )
    {
        if( x == 0x0 )
            throw std::runtime_error("Gadgetron::clear(): Invalid input array");

        if ( x->get_number_of_elements() > 0 )
        {
            memset( x->get_data_ptr(), 0, x->get_number_of_elements()*sizeof(T));
        }
    }

    template<class T> inline void clear( hoNDArray<T>& x )
    {
        if ( x.get_number_of_elements() > 0 )
        {
            memset( x.get_data_ptr(), 0, x.get_number_of_elements()*sizeof(T));
        }
    }

    template<class T> void fill( hoNDArray<T> *x, T val )
    {
        if( x == 0x0 )
            throw std::runtime_error("Gadgetron::fill(): Invalid input array");

        arma::Col<typename stdType<T>::Type> aRes = as_arma_col(x);
        aRes.fill(*((typename stdType<T>::Type*)&val));
    }

    //
    // TODO:
    // The clamp functions could (probably) be implemented much like we use Thrust for the device versions
    // - i.e. using Armadillo's transform on the array.
    // However this requires a newer version of Armadillo as current Linux distributions provide...
    //

    template<typename T> struct hoNDA_clamp //: public thrust::unary_function<T,T>
    {
      hoNDA_clamp( T _min, T _max, T _min_val, T _max_val ) : min(_min), max(_max), min_val(_min_val), max_val(_max_val) {}
        T operator()(const T &x) const 
        {
            if( x < min ) return min_val;
            else if ( x >= max) return max_val;
            else return x;
        }
      T min, max;
      T min_val, max_val;
    };

    template<typename T> struct hoNDA_clamp< std::complex<T> > //: public thrust::unary_function< std::complex<T>, std::complex<T> >
    {
      hoNDA_clamp( T _min, T _max, std::complex<T> _min_val, std::complex<T> _max_val ) : min(_min), max(_max), min_val(_min_val), max_val(_max_val) {}
        std::complex<T> operator()(const std::complex<T> &x) const 
        {
            if( real(x) < min ) return min_val;
            else if ( real(x) >= max) return max_val;
            else return std::complex<T>(real(x));
        }
      T min, max;
      std::complex<T> min_val, max_val;
    };

    template<typename T> struct hoNDA_clamp< complext<T> > //: public thrust::unary_function< complext<T>, complext<T> >
    {
        hoNDA_clamp( T _min, T _max, complext<T> _min_val, complext<T> _max_val ) : min(_min), max(_max), min_val(_min_val), max_val(_max_val) {}
        complext<T> operator()(const complext<T> &x) const 
        {
            if( real(x) < min ) return min_val;
            else if ( real(x) >= max) return max_val;
            else return complext<T>(real(x));
        }
        T min, max;
        complext<T> min_val, max_val;
    };

    template<class T> void clamp( hoNDArray<T> *x, 
                                  typename realType<T>::Type min, typename realType<T>::Type max, T min_val, T max_val )
    { 
        if( x == 0x0 )
            throw std::runtime_error("Gadgetron::clamp(): Invalid input array");

        hoNDA_clamp<T> functor(min, max, min_val, max_val);
        std::transform(x->begin(),x->end(),x->begin(),functor);
    }  

    template<class T> void clamp( hoNDArray<T> *x, typename realType<T>::Type min, typename realType<T>::Type max )
    {
        clamp(x,min,max,T(min),T(max));
    }

    template<typename T> struct hoNDA_clamp_min //: public thrust::unary_function<T,T>
    {
        hoNDA_clamp_min( T _min ) : min(_min) {}
        T operator()(const T &x) const 
        {
            if( x < min ) return min;
            else return x;
        }
        T min;
    };

    template<typename T> struct hoNDA_clamp_min< std::complex<T> > //: public thrust::unary_function< std::complex<T>, std::complex<T> >
    {
        hoNDA_clamp_min( T _min ) : min(_min) {}
        std::complex<T> operator()(const std::complex<T> &x) const 
        {
            if( real(x) < min ) return std::complex<T>(min);
            else return std::complex<T>(real(x));
        }
        T min;
    };

    template<typename T> struct hoNDA_clamp_min< complext<T> > //: public thrust::unary_function< complext<T>, complext<T> >
    {
        hoNDA_clamp_min( T _min ) : min(_min) {}
        complext<T> operator()(const complext<T> &x) const 
        {
            if( real(x) < min ) return complext<T>(min);
            else return complext<T>(real(x));
        }
        T min;
    };

    template<class T> void clamp_min( hoNDArray<T> *x, typename realType<T>::Type min )
    { 
        if( x == 0x0 )
            throw std::runtime_error("Gadgetron::clamp_min(): Invalid input array");

        hoNDA_clamp_min<T> functor(min);
        std::transform(x->begin(),x->end(),x->begin(),functor);
    }  

    template<typename T> struct hoNDA_clamp_max //: public thrust::unary_function<T,T>
    {
        hoNDA_clamp_max( T _max ) : max(_max) {}
        T operator()(const T &x) const 
        {
            if( x > max ) return max;
            else return x;
        }
        T max;
    };

    template<typename T> struct hoNDA_clamp_max< std::complex<T> > //: public thrust::unary_function< std::complex<T>, std::complex<T> >
    {
        hoNDA_clamp_max( T _max ) : max(_max) {}
        std::complex<T> operator()(const std::complex<T> &x) const 
        {
            if( real(x) > max ) return std::complex<T>(max);
            else return std::complex<T>(real(x));
        }
        T max;
    };

    template<typename T> struct hoNDA_clamp_max< complext<T> > //: public thrust::unary_function< complext<T>, complext<T> >
    {
        hoNDA_clamp_max( T _max ) : max(_max) {}
        complext<T> operator()(const complext<T> &x) const 
        {
            if( real(x) > max ) return complext<T>(max);
            else return complext<T>(real(x));
        }
        T max;
    };

    template<class T> void clamp_max( hoNDArray<T> *x, typename realType<T>::Type max )
    { 
        if( x == 0x0 )
            throw std::runtime_error("Gadgetron::clamp_max(): Invalid input array");

        hoNDA_clamp_max<T> functor(max);
        std::transform(x->begin(),x->end(),x->begin(),functor);
    }

    template<class T> void normalize( hoNDArray<T> *x, typename realType<T>::Type val )
    {
        if( x == 0x0 )
            throw std::runtime_error("Gadgetron::normalize(): Invalid input array");

        size_t max_idx = amax(x);
        T max_val_before = x->get_data_ptr()[max_idx];
        typename realType<T>::Type scale = val/abs(max_val_before);
        *x *= scale;
    }

    template<class T> void shrink1( hoNDArray<T> *x, typename realType<T>::Type gamma, hoNDArray<T> *out )
    {
        if( x == 0x0 )
            throw std::runtime_error("Gadgetron::shrink1(): Invalid input array");

        T *outPtr = (out==0x0) ? x->get_data_ptr() : out->get_data_ptr();

#ifdef USE_OMP
#pragma omp parallel for
#endif
        for( long long i = 0; i < x->get_number_of_elements(); i++ ) {
            T prev = x->get_data_ptr()[i];
            typename realType<T>::Type absPrev = abs(prev);
            T sgnPrev = (absPrev <= typename realType<T>::Type(0)) ? T(0) : prev/absPrev;
            outPtr[i] = sgnPrev*std::max(absPrev-gamma, typename realType<T>::Type(0));
        } 
    }

    template<class T> void pshrink( hoNDArray<T> *x, typename realType<T>::Type gamma,typename realType<T>::Type p, hoNDArray<T> *out )
    {
        if( x == 0x0 )
            throw std::runtime_error("Gadgetron::pshrink(): Invalid input array");

        T *outPtr = (out==0x0) ? x->get_data_ptr() : out->get_data_ptr();

#ifdef USE_OMP
#pragma omp parallel for
#endif
        for( long long i = 0; i < x->get_number_of_elements(); i++ ) {
            T prev = x->get_data_ptr()[i];
            typename realType<T>::Type absPrev = abs(prev);
            T sgnPrev = (absPrev <= typename realType<T>::Type(0)) ? T(0) : prev/absPrev;
            outPtr[i] = sgnPrev*std::max(absPrev-gamma*std::pow(absPrev,p-1), typename realType<T>::Type(0));
        }
    }

    template<class T> void shrinkd ( hoNDArray<T> *_x, hoNDArray<typename realType<T>::Type> *_s, typename realType<T>::Type gamma, hoNDArray<T> *out )
    {
        if( _x == 0x0  || _s == 0 )
            throw std::runtime_error("Gadgetron::shrinkd(): Invalid input array");

        T *outPtr = (out==0x0) ? _x->get_data_ptr() : out->get_data_ptr();

#ifdef USE_OMP
#pragma omp parallel for
#endif
        for( long long i = 0; i < _x->get_number_of_elements(); i++ ) {
            T x = _x->get_data_ptr()[i];
            typename realType<T>::Type s = _s->get_data_ptr()[i];
            if (s > gamma)
            	outPtr[i] = x/s*(s-gamma);
            else
            	outPtr[i] = 0;
        } 
    }

    template<class T> void pshrinkd( hoNDArray<T> *_x, hoNDArray<typename realType<T>::Type> *_s, typename realType<T>::Type gamma,typename realType<T>::Type p, hoNDArray<T> *out )
    {
        if( _x == 0x0 )
            throw std::runtime_error("Gadgetron::pshrinkd(): Invalid input array");

        T *outPtr = (out==0x0) ? _x->get_data_ptr() : out->get_data_ptr();

#ifdef USE_OMP
#pragma omp parallel for
#endif
        for( long long i = 0; i < _x->get_number_of_elements(); i++ )
        {
            T x = _x->get_data_ptr()[i];
            typename realType<T>::Type s = _s->get_data_ptr()[i];
            outPtr[i] = x/s*std::max(s-gamma*std::pow(s,p-1),typename realType<T>::Type(0));
        }
    }

    #ifdef USE_MKL

    // ----------------------------------------------------------------------------------------
    // float
    // ----------------------------------------------------------------------------------------

    bool add(const hoNDArray<float>& x, const hoNDArray<float>& y, hoNDArray<float>& r)
    {
        GADGET_DEBUG_CHECK_RETURN_FALSE(x.get_number_of_elements()==y.get_number_of_elements());
        if ( r.get_number_of_elements()!=x.get_number_of_elements())
        {
            r = x;
        }

        vsAdd(x.get_number_of_elements(), x.begin(), y.begin(), r.begin());
        GADGET_CHECK_RETURN_FALSE(vmlGetErrStatus()==0);

        return true;
    }

    bool subtract(const hoNDArray<float>& x, const hoNDArray<float>& y, hoNDArray<float>& r)
    {
        GADGET_DEBUG_CHECK_RETURN_FALSE(x.get_number_of_elements()==y.get_number_of_elements());
        if ( r.get_number_of_elements()!=x.get_number_of_elements())
        {
            r = x;
        }

        vsSub(x.get_number_of_elements(), x.begin(), y.begin(), r.begin());
        GADGET_CHECK_RETURN_FALSE(vmlGetErrStatus()==0);

        return true;
    }

    bool multiply(const hoNDArray<float>& x, const hoNDArray<float>& y, hoNDArray<float>& r)
    {
        GADGET_DEBUG_CHECK_RETURN_FALSE(x.get_number_of_elements()==y.get_number_of_elements());
        if ( r.get_number_of_elements()!=x.get_number_of_elements())
        {
            r = x;
        }

        vsMul(x.get_number_of_elements(), x.begin(), y.begin(), r.begin());
        GADGET_CHECK_RETURN_FALSE(vmlGetErrStatus()==0);

        return true;
    }

    bool divide(const hoNDArray<float>& x, const hoNDArray<float>& y, hoNDArray<float>& r)
    {
        GADGET_DEBUG_CHECK_RETURN_FALSE(x.get_number_of_elements()==y.get_number_of_elements());
        if ( r.get_number_of_elements()!=x.get_number_of_elements())
        {
            r = x;
        }

        vsDiv(x.get_number_of_elements(), x.begin(), y.begin(), r.begin());
        GADGET_CHECK_RETURN_FALSE(vmlGetErrStatus()==0);

        return true;
    }

    bool absolute(const hoNDArray<float>& x, hoNDArray<float>& r)
    {
        if ( r.get_number_of_elements()!=x.get_number_of_elements())
        {
            r = x;
        }

        vsAbs(x.get_number_of_elements(), x.begin(), r.begin());
        GADGET_CHECK_RETURN_FALSE(vmlGetErrStatus()==0);

        return true;
    }

    bool argument(const hoNDArray<float>& x, hoNDArray<float>& r)
    {
        if ( r.get_number_of_elements()!=x.get_number_of_elements())
        {
            r = x;
        }

        memset(r.begin(), 0, r.get_number_of_bytes());

        return true;
    }

    bool sqrt(const hoNDArray<float>& x, hoNDArray<float>& r)
    {
        if ( r.get_number_of_elements()!=x.get_number_of_elements())
        {
            r = x;
        }

        vsSqrt(x.get_number_of_elements(), x.begin(), r.begin());
        GADGET_CHECK_RETURN_FALSE(vmlGetErrStatus()==0);

        return true;
    }

    bool minAbsolute(const hoNDArray<float>& x, float& r, size_t& ind)
    {
        try
        {
            MKL_INT n = x.get_number_of_elements();
            MKL_INT incx = 1;
            ind = (size_t)(isamin(&n, x.begin(), &incx));
            r = x.at(ind);
        }
        catch(...)
        {
            return false;
        }

        return true;
    }

    bool maxAbsolute(const hoNDArray<float>& x, float& r, size_t& ind)
    {
        try
        {
            MKL_INT n = x.get_number_of_elements();
            MKL_INT incx = 1;
            ind = (size_t)(isamax(&n, x.begin(), &incx));
            r = x.at(ind);
        }
        catch(...)
        {
            return false;
        }

        return true;
    }

    bool addEpsilon(hoNDArray<float>& x)
    {
        try
        {
            size_t n = x.get_number_of_elements();
            float* pX = x.begin();

            long long i;

            #pragma omp parallel for default(none) private(i) shared(n, pX)
            for (i=0; i<(long long)n; i++ )
            {
                if ( GT_ABS(pX[i]) < FLT_EPSILON )
                {
                    pX[i] += GT_SGN(pX[i])*FLT_EPSILON;
                }
            }
        }
        catch(...)
        {
            return false;
        }

        return true;
    }

    bool norm2(const hoNDArray<float>& x, float& r)
    {
        try
        {
            MKL_INT incx = 1;
            MKL_INT n = x.get_number_of_elements();
            r = snrm2(&n, x.begin(), &incx);
        }
        catch(...)
        {
            return false;
        }

        return true;
    }

    bool norm1(const hoNDArray<float>& x, float& r)
    {
        try
        {
            MKL_INT incx = 1;
            MKL_INT n = x.get_number_of_elements();
            r = sasum(&n, x.begin(), &incx);
        }
        catch(...)
        {
            return false;
        }

        return true;
    }

    bool conv2(const hoNDArray<float>& x, const hoNDArray<float>& ker, hoNDArray<float>& z)
    {
        try
        {
            if ( !z.dimensions_equal(&x) )
            {
                z = x;
            }

            size_t RO = x.get_size(0);
            size_t E1 = x.get_size(1);

            size_t kerRO = ker.get_size(0);
            size_t kerE1 = ker.get_size(1);

            size_t num = x.get_number_of_elements()/(RO*E1);

            int status;
            VSLConvTaskPtr task;

            MKL_INT kerShape[2];
            kerShape[0] = kerRO; kerShape[1] = kerE1;

            MKL_INT xshape[2];
            xshape[0] = RO; xshape[1] = E1;

            MKL_INT start[2];
            start[0] = kerRO/2;
            start[1] = kerE1/2;

            MKL_INT kerStride[2], xstride[2], zstride[2];
            kerStride[0] = 1; kerStride[1] = kerRO;
            xstride[0] = 1; xstride[1] = RO;
            zstride[0] = 1; zstride[1] = RO;

            const float* pX = x.begin();
            const float* pKer = ker.begin();
            float* pZ = z.begin();

            if ( num == 1 )
            {
                status = vslsConvNewTask(&task, VSL_CONV_MODE_AUTO, 2, kerShape, xshape, xshape);
                GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

                status = vslConvSetStart(task, start);
                GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

                status = vslsConvExec(task, pKer, kerStride, pX, xstride, pZ, zstride);
                GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

                 vslConvDeleteTask(&task);
            }
            else
            {
                status = vslsConvNewTaskX(&task, VSL_CONV_MODE_AUTO, 2, kerShape, xshape, xshape, pKer, kerStride);
                GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

                status = vslConvSetStart(task, start);
                GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

                long long n;

                #pragma omp parallel for default(none) private(n) shared(num, task, pX, RO, E1, status, xstride, pZ, zstride)
                for ( n=0; n<(long long)num; n++ )
                {
                    status = vslsConvExecX(task, pX+n*RO*E1, xstride, pZ+n*RO*E1, zstride);
                }
                GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

                vslConvDeleteTask(&task);
            }
        }
        catch(...)
        {
            GADGET_ERROR_MSG("Errors happened in conv2(const hoNDArray<float>& x, const hoNDArray<float>& ker, hoNDArray<float>& z) ... ");
            return false;
        }

        return true;
    }

    bool conv3(const hoNDArray<float>& x, const hoNDArray<float>& ker, hoNDArray<float>& z)
    {
        try
        {
            if ( !z.dimensions_equal(&x) )
            {
                z = x;
            }

            size_t RO = x.get_size(0);
            size_t E1 = x.get_size(1);
            size_t E2 = x.get_size(2);

            size_t kerRO = ker.get_size(0);
            size_t kerE1 = ker.get_size(1);
            size_t kerE2 = ker.get_size(2);

            size_t num = x.get_number_of_elements()/(RO*E1*E2);

            int status;
            VSLConvTaskPtr task;

            MKL_INT kerShape[3];
            kerShape[0] = kerRO; kerShape[1] = kerE1; kerShape[2] = kerE2;

            MKL_INT xshape[3];
            xshape[0] = RO; xshape[1] = E1; xshape[2] = E2;

            MKL_INT start[3];
            start[0] = kerRO/2;
            start[1] = kerE1/2;
            start[2] = kerE2/2;

            MKL_INT kerStride[3], xstride[3], zstride[3];
            kerStride[0] = 1; kerStride[1] = kerRO; kerStride[2] = kerRO*kerE1;
            xstride[0] = 1; xstride[1] = RO; xstride[2] = RO*E1;
            zstride[0] = 1; zstride[1] = RO; zstride[2] = RO*E1;

            const float* pX = x.begin();
            const float* pKer = ker.begin();
            float* pZ = z.begin();

            if ( num == 1 )
            {
                status = vslsConvNewTask(&task, VSL_CONV_MODE_AUTO, 3, kerShape, xshape, xshape);
                GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

                status = vslConvSetStart(task, start);
                GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

                status = vslsConvExec(task, pKer, kerStride, pX, xstride, pZ, zstride);
                GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

                 vslConvDeleteTask(&task);
            }
            else
            {
                status = vslsConvNewTaskX(&task, VSL_CONV_MODE_AUTO, 3, kerShape, xshape, xshape, pKer, kerStride);
                GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

                status = vslConvSetStart(task, start);
                GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

                long long n;

                #pragma omp parallel for default(none) private(n) shared(num, task, pX, RO, E1, E2, status, xstride, pZ, zstride)
                for ( n=0; n<(long long)num; n++ )
                {
                    status = vslsConvExecX(task, pX+n*RO*E1*E2, xstride, pZ+n*RO*E1*E2, zstride);
                }
                GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

                vslConvDeleteTask(&task);
            }
        }
        catch(...)
        {
            GADGET_ERROR_MSG("Errors happened in conv3(const hoNDArray<float>& x, const hoNDArray<float>& ker, hoNDArray<float>& z) ... ");
            return false;
        }

        return true;
    }

    bool inv(const hoNDArray<float>& x, hoNDArray<float>& r)
    {
        try
        {
            if ( !r.dimensions_equal(&x) )
            {
                r = x;
            }

            long long n = x.get_number_of_elements();
            vsInv(n, x.begin(), r.begin());
            GADGET_CHECK_RETURN_FALSE(vmlGetErrStatus()==0);
        }
        catch(...)
        {
            GADGET_ERROR_MSG("Errors happened in inv(const hoNDArray<float>& x, hoNDArray<float>& r) ... ");
            return false;
        }

        return true;
    }

    // ----------------------------------------------------------------------------------------
    // double
    // ----------------------------------------------------------------------------------------

    bool add(const hoNDArray<double>& x, const hoNDArray<double>& y, hoNDArray<double>& r)
    {
        GADGET_DEBUG_CHECK_RETURN_FALSE(x.get_number_of_elements()==y.get_number_of_elements());
        if ( r.get_number_of_elements()!=x.get_number_of_elements())
        {
            r = x;
        }

        vdAdd(x.get_number_of_elements(), x.begin(), y.begin(), r.begin());
        GADGET_CHECK_RETURN_FALSE(vmlGetErrStatus()==0);

        return true;
    }

    bool subtract(const hoNDArray<double>& x, const hoNDArray<double>& y, hoNDArray<double>& r)
    {
        GADGET_DEBUG_CHECK_RETURN_FALSE(x.get_number_of_elements()==y.get_number_of_elements());
        if ( r.get_number_of_elements()!=x.get_number_of_elements())
        {
            r = x;
        }

        vdSub(x.get_number_of_elements(), x.begin(), y.begin(), r.begin());
        GADGET_CHECK_RETURN_FALSE(vmlGetErrStatus()==0);

        return true;
    }

    bool multiply(const hoNDArray<double>& x, const hoNDArray<double>& y, hoNDArray<double>& r)
    {
        GADGET_DEBUG_CHECK_RETURN_FALSE(x.get_number_of_elements()==y.get_number_of_elements());
        if ( r.get_number_of_elements()!=x.get_number_of_elements())
        {
            r = x;
        }

        vdMul(x.get_number_of_elements(), x.begin(), y.begin(), r.begin());
        GADGET_CHECK_RETURN_FALSE(vmlGetErrStatus()==0);

        return true;
    }

    bool divide(const hoNDArray<double>& x, const hoNDArray<double>& y, hoNDArray<double>& r)
    {
        GADGET_DEBUG_CHECK_RETURN_FALSE(x.get_number_of_elements()==y.get_number_of_elements());
        if ( r.get_number_of_elements()!=x.get_number_of_elements())
        {
            r = x;
        }

        vdDiv(x.get_number_of_elements(), x.begin(), y.begin(), r.begin());
        GADGET_CHECK_RETURN_FALSE(vmlGetErrStatus()==0);

        return true;
    }

    bool absolute(const hoNDArray<double>& x, hoNDArray<double>& r)
    {
        if ( r.get_number_of_elements()!=x.get_number_of_elements())
        {
            r = x;
        }

        vdAbs(x.get_number_of_elements(), x.begin(), r.begin());
        GADGET_CHECK_RETURN_FALSE(vmlGetErrStatus()==0);

        return true;
    }

    bool argument(const hoNDArray<double>& x, hoNDArray<double>& r)
    {
        if ( r.get_number_of_elements()!=x.get_number_of_elements())
        {
            r = x;
        }

        memset(r.begin(), 0, r.get_number_of_bytes());

        return true;
    }

    bool sqrt(const hoNDArray<double>& x, hoNDArray<double>& r)
    {
        if ( r.get_number_of_elements()!=x.get_number_of_elements())
        {
            r = x;
        }

        vdSqrt(x.get_number_of_elements(), x.begin(), r.begin());
        GADGET_CHECK_RETURN_FALSE(vmlGetErrStatus()==0);

        return true;
    }

    bool minAbsolute(const hoNDArray<double>& x, double& r, size_t& ind)
    {
        try
        {
            MKL_INT n = x.get_number_of_elements();
            MKL_INT incx = 1;
            ind = (size_t)(idamin(&n, x.begin(), &incx));
            r = x.at(ind);
        }
        catch(...)
        {
            return false;
        }

        return true;
    }

    bool maxAbsolute(const hoNDArray<double>& x, double& r, size_t& ind)
    {
        try
        {
            MKL_INT n = x.get_number_of_elements();
            MKL_INT incx = 1;
            ind = (size_t)(idamax(&n, x.begin(), &incx));
            r = x.at(ind);
        }
        catch(...)
        {
            return false;
        }

        return true;
    }

    bool addEpsilon(hoNDArray<double>& x)
    {
        try
        {
            size_t n = x.get_number_of_elements();
            double* pX = x.begin();

            long long i;

            #pragma omp parallel for default(none) private(i) shared(n, pX)
            for (i=0; i<(long long)n; i++ )
            {
                if ( GT_ABS(pX[i]) < DBL_EPSILON )
                {
                    pX[i] += GT_SGN(pX[i])*DBL_EPSILON;
                }
            }
        }
        catch(...)
        {
            return false;
        }

        return true;
    }

    bool norm2(const hoNDArray<double>& x, double& r)
    {
        try
        {
            MKL_INT incx = 1;
            MKL_INT n = x.get_number_of_elements();
            r = dnrm2(&n, x.begin(), &incx);
        }
        catch(...)
        {
            return false;
        }

        return true;
    }

    bool norm1(const hoNDArray<double>& x, double& r)
    {
        try
        {
            MKL_INT incx = 1;
            MKL_INT n = x.get_number_of_elements();
            r = dasum(&n, x.begin(), &incx);
        }
        catch(...)
        {
            return false;
        }

        return true;
    }

    bool conv2(const hoNDArray<double>& x, const hoNDArray<double>& ker, hoNDArray<double>& z)
    {
        try
        {
            if ( !z.dimensions_equal(&x) )
            {
                z = x;
            }

            size_t RO = x.get_size(0);
            size_t E1 = x.get_size(1);

            size_t kerRO = ker.get_size(0);
            size_t kerE1 = ker.get_size(1);

            size_t num = x.get_number_of_elements()/(RO*E1);

            int status;
            VSLConvTaskPtr task;

            MKL_INT kerShape[2];
            kerShape[0] = kerRO; kerShape[1] = kerE1;

            MKL_INT xshape[2];
            xshape[0] = RO; xshape[1] = E1;

            MKL_INT start[2];
            start[0] = kerRO/2;
            start[1] = kerE1/2;

            MKL_INT kerStride[2], xstride[2], zstride[2];
            kerStride[0] = 1; kerStride[1] = kerRO;
            xstride[0] = 1; xstride[1] = RO;
            zstride[0] = 1; zstride[1] = RO;

            const double* pX = x.begin();
            const double* pKer = ker.begin();
            double* pZ = z.begin();

            if ( num == 1 )
            {
                status = vsldConvNewTask(&task, VSL_CONV_MODE_AUTO, 2, kerShape, xshape, xshape);
                GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

                status = vslConvSetStart(task, start);
                GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

                status = vsldConvExec(task, pKer, kerStride, pX, xstride, pZ, zstride);
                GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

                 vslConvDeleteTask(&task);
            }
            else
            {
                status = vsldConvNewTaskX(&task, VSL_CONV_MODE_AUTO, 2, kerShape, xshape, xshape, pKer, kerStride);
                GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

                status = vslConvSetStart(task, start);
                GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

                long long n;

                #pragma omp parallel for default(none) private(n) shared(num, task, pX, RO, E1, status, xstride, pZ, zstride)
                for ( n=0; n<(long long)num; n++ )
                {
                    status = vsldConvExecX(task, pX+n*RO*E1, xstride, pZ+n*RO*E1, zstride);
                }
                GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

                vslConvDeleteTask(&task);
            }
        }
        catch(...)
        {
            GADGET_ERROR_MSG("Errors happened in conv2(const hoNDArray<double>& x, const hoNDArray<double>& ker, hoNDArray<double>& z) ... ");
            return false;
        }

        return true;
    }

    bool conv3(const hoNDArray<double>& x, const hoNDArray<double>& ker, hoNDArray<double>& z)
    {
        try
        {
            if ( !z.dimensions_equal(&x) )
            {
                z = x;
            }

            size_t RO = x.get_size(0);
            size_t E1 = x.get_size(1);
            size_t E2 = x.get_size(2);

            size_t kerRO = ker.get_size(0);
            size_t kerE1 = ker.get_size(1);
            size_t kerE2 = ker.get_size(2);

            size_t num = x.get_number_of_elements()/(RO*E1*E2);

            int status;
            VSLConvTaskPtr task;

            MKL_INT kerShape[3];
            kerShape[0] = kerRO; kerShape[1] = kerE1; kerShape[2] = kerE2;

            MKL_INT xshape[3];
            xshape[0] = RO; xshape[1] = E1; xshape[2] = E2;

            MKL_INT start[3];
            start[0] = kerRO/2;
            start[1] = kerE1/2;
            start[2] = kerE2/2;

            MKL_INT kerStride[3], xstride[3], zstride[3];
            kerStride[0] = 1; kerStride[1] = kerRO; kerStride[2] = kerRO*kerE1;
            xstride[0] = 1; xstride[1] = RO; xstride[2] = RO*E1;
            zstride[0] = 1; zstride[1] = RO; zstride[2] = RO*E1;

            const double* pX = x.begin();
            const double* pKer = ker.begin();
            double* pZ = z.begin();

            if ( num == 1 )
            {
                status = vsldConvNewTask(&task, VSL_CONV_MODE_AUTO, 3, kerShape, xshape, xshape);
                GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

                status = vslConvSetStart(task, start);
                GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

                status = vsldConvExec(task, pKer, kerStride, pX, xstride, pZ, zstride);
                GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

                 vslConvDeleteTask(&task);
            }
            else
            {
                status = vsldConvNewTaskX(&task, VSL_CONV_MODE_AUTO, 3, kerShape, xshape, xshape, pKer, kerStride);
                GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

                status = vslConvSetStart(task, start);
                GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

                long long n;

                #pragma omp parallel for default(none) private(n) shared(num, task, pX, RO, E1, E2, status, xstride, pZ, zstride)
                for ( n=0; n<(long long)num; n++ )
                {
                    status = vsldConvExecX(task, pX+n*RO*E1*E2, xstride, pZ+n*RO*E1*E2, zstride);
                }
                GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

                vslConvDeleteTask(&task);
            }
        }
        catch(...)
        {
            GADGET_ERROR_MSG("Errors happened in conv3(const hoNDArray<double>& x, const hoNDArray<double>& ker, hoNDArray<double>& z) ... ");
            return false;
        }

        return true;
    }

    bool inv(const hoNDArray<double>& x, hoNDArray<double>& r)
    {
        try
        {
            if ( !r.dimensions_equal(&x) )
            {
                r = x;
            }

            long long n = x.get_number_of_elements();
            vdInv(n, x.begin(), r.begin());
            GADGET_CHECK_RETURN_FALSE(vmlGetErrStatus()==0);
        }
        catch(...)
        {
            GADGET_ERROR_MSG("Errors happened in inv(const hoNDArray<double>& x, hoNDArray<double>& r) ... ");
            return false;
        }

        return true;
    }

    // ----------------------------------------------------------------------------------------
    // GT_Complex8
    // ----------------------------------------------------------------------------------------

    bool add(const hoNDArray<GT_Complex8>& x, const hoNDArray<GT_Complex8>& y, hoNDArray<GT_Complex8>& r)
    {
        GADGET_DEBUG_CHECK_RETURN_FALSE(x.get_number_of_elements()==y.get_number_of_elements());
        if ( r.get_number_of_elements()!=x.get_number_of_elements())
        {
            r = x;
        }

        vcAdd(x.get_number_of_elements(), reinterpret_cast<const MKL_Complex8*>(x.begin()), reinterpret_cast<const MKL_Complex8*>(y.begin()), reinterpret_cast<MKL_Complex8*>(r.begin()));
        GADGET_CHECK_RETURN_FALSE(vmlGetErrStatus()==0);

        return true;
    }

    bool subtract(const hoNDArray<GT_Complex8>& x, const hoNDArray<GT_Complex8>& y, hoNDArray<GT_Complex8>& r)
    {
        GADGET_DEBUG_CHECK_RETURN_FALSE(x.get_number_of_elements()==y.get_number_of_elements());
        if ( r.get_number_of_elements()!=x.get_number_of_elements())
        {
            r = x;
        }

        vcSub(x.get_number_of_elements(), reinterpret_cast<const MKL_Complex8*>(x.begin()), reinterpret_cast<const MKL_Complex8*>(y.begin()), reinterpret_cast<MKL_Complex8*>(r.begin()));
        GADGET_CHECK_RETURN_FALSE(vmlGetErrStatus()==0);

        return true;
    }

    bool multiply(const hoNDArray<GT_Complex8>& x, const hoNDArray<GT_Complex8>& y, hoNDArray<GT_Complex8>& r)
    {
        GADGET_DEBUG_CHECK_RETURN_FALSE(x.get_number_of_elements()==y.get_number_of_elements());
        if ( r.get_number_of_elements()!=x.get_number_of_elements())
        {
            r = x;
        }

        vcMul(x.get_number_of_elements(), reinterpret_cast<const MKL_Complex8*>(x.begin()), reinterpret_cast<const MKL_Complex8*>(y.begin()), reinterpret_cast<MKL_Complex8*>(r.begin()));
        GADGET_CHECK_RETURN_FALSE(vmlGetErrStatus()==0);

        return true;
    }

    bool divide(const hoNDArray<GT_Complex8>& x, const hoNDArray<GT_Complex8>& y, hoNDArray<GT_Complex8>& r)
    {
        GADGET_DEBUG_CHECK_RETURN_FALSE(x.get_number_of_elements()==y.get_number_of_elements());
        if ( r.get_number_of_elements()!=x.get_number_of_elements())
        {
            r = x;
        }

        vcDiv(x.get_number_of_elements(), reinterpret_cast<const MKL_Complex8*>(x.begin()), reinterpret_cast<const MKL_Complex8*>(y.begin()), reinterpret_cast<MKL_Complex8*>(r.begin()));
        GADGET_CHECK_RETURN_FALSE(vmlGetErrStatus()==0);

        return true;
    }

    bool absolute(const hoNDArray<GT_Complex8>& x, hoNDArray<float>& r)
    {
        if ( r.get_number_of_elements()!=x.get_number_of_elements())
        {
            r.create(x.get_dimensions());
        }

        vcAbs(x.get_number_of_elements(), reinterpret_cast<const MKL_Complex8*>(x.begin()), r.begin());
        GADGET_CHECK_RETURN_FALSE(vmlGetErrStatus()==0);

        return true;
    }

    bool absolute(const hoNDArray<GT_Complex8>& x, hoNDArray<GT_Complex8>& r)
    {
        if ( r.get_number_of_elements()!=x.get_number_of_elements())
        {
            r.create(x.get_dimensions());
        }

        hoNDArray<float> rTmp;
        rTmp.create(x.get_dimensions());

        vcAbs(x.get_number_of_elements(), reinterpret_cast<const MKL_Complex8*>(x.begin()), rTmp.begin());
        GADGET_CHECK_RETURN_FALSE(vmlGetErrStatus()==0);

        //GADGET_CHECK_RETURN_FALSE(r.copyFrom(rTmp));
	r.copyFrom(rTmp);

        return true;
    }

    bool sqrt(const hoNDArray<GT_Complex8>& x, hoNDArray<GT_Complex8>& r)
    {
        if ( r.get_number_of_elements()!=x.get_number_of_elements())
        {
            r.create(x.get_dimensions());
        }

        vcSqrt(x.get_number_of_elements(), reinterpret_cast<const MKL_Complex8*>(x.begin()), reinterpret_cast<MKL_Complex8*>(r.begin()));
        GADGET_CHECK_RETURN_FALSE(vmlGetErrStatus()==0);

        return true;
    }

    bool minAbsolute(const hoNDArray<GT_Complex8>& x, GT_Complex8& r, size_t& ind)
    {
        try
        {
            MKL_INT n = x.get_number_of_elements();
            MKL_INT incx = 1;
            ind = (size_t)(icamin(&n, reinterpret_cast<const MKL_Complex8*>(x.begin()), &incx));
            r = x.at(ind);
        }
        catch(...)
        {
            return false;
        }

        return true;
    }

    bool maxAbsolute(const hoNDArray<GT_Complex8>& x, GT_Complex8& r, size_t& ind)
    {
        try
        {
            MKL_INT n = x.get_number_of_elements();
            MKL_INT incx = 1;
            ind = (size_t)(icamax(&n, reinterpret_cast<const MKL_Complex8*>(x.begin()), &incx));
            r = x.at(ind);
        }
        catch(...)
        {
            return false;
        }

        return true;
    }

    bool multiplyConj(const hoNDArray<GT_Complex8>& x, const hoNDArray<GT_Complex8>& y, hoNDArray<GT_Complex8>& r)
    {
        GADGET_DEBUG_CHECK_RETURN_FALSE(x.get_number_of_elements()==y.get_number_of_elements());
        if ( r.get_number_of_elements()!=x.get_number_of_elements())
        {
            r = x;
        }

        vcMulByConj(x.get_number_of_elements(), reinterpret_cast<const MKL_Complex8*>(x.begin()), reinterpret_cast<const MKL_Complex8*>(y.begin()), reinterpret_cast<MKL_Complex8*>(r.begin()));
        GADGET_CHECK_RETURN_FALSE(vmlGetErrStatus()==0);

        return true;
    }

    bool argument(const hoNDArray<GT_Complex8>& x, hoNDArray<float>& r)
    {
        if ( r.get_number_of_elements()!=x.get_number_of_elements())
        {
            r.create(x.get_dimensions());
        }

        vcArg(x.get_number_of_elements(), reinterpret_cast<const MKL_Complex8*>(x.begin()), r.begin());
        GADGET_CHECK_RETURN_FALSE(vmlGetErrStatus()==0);

        return true;
    }

    bool conjugate(const hoNDArray<GT_Complex8>& x, hoNDArray<GT_Complex8>& r)
    {
        if ( r.get_number_of_elements()!=x.get_number_of_elements())
        {
            r.create(x.get_dimensions());
        }

        vcConj(x.get_number_of_elements(), reinterpret_cast<const MKL_Complex8*>(x.begin()), reinterpret_cast<MKL_Complex8*>(r.begin()));
        GADGET_CHECK_RETURN_FALSE(vmlGetErrStatus()==0);

        return true;
    }

    bool addEpsilon(hoNDArray<GT_Complex8>& x)
    {
        try
        {
            size_t n = x.get_number_of_elements();
            GT_Complex8* pX = x.begin();

            long long i;

            #pragma omp parallel for default(none) private(i) shared(n, pX)
            for (i=0; i<(long long)n; i++ )
            {
                if ( std::abs(pX[i]) < FLT_EPSILON )
                {
                    pX[i] += FLT_EPSILON;
                }
            }
        }
        catch(...)
        {
            return false;
        }

        return true;
    }

    bool norm2(const hoNDArray<GT_Complex8>& x, float& r)
    {
        try
        {
            MKL_INT incx = 1;
            MKL_INT n = x.get_number_of_elements();
            r = scnrm2(&n, reinterpret_cast<const MKL_Complex8*>(x.begin()), &incx);
        }
        catch(...)
        {
            return false;
        }

        return true;
    }

    bool norm1(const hoNDArray<GT_Complex8>& x, float& r)
    {
        try
        {
            hoNDArray<float> a;
            GADGET_CHECK_RETURN_FALSE(absolute(x, a));
            GADGET_CHECK_RETURN_FALSE(norm1(a, r));
        }
        catch(...)
        {
            return false;
        }

        return true;
    }

    bool dotc(const hoNDArray<GT_Complex8>& x, const hoNDArray<GT_Complex8>& y, GT_Complex8& r)
    {
        try
        {
            GADGET_DEBUG_CHECK_RETURN_FALSE(x.get_number_of_elements()==y.get_number_of_elements());

            MKL_INT N = x.get_number_of_elements();
            MKL_INT incx(1), incy(1);
            cdotc(reinterpret_cast<MKL_Complex8*>(&r), &N, reinterpret_cast<const MKL_Complex8*>(x.begin()), &incx, 
                    reinterpret_cast<const MKL_Complex8*>(y.begin()), &incy);
        }
        catch(...)
        {
            return false;
        }

        return true;
    }

    bool conv2(const hoNDArray<GT_Complex8>& x, const hoNDArray<GT_Complex8>& ker, hoNDArray<GT_Complex8>& z)
    {
        try
        {
            if ( !z.dimensions_equal(&x) )
            {
                z = x;
            }

            size_t RO = x.get_size(0);
            size_t E1 = x.get_size(1);

            size_t kerRO = ker.get_size(0);
            size_t kerE1 = ker.get_size(1);

            size_t num = x.get_number_of_elements()/(RO*E1);

            int status;
            VSLConvTaskPtr task;

            MKL_INT kerShape[2];
            kerShape[0] = kerRO; kerShape[1] = kerE1;

            MKL_INT xshape[2];
            xshape[0] = RO; xshape[1] = E1;

            MKL_INT start[2];
            start[0] = kerRO/2;
            start[1] = kerE1/2;

            MKL_INT kerStride[2], xstride[2], zstride[2];
            kerStride[0] = 1; kerStride[1] = kerRO;
            xstride[0] = 1; xstride[1] = RO;
            zstride[0] = 1; zstride[1] = RO;

            const MKL_Complex8* pX = reinterpret_cast<const MKL_Complex8*>(x.begin());
            const MKL_Complex8* pKer = reinterpret_cast<const MKL_Complex8*>(ker.begin());
            MKL_Complex8* pZ = reinterpret_cast<MKL_Complex8*>(z.begin());

            if ( num == 1 )
            {
                status = vslcConvNewTask(&task, VSL_CONV_MODE_AUTO, 2, kerShape, xshape, xshape);
                GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

                status = vslConvSetStart(task, start);
                GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

                status = vslcConvExec(task, pKer, kerStride, pX, xstride, pZ, zstride);
                GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

                 vslConvDeleteTask(&task);
            }
            else
            {
                status = vslcConvNewTaskX(&task, VSL_CONV_MODE_AUTO, 2, kerShape, xshape, xshape, pKer, kerStride);
                GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

                status = vslConvSetStart(task, start);
                GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

                long long n;

                #pragma omp parallel for default(none) private(n) shared(num, task, pX, RO, E1, status, xstride, pZ, zstride)
                for ( n=0; n<(long long)num; n++ )
                {
                    status = vslcConvExecX(task, pX+n*RO*E1, xstride, pZ+n*RO*E1, zstride);
                }
                GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

                vslConvDeleteTask(&task);
            }
        }
        catch(...)
        {
            GADGET_ERROR_MSG("Errors happened in conv2(const hoNDArray<GT_Complex8>& x, const hoNDArray<GT_Complex8>& ker, hoNDArray<GT_Complex8>& z) ... ");
            return false;
        }

        return true;
    }

    bool conv3(const hoNDArray<GT_Complex8>& x, const hoNDArray<GT_Complex8>& ker, hoNDArray<GT_Complex8>& z)
    {
        try
        {
            if ( !z.dimensions_equal(&x) )
            {
                z = x;
            }

            size_t RO = x.get_size(0);
            size_t E1 = x.get_size(1);
            size_t E2 = x.get_size(2);

            size_t kerRO = ker.get_size(0);
            size_t kerE1 = ker.get_size(1);
            size_t kerE2 = ker.get_size(2);

            size_t num = x.get_number_of_elements()/(RO*E1*E2);

            int status;
            VSLConvTaskPtr task;

            MKL_INT kerShape[3];
            kerShape[0] = kerRO; kerShape[1] = kerE1; kerShape[2] = kerE2;

            MKL_INT xshape[3];
            xshape[0] = RO; xshape[1] = E1; xshape[2] = E2;

            MKL_INT start[3];
            start[0] = kerRO/2;
            start[1] = kerE1/2;
            start[2] = kerE2/2;

            MKL_INT kerStride[3], xstride[3], zstride[3];
            kerStride[0] = 1; kerStride[1] = kerRO; kerStride[2] = kerRO*kerE1;
            xstride[0] = 1; xstride[1] = RO; xstride[2] = RO*E1;
            zstride[0] = 1; zstride[1] = RO; zstride[2] = RO*E1;

            const MKL_Complex8* pX = reinterpret_cast<const MKL_Complex8*>(x.begin());
            const MKL_Complex8* pKer = reinterpret_cast<const MKL_Complex8*>(ker.begin());
            MKL_Complex8* pZ = reinterpret_cast<MKL_Complex8*>(z.begin());

            if ( num == 1 )
            {
                status = vslcConvNewTask(&task, VSL_CONV_MODE_AUTO, 3, kerShape, xshape, xshape);
                GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

                status = vslConvSetStart(task, start);
                GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

                status = vslcConvExec(task, pKer, kerStride, pX, xstride, pZ, zstride);
                GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

                 vslConvDeleteTask(&task);
            }
            else
            {
                status = vslcConvNewTaskX(&task, VSL_CONV_MODE_AUTO, 3, kerShape, xshape, xshape, pKer, kerStride);
                GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

                status = vslConvSetStart(task, start);
                GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

                long long n;

                #pragma omp parallel for default(none) private(n) shared(num, task, pX, RO, E1, E2, status, xstride, pZ, zstride)
                for ( n=0; n<(long long)num; n++ )
                {
                    status = vslcConvExecX(task, pX+n*RO*E1*E2, xstride, pZ+n*RO*E1*E2, zstride);
                }
                GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

                vslConvDeleteTask(&task);
            }
        }
        catch(...)
        {
            GADGET_ERROR_MSG("Errors happened in conv3(const hoNDArray<GT_Complex8>& x, const hoNDArray<GT_Complex8>& ker, hoNDArray<GT_Complex8>& z) ... ");
            return false;
        }

        return true;
    }

    bool inv(const hoNDArray<GT_Complex8>& x, hoNDArray<GT_Complex8>& r)
    {
        try
        {
            if ( !r.dimensions_equal(&x) )
            {
                r = x;
            }

            const GT_Complex8* pX = x.begin();
            GT_Complex8* pR = r.begin();

            GT_Complex8 v(1.0);
            long long n = x.get_number_of_elements();
            long long ii;

            #pragma omp parallel for default(none) private(ii) shared(n, pX, pR, v)
            for ( ii=0; ii<n; ii++ )
            {
                pR[ii] = v/pX[ii];
            }
        }
        catch(...)
        {
            GADGET_ERROR_MSG("Errors happened in inv(const hoNDArray<GT_Complex8>& x, hoNDArray<GT_Complex8>& r) ... ");
            return false;
        }

        return true;
    }

    // ----------------------------------------------------------------------------------------
    // GT_Complex16
    // ----------------------------------------------------------------------------------------

    bool add(const hoNDArray<GT_Complex16>& x, const hoNDArray<GT_Complex16>& y, hoNDArray<GT_Complex16>& r)
    {
        GADGET_DEBUG_CHECK_RETURN_FALSE(x.get_number_of_elements()==y.get_number_of_elements());
        if ( r.get_number_of_elements()!=x.get_number_of_elements())
        {
            r = x;
        }

        vzAdd(x.get_number_of_elements(), reinterpret_cast<const MKL_Complex16*>(x.begin()), reinterpret_cast<const MKL_Complex16*>(y.begin()), reinterpret_cast<MKL_Complex16*>(r.begin()));
        GADGET_CHECK_RETURN_FALSE(vmlGetErrStatus()==0);

        return true;
    }

    bool subtract(const hoNDArray<GT_Complex16>& x, const hoNDArray<GT_Complex16>& y, hoNDArray<GT_Complex16>& r)
    {
        GADGET_DEBUG_CHECK_RETURN_FALSE(x.get_number_of_elements()==y.get_number_of_elements());
        if ( r.get_number_of_elements()!=x.get_number_of_elements())
        {
            r = x;
        }

        vzSub(x.get_number_of_elements(), reinterpret_cast<const MKL_Complex16*>(x.begin()), reinterpret_cast<const MKL_Complex16*>(y.begin()), reinterpret_cast<MKL_Complex16*>(r.begin()));
        GADGET_CHECK_RETURN_FALSE(vmlGetErrStatus()==0);

        return true;
    }

    bool multiply(const hoNDArray<GT_Complex16>& x, const hoNDArray<GT_Complex16>& y, hoNDArray<GT_Complex16>& r)
    {
        GADGET_DEBUG_CHECK_RETURN_FALSE(x.get_number_of_elements()==y.get_number_of_elements());
        if ( r.get_number_of_elements()!=x.get_number_of_elements())
        {
            r = x;
        }

        vzMul(x.get_number_of_elements(), reinterpret_cast<const MKL_Complex16*>(x.begin()), reinterpret_cast<const MKL_Complex16*>(y.begin()), reinterpret_cast<MKL_Complex16*>(r.begin()));
        GADGET_CHECK_RETURN_FALSE(vmlGetErrStatus()==0);

        return true;
    }

    bool divide(const hoNDArray<GT_Complex16>& x, const hoNDArray<GT_Complex16>& y, hoNDArray<GT_Complex16>& r)
    {
        GADGET_DEBUG_CHECK_RETURN_FALSE(x.get_number_of_elements()==y.get_number_of_elements());
        if ( r.get_number_of_elements()!=x.get_number_of_elements())
        {
            r = x;
        }

        vzDiv(x.get_number_of_elements(), reinterpret_cast<const MKL_Complex16*>(x.begin()), reinterpret_cast<const MKL_Complex16*>(y.begin()), reinterpret_cast<MKL_Complex16*>(r.begin()));
        GADGET_CHECK_RETURN_FALSE(vmlGetErrStatus()==0);

        return true;
    }

    bool absolute(const hoNDArray<GT_Complex16>& x, hoNDArray<double>& r)
    {
        if ( r.get_number_of_elements()!=x.get_number_of_elements())
        {
            r.create(x.get_dimensions());
        }

        vzAbs(x.get_number_of_elements(), reinterpret_cast<const MKL_Complex16*>(x.begin()), r.begin());
        GADGET_CHECK_RETURN_FALSE(vmlGetErrStatus()==0);

        return true;
    }

    bool absolute(const hoNDArray<GT_Complex16>& x, hoNDArray<GT_Complex16>& r)
    {
        if ( r.get_number_of_elements()!=x.get_number_of_elements())
        {
            r.create(x.get_dimensions());
        }

        hoNDArray<double> rTmp;
        rTmp.create(x.get_dimensions());

        vzAbs(x.get_number_of_elements(), reinterpret_cast<const MKL_Complex16*>(x.begin()), rTmp.begin());
        GADGET_CHECK_RETURN_FALSE(vmlGetErrStatus()==0);

        //GADGET_CHECK_RETURN_FALSE(r.copyFrom(rTmp));
	r.copyFrom(rTmp);

        return true;
    }

    bool sqrt(const hoNDArray<GT_Complex16>& x, hoNDArray<GT_Complex16>& r)
    {
        if ( r.get_number_of_elements()!=x.get_number_of_elements())
        {
            r.create(x.get_dimensions());
        }

        vzSqrt(x.get_number_of_elements(), reinterpret_cast<const MKL_Complex16*>(x.begin()), reinterpret_cast<MKL_Complex16*>(r.begin()));
        GADGET_CHECK_RETURN_FALSE(vmlGetErrStatus()==0);

        return true;
    }

    bool minAbsolute(const hoNDArray<GT_Complex16>& x, GT_Complex16& r, size_t& ind)
    {
        try
        {
            MKL_INT n = x.get_number_of_elements();
            MKL_INT incx = 1;
            ind = (size_t)(izamin(&n, reinterpret_cast<const MKL_Complex16*>(x.begin()), &incx));
            r = x.at(ind);
        }
        catch(...)
        {
            return false;
        }

        return true;
    }

    bool maxAbsolute(const hoNDArray<GT_Complex16>& x, GT_Complex16& r, size_t& ind)
    {
        try
        {
            MKL_INT n = x.get_number_of_elements();
            MKL_INT incx = 1;
            ind = (size_t)(izamax(&n, reinterpret_cast<const MKL_Complex16*>(x.begin()), &incx));
            r = x.at(ind);
        }
        catch(...)
        {
            return false;
        }

        return true;
    }

    bool multiplyConj(const hoNDArray<GT_Complex16>& x, const hoNDArray<GT_Complex16>& y, hoNDArray<GT_Complex16>& r)
    {
        GADGET_DEBUG_CHECK_RETURN_FALSE(x.get_number_of_elements()==y.get_number_of_elements());
        if ( r.get_number_of_elements()!=x.get_number_of_elements())
        {
            r = x;
        }

        vzMulByConj(x.get_number_of_elements(), reinterpret_cast<const MKL_Complex16*>(x.begin()), reinterpret_cast<const MKL_Complex16*>(y.begin()), reinterpret_cast<MKL_Complex16*>(r.begin()));
        GADGET_CHECK_RETURN_FALSE(vmlGetErrStatus()==0);

        return true;
    }

    bool argument(const hoNDArray<GT_Complex16>& x, hoNDArray<double>& r)
    {
        if ( r.get_number_of_elements()!=x.get_number_of_elements())
        {
            r.create(x.get_dimensions());
        }

        vzArg(x.get_number_of_elements(), reinterpret_cast<const MKL_Complex16*>(x.begin()), r.begin());
        GADGET_CHECK_RETURN_FALSE(vmlGetErrStatus()==0);

        return true;
    }

    bool conjugate(const hoNDArray<GT_Complex16>& x, hoNDArray<GT_Complex16>& r)
    {
        if ( r.get_number_of_elements()!=x.get_number_of_elements())
        {
            r.create(x.get_dimensions());
        }

        vzConj(x.get_number_of_elements(), reinterpret_cast<const MKL_Complex16*>(x.begin()), reinterpret_cast<MKL_Complex16*>(r.begin()));
        GADGET_CHECK_RETURN_FALSE(vmlGetErrStatus()==0);

        return true;
    }

    bool addEpsilon(hoNDArray<GT_Complex16>& x)
    {
        try
        {
            size_t n = x.get_number_of_elements();
            GT_Complex16* pX = x.begin();

            long long i;

            #pragma omp parallel for default(none) private(i) shared(n, pX)
            for (i=0; i<(long long)n; i++ )
            {
                if ( std::abs(pX[i]) < DBL_EPSILON )
                {
                    pX[i] += DBL_EPSILON;
                }
            }
        }
        catch(...)
        {
            return false;
        }

        return true;
    }

    bool norm2(const hoNDArray<GT_Complex16>& x, double& r)
    {
        try
        {
            MKL_INT incx = 1;
            MKL_INT n = x.get_number_of_elements();
            r = dznrm2(&n, reinterpret_cast<const MKL_Complex16*>(x.begin()), &incx);
        }
        catch(...)
        {
            return false;
        }

        return true;
    }

    bool norm1(const hoNDArray<GT_Complex16>& x, double& r)
    {
        try
        {
            hoNDArray<double> a;
            GADGET_CHECK_RETURN_FALSE(absolute(x, a));
            GADGET_CHECK_RETURN_FALSE(norm1(a, r));
        }
        catch(...)
        {
            return false;
        }

        return true;
    }

    bool dotc(const hoNDArray<GT_Complex16>& x, const hoNDArray<GT_Complex16>& y, GT_Complex16& r)
    {
        try
        {
            GADGET_DEBUG_CHECK_RETURN_FALSE(x.get_number_of_elements()==y.get_number_of_elements());

            MKL_INT N = x.get_number_of_elements();
            MKL_INT incx(1), incy(1);
            zdotc(reinterpret_cast<MKL_Complex16*>(&r), &N, reinterpret_cast<const MKL_Complex16*>(x.begin()), &incx, 
                    reinterpret_cast<const MKL_Complex16*>(y.begin()), &incy);
        }
        catch(...)
        {
            return false;
        }

        return true;
    }

    bool conv2(const hoNDArray<GT_Complex16>& x, const hoNDArray<GT_Complex16>& ker, hoNDArray<GT_Complex16>& z)
    {
        try
        {
            if ( !z.dimensions_equal(&x) )
            {
                z = x;
            }

            size_t RO = x.get_size(0);
            size_t E1 = x.get_size(1);

            size_t kerRO = ker.get_size(0);
            size_t kerE1 = ker.get_size(1);

            size_t num = x.get_number_of_elements()/(RO*E1);

            int status;
            VSLConvTaskPtr task;

            MKL_INT kerShape[2];
            kerShape[0] = kerRO; kerShape[1] = kerE1;

            MKL_INT xshape[2];
            xshape[0] = RO; xshape[1] = E1;

            MKL_INT start[2];
            start[0] = kerRO/2;
            start[1] = kerE1/2;

            MKL_INT kerStride[2], xstride[2], zstride[2];
            kerStride[0] = 1; kerStride[1] = kerRO;
            xstride[0] = 1; xstride[1] = RO;
            zstride[0] = 1; zstride[1] = RO;

            const MKL_Complex16* pX = reinterpret_cast<const MKL_Complex16*>(x.begin());
            const MKL_Complex16* pKer = reinterpret_cast<const MKL_Complex16*>(ker.begin());
            MKL_Complex16* pZ = reinterpret_cast<MKL_Complex16*>(z.begin());

            if ( num == 1 )
            {
                status = vslzConvNewTask(&task, VSL_CONV_MODE_AUTO, 2, kerShape, xshape, xshape);
                GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

                status = vslConvSetStart(task, start);
                GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

                status = vslzConvExec(task, pKer, kerStride, pX, xstride, pZ, zstride);
                GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

                 vslConvDeleteTask(&task);
            }
            else
            {
                status = vslzConvNewTaskX(&task, VSL_CONV_MODE_AUTO, 2, kerShape, xshape, xshape, pKer, kerStride);
                GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

                status = vslConvSetStart(task, start);
                GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

                long long n;

                #pragma omp parallel for default(none) private(n) shared(num, task, pX, RO, E1, status, xstride, pZ, zstride)
                for ( n=0; n<(long long)num; n++ )
                {
                    status = vslzConvExecX(task, pX+n*RO*E1, xstride, pZ+n*RO*E1, zstride);
                }
                GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

                vslConvDeleteTask(&task);
            }
        }
        catch(...)
        {
            GADGET_ERROR_MSG("Errors happened in conv2(const hoNDArray<GT_Complex16>& x, const hoNDArray<GT_Complex16>& ker, hoNDArray<GT_Complex16>& z) ... ");
            return false;
        }

        return true;
    }

    bool conv3(const hoNDArray<GT_Complex16>& x, const hoNDArray<GT_Complex16>& ker, hoNDArray<GT_Complex16>& z)
    {
        try
        {
            if ( !z.dimensions_equal(&x) )
            {
                z = x;
            }

            size_t RO = x.get_size(0);
            size_t E1 = x.get_size(1);
            size_t E2 = x.get_size(2);

            size_t kerRO = ker.get_size(0);
            size_t kerE1 = ker.get_size(1);
            size_t kerE2 = ker.get_size(2);

            size_t num = x.get_number_of_elements()/(RO*E1*E2);

            int status;
            VSLConvTaskPtr task;

            MKL_INT kerShape[3];
            kerShape[0] = kerRO; kerShape[1] = kerE1; kerShape[2] = kerE2;

            MKL_INT xshape[3];
            xshape[0] = RO; xshape[1] = E1; xshape[2] = E2;

            MKL_INT start[3];
            start[0] = kerRO/2;
            start[1] = kerE1/2;
            start[2] = kerE2/2;

            MKL_INT kerStride[3], xstride[3], zstride[3];
            kerStride[0] = 1; kerStride[1] = kerRO; kerStride[2] = kerRO*kerE1;
            xstride[0] = 1; xstride[1] = RO; xstride[2] = RO*E1;
            zstride[0] = 1; zstride[1] = RO; zstride[2] = RO*E1;

            const MKL_Complex16* pX = reinterpret_cast<const MKL_Complex16*>(x.begin());
            const MKL_Complex16* pKer = reinterpret_cast<const MKL_Complex16*>(ker.begin());
            MKL_Complex16* pZ = reinterpret_cast<MKL_Complex16*>(z.begin());

            if ( num == 1 )
            {
                status = vslzConvNewTask(&task, VSL_CONV_MODE_AUTO, 3, kerShape, xshape, xshape);
                GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

                status = vslConvSetStart(task, start);
                GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

                status = vslzConvExec(task, pKer, kerStride, pX, xstride, pZ, zstride);
                GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

                 vslConvDeleteTask(&task);
            }
            else
            {
                status = vslzConvNewTaskX(&task, VSL_CONV_MODE_AUTO, 3, kerShape, xshape, xshape, pKer, kerStride);
                GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

                status = vslConvSetStart(task, start);
                GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

                long long n;

                #pragma omp parallel for default(none) private(n) shared(num, task, pX, RO, E1, E2, status, xstride, pZ, zstride)
                for ( n=0; n<(long long)num; n++ )
                {
                    status = vslzConvExecX(task, pX+n*RO*E1*E2, xstride, pZ+n*RO*E1*E2, zstride);
                }
                GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

                vslConvDeleteTask(&task);
            }
        }
        catch(...)
        {
            GADGET_ERROR_MSG("Errors happened in conv3(const hoNDArray<GT_Complex16>& x, const hoNDArray<GT_Complex16>& ker, hoNDArray<GT_Complex16>& z) ... ");
            return false;
        }

        return true;
    }

    bool inv(const hoNDArray<GT_Complex16>& x, hoNDArray<GT_Complex16>& r)
    {
        try
        {
            if ( !r.dimensions_equal(&x) )
            {
                r = x;
            }

            const GT_Complex16* pX = x.begin();
            GT_Complex16* pR = r.begin();

            GT_Complex16 v(1.0);
            long long n = x.get_number_of_elements();
            long long ii;

            #pragma omp parallel for default(none) private(ii) shared(n, pX, pR, v)
            for ( ii=0; ii<n; ii++ )
            {
                pR[ii] = v/pX[ii];
            }
        }
        catch(...)
        {
            GADGET_ERROR_MSG("Errors happened in inv(const hoNDArray<GT_Complex16>& x, hoNDArray<GT_Complex16>& r) ... ");
            return false;
        }

        return true;
    }

    // ----------------------------------------------------------------------------------------
    // templated functions
    // ----------------------------------------------------------------------------------------

    template<typename T> 
    bool sumOverLastDimension(const hoNDArray<T>& x, hoNDArray<T>& r)
    {
        try
        {
            boost::shared_ptr< std::vector<size_t> > dim = x.get_dimensions();
            size_t NDim = dim->size();

            std::vector<size_t> dimR(NDim-1);

            size_t d;
            for ( d=0; d<NDim-1; d++ )
            {
                dimR[d] = (*dim)[d];
            }

            if ( !r.dimensions_equal(&dimR) )
            {
                r.create(&dimR);
            }

            // Gadgetron::clear(&r);

            if ( x.get_size(NDim-1) <= 1 )
            {
                memcpy(r.begin(), x.begin(), x.get_number_of_bytes());
                return true;
            }

            size_t lastDim = x.get_size(NDim-1);
            size_t NR = r.get_number_of_elements();
            T* pA = const_cast<T*>(x.begin());
            T* pR = r.begin();

            memcpy(pR, pA, sizeof(T)*NR);

            // sum over the last dim
            hoNDArray<T> tmp;
            for ( d=1; d<lastDim; d++ )
            {
                tmp.create(&dimR, pA+d*NR);
                add(tmp, r, r);
            }
        }
        catch (...)
        {
            GADGET_ERROR_MSG("Errors in sumOverLastDimension(const hoNDArray<T>& x, hoNDArray<T>& r) ... ");
            return false;
        }
        return true;
    }

    template<typename T> 
    bool sumOverSecondLastDimension(const hoNDArray<T>& x, hoNDArray<T>& r)
    {
        try
        {
            boost::shared_ptr< std::vector<size_t> > dim = x.get_dimensions();
            size_t NDim = dim->size();

            if ( NDim < 2 ) return true;

            std::vector<size_t> dimR(NDim-1);
            std::vector<size_t> dimRInternal(NDim-2);

            size_t d;
            for ( d=0; d<NDim-2; d++ )
            {
                dimR[d] = (*dim)[d];
                dimRInternal[d] = (*dim)[d];
            }
            dimR[NDim-2] = (*dim)[NDim-1];

            if ( !r.dimensions_equal(&dimR) )
            {
                r.create(&dimR);
            }

            if ( x.get_size(NDim-2) <= 1 )
            {
                memcpy(r.begin(), x.begin(), x.get_number_of_bytes());
                return true;
            }

            size_t lastDim = x.get_size(NDim-1);
            size_t secondLastDim = x.get_size(NDim-2);
            size_t NS = x.get_number_of_elements()/lastDim;
            size_t NR = r.get_number_of_elements()/lastDim;
            T* pA = const_cast<T*>(x.begin());
            T* pR = r.begin();

            //int l;
            //#pragma omp parallel default(none) private(l) shared(lastDim, secondLastDim, NR, pA, pR, dimRInternal)
            //{
            //    hoNDArray<T> tmp(&dimRInternal);

            //    #pragma omp for
            //    for ( l=0; l<(int)lastDim; l++ )
            //    {
            //        memcpy(tmp.begin(), pA+l*NR*secondLastDim, sizeof(T)*NR);
            //        for ( size_t s=1; s<secondLastDim; s++ )
            //        {
            //            hoNDArray<T> tmp2;
            //            tmp2.create(&dimRInternal, pA+l*NR*secondLastDim+s*NR);
            //            add(tmp, tmp2, tmp);
            //        }

            //        memcpy(pR+l*NR, tmp.begin(), sizeof(T)*NR);
            //    }
            //}

            int l;
            #pragma omp parallel default(none) private(l) shared(lastDim, secondLastDim, NS, NR, pA, pR, dimRInternal)
            {
                hoNDArray<T> tmp, tmp2;

                #pragma omp for
                for ( l=0; l<(int)lastDim; l++ )
                {
                    memcpy(pR+l*NR, pA+l*NS, sizeof(T)*NR);
                    tmp.create(&dimRInternal, pR+l*NR);
                    for ( size_t s=1; s<secondLastDim; s++ )
                    {
                        tmp2.create(&dimRInternal, pA+l*NS+s*NR);
                        add(tmp, tmp2, tmp);
                    }
                }
            }
        }
        catch (...)
        {
            GADGET_ERROR_MSG("Errors in sumOverSecondLastDimension(const hoNDArray<T>& x, hoNDArray<T>& r) ... ");
            return false;
        }
        return true;
    }

    // e.g. x is 3D and y is 4D array, r(:,:,:,n) = y(:,:,:,n) .* x
    template<typename T> 
    bool multiplyOverLastDimension(const hoNDArray<T>& x, const hoNDArray<T>& y, hoNDArray<T>& r)
    {
        try
        {
            boost::shared_ptr< std::vector<size_t> > dimX = x.get_dimensions();
            boost::shared_ptr< std::vector<size_t> > dimY = y.get_dimensions();

            size_t NDim = dimY->size();

            GADGET_CHECK_RETURN_FALSE(dimX->size()==NDim-1);

            if ( !r.dimensions_equal(dimY.get()) )
            {
                r.create(dimY);
            }

            if ( y.get_size(NDim-1) <= 1 )
            {
                GADGET_CHECK_RETURN_FALSE(multiply(x, y, r));
                return true;
            }

            size_t lastDim = y.get_size(NDim-1);
            size_t N = x.get_number_of_elements();
            const T* pX = x.begin();
            const T* pY = y.begin();
            T* pR = r.begin();

            int d;

            #ifdef GCC_OLD_FLAG
                #pragma omp parallel default(none) private(d) shared(dimX, lastDim, N, pY, pR)
            #else
                #pragma omp parallel default(none) private(d) shared(x, dimX, lastDim, N, pY, pR)
            #endif
            {
                hoNDArray<T> tmpY, tmpR;

                #pragma omp for
                for ( d=0; d<(int)lastDim; d++ )
                {
                    tmpY.create(dimX.get(), const_cast<T*>(pY+d*N));
                    tmpR.create(dimX.get(), pR+d*N);
                    multiply(x, tmpY, tmpR);
                }
            }
        }
        catch (...)
        {
            GADGET_ERROR_MSG("Errors in multiplyOverLastDimension(const hoNDArray<float>& x, const hoNDArray<float>& y, hoNDArray<float>& r) ... ");
            return false;
        }
        return true;
    }

    // e.g. x is 3D and y is 4D array, r(:,:,:,n) = y(:,:,:,n) ./ x
    template<typename T> 
    bool divideOverLastDimension(const hoNDArray<T>& x, const hoNDArray<T>& y, hoNDArray<T>& r)
    {
        try
        {
            boost::shared_ptr< std::vector<size_t> > dimX = x.get_dimensions();
            boost::shared_ptr< std::vector<size_t> > dimY = y.get_dimensions();

            size_t NDim = dimY->size();

            GADGET_CHECK_RETURN_FALSE(dimX->size()==NDim-1);

            if ( !r.dimensions_equal(dimY.get()) )
            {
                r.create(dimY);
            }

            if ( y.get_size(NDim-1) <= 1 )
            {
                GADGET_CHECK_RETURN_FALSE(divide(y, x, r));
                return true;
            }

            size_t lastDim = y.get_size(NDim-1);
            size_t N = x.get_number_of_elements();
            T* pY = const_cast<T*>(y.begin());
            T* pR = r.begin();

            int d;

            #ifdef GCC_OLD_FLAG
                #pragma omp parallel default(none) private(d) shared(dimX, lastDim, N, pY, pR)
            #else
                #pragma omp parallel default(none) private(d) shared(x, dimX, lastDim, N, pY, pR)
            #endif
            {
                hoNDArray<T> tmpY, tmpR;

                #pragma omp for
                for ( d=0; d<(int)lastDim; d++ )
                {
                    tmpY.create(dimX, pY+d*N);
                    tmpR.create(dimX, pR+d*N);
                    divide(tmpY, x, tmpR);
                }
            }
        }
        catch (...)
        {
            GADGET_ERROR_MSG("Errors in divideOverLastDimension(const hoNDArray<T>& x, const hoNDArray<T>& y, hoNDArray<T>& r) ... ");
            return false;
        }
        return true;
    }

    template<typename T> 
    bool sumOver1stDimension(const hoNDArray<T>& x, hoNDArray<T>& r)
    {
        try
        {
            size_t RO = x.get_size(0);
            size_t num = x.get_number_of_elements()/(RO);

            boost::shared_ptr< std::vector<size_t> > dim = x.get_dimensions();

            std::vector<size_t> dimAve(*dim);
            dimAve[0] = 1;
            r.create(&dimAve);

            const T* pX = x.begin();
            T* pR = r.begin();

            int n;
            #pragma omp parallel for default(none) private(n) shared(RO, num, pX, pR)
            for ( n=0; n<(int)num; n++ )
            {
                T xsum = pX[n*RO];
                for (size_t ro=1; ro<RO; ro++ )
                {
                    xsum += pX[n*RO+ro];
                }

                pR[n] = xsum;
            }
        }
        catch(...)
        {
            GADGET_ERROR_MSG("Errors in sumOver1stDimension(...) ... ");
            return false;
        }

        return true;
    }

    template<typename T> 
    bool sumOver2ndDimension(const hoNDArray<T>& x, hoNDArray<T>& r)
    {
        try
        {
            size_t NDim = x.get_number_of_dimensions();

            if ( NDim < 2 )
            {
                r = x;
                return true;
            }

            size_t RO = x.get_size(0);
            size_t E1 = x.get_size(1);

            size_t num = x.get_number_of_elements()/(RO*E1);

            boost::shared_ptr< std::vector<size_t> > dim = x.get_dimensions();

            std::vector<size_t> dimAve(*dim);
            dimAve[1] = 1;
            r.create(&dimAve);

            int n;
            #ifdef GCC_OLD_FLAG
                #pragma omp parallel for default(none) private(n) shared(RO, E1, num)
            #else
                #pragma omp parallel for default(none) private(n) shared(RO, E1, num, x, r)
            #endif
            for ( n=0; n<(int)num; n++ )
            {
                hoNDArray<T> xsum(RO, const_cast<T*>(r.begin()+n*RO));
                memcpy(xsum.begin(), x.begin()+n*RO*E1, xsum.get_number_of_bytes());

                for (size_t e1=1; e1<E1; e1++ )
                {
                    hoNDArray<T> x1D(RO, const_cast<T*>(x.begin()+n*RO*E1+e1*RO));
                    Gadgetron::add(x1D, xsum, xsum);
                }
            }
        }
        catch(...)
        {
            GADGET_ERROR_MSG("Errors in sumOver2ndDimension(...) ... ");
            return false;
        }

        return true;
    }

    template<typename T> 
    bool sumOver3rdDimension(const hoNDArray<T>& x, hoNDArray<T>& r)
    {
        try
        {
            size_t NDim = x.get_number_of_dimensions();

            if ( NDim < 3 )
            {
                r = x;
                return true;
            }

            size_t RO = x.get_size(0);
            size_t E1 = x.get_size(1);
            size_t CHA = x.get_size(2);

            size_t num = x.get_number_of_elements()/(RO*E1*CHA);

            boost::shared_ptr< std::vector<size_t> > dim = x.get_dimensions();

            std::vector<size_t> dimAve(*dim);
            dimAve[2] = 1;
            r.create(&dimAve);

            int n;
            #ifdef GCC_OLD_FLAG
                #pragma omp parallel for default(none) private(n) shared(RO, E1, CHA, num)
            #else
                #pragma omp parallel for default(none) private(n) shared(RO, E1, CHA, num, x, r)
            #endif 
            for ( n=0; n<(int)num; n++ )
            {
                hoNDArray<T> xsum(RO, E1, const_cast<T*>(r.begin()+n*RO*E1));
                memcpy(xsum.begin(), x.begin()+n*RO*E1*CHA, xsum.get_number_of_bytes());

                for (size_t cha=1; cha<CHA; cha++ )
                {
                    hoNDArray<T> x2D(RO, E1, const_cast<T*>(x.begin()+n*RO*E1*CHA+cha*RO*E1));
                    Gadgetron::add(x2D, xsum, xsum);
                }
            }
        }
        catch(...)
        {
            GADGET_ERROR_MSG("Errors in sumOver3rdDimension(...) ... ");
            return false;
        }

        return true;
    }

    template<typename T> bool sumOver4thDimension(const hoNDArray<T>& x, hoNDArray<T>& r)
    {
        try
        {
            size_t NDim = x.get_number_of_dimensions();

            if ( NDim < 4 )
            {
                r = x;
                return true;
            }

            size_t RO = x.get_size(0);
            size_t E1 = x.get_size(1);
            size_t CHA = x.get_size(2);
            size_t N = x.get_size(3);

            size_t num = x.get_number_of_elements()/(RO*E1*CHA*N);

            boost::shared_ptr< std::vector<size_t> > dim = x.get_dimensions();

            std::vector<size_t> dimAve(*dim);
            dimAve[3] = 1;
            r.create(&dimAve);

            int n;
            #ifdef GCC_OLD_FLAG
                #pragma omp parallel for default(none) private(n) shared(RO, E1, CHA, N, num)
            #else
                #pragma omp parallel for default(none) private(n) shared(RO, E1, CHA, N, num, x, r)
            #endif
            for ( n=0; n<(int)num; n++ )
            {
                hoNDArray<T> xsum(RO, E1, CHA, const_cast<T*>(r.begin()+n*RO*E1*CHA));
                memcpy(xsum.begin(), x.begin()+n*RO*E1*CHA*N, xsum.get_number_of_bytes());

                for (size_t nn=1; nn<N; nn++ )
                {
                    hoNDArray<T> x3D(RO, E1, CHA, const_cast<T*>(x.begin()+n*RO*E1*CHA*N+nn*RO*E1*CHA));
                    Gadgetron::add(x3D, xsum, xsum);
                }
            }
        }
        catch(...)
        {
            GADGET_ERROR_MSG("Errors in sumOver4thDimension(const hoNDArray<T>& x, hoNDArray<T>& r) ... ");
            return false;
        }

        return true;
    }

    template<typename T> bool sumOver5thDimension(const hoNDArray<T>& x, hoNDArray<T>& r)
    {
        try
        {
            size_t NDim = x.get_number_of_dimensions();

            if ( NDim < 5 )
            {
                r = x;
                return true;
            }

            size_t RO = x.get_size(0);
            size_t E1 = x.get_size(1);
            size_t CHA = x.get_size(2);
            size_t N = x.get_size(3);
            size_t S = x.get_size(4);

            size_t num = x.get_number_of_elements()/(RO*E1*CHA*N*S);

            boost::shared_ptr< std::vector<size_t> > dim = x.get_dimensions();

            std::vector<size_t> dimAve(*dim);
            dimAve[4] = 1;
            r.create(&dimAve);

            int n;
            #ifdef GCC_OLD_FLAG
                #pragma omp parallel for default(none) private(n) shared(RO, E1, CHA, N, S, num) if (num > 4)
            #else
                #pragma omp parallel for default(none) private(n) shared(RO, E1, CHA, N, S, num, x, r) if (num > 4)
            #endif
            for ( n=0; n<(int)num; n++ )
            {
                hoNDArray<T> xsum(RO, E1, CHA, N, const_cast<T*>(r.begin()+n*RO*E1*CHA*N));
                memcpy(xsum.begin(), x.begin()+n*RO*E1*CHA*N*S, xsum.get_number_of_bytes());

                for (size_t s=1; s<S; s++ )
                {
                    hoNDArray<T> x4D(RO, E1, CHA, N, const_cast<T*>(x.begin()+n*RO*E1*CHA*N*S+s*RO*E1*CHA*N));
                    Gadgetron::add(x4D, xsum, xsum);
                }
            }
        }
        catch(...)
        {
            GADGET_ERROR_MSG("Errors in sumOver5thDimension(const hoNDArray<T>& x, hoNDArray<T>& r) ... ");
            return false;
        }

        return true;
    }

    // e.g. x is 3D and y is 4D array, r(:,:,n,:) = y(:,:,n,:) .* x3D
    template<typename T> 
    bool multiplyOver3rdDimension(const hoNDArray<T>& x3D, const hoNDArray<T>& y4D, hoNDArray<T>& r)
    {
        try
        {
            boost::shared_ptr< std::vector<size_t> > dimX = x3D.get_dimensions();
            boost::shared_ptr< std::vector<size_t> > dimY = y4D.get_dimensions();

            size_t NDim = dimY->size();

            GADGET_CHECK_RETURN_FALSE(dimX->size()>=3);
            GADGET_CHECK_RETURN_FALSE(NDim>=4);
            GADGET_CHECK_RETURN_FALSE((*dimX)[0]==(*dimY)[0]);
            GADGET_CHECK_RETURN_FALSE((*dimX)[1]==(*dimY)[1]);

            if ( !r.dimensions_equal(dimY.get()) )
            {
                r.create(dimY);
            }

            int t, N2D = x3D.get_size(0)*x3D.get_size(1);
            int sz = y4D.get_size(2);
            int st = y4D.get_number_of_elements()/(N2D*sz);

            if ( sz == 1 )
            {
                GADGET_CHECK_RETURN_FALSE(multiply(x3D, y4D, r));
                return true;
            }

            const T* pX = x3D.begin();
            const T* pY = y4D.begin();
            T* pR = r.begin();

            std::vector<size_t> dim2D(2);
            dim2D[0] = (*dimY)[0];
            dim2D[1] = (*dimY)[1];

            #pragma omp parallel for default(none) private(t) shared(N2D, sz, st, dim2D, pX, pY, pR)
            for ( t=0; t<st; t++ )
            {
                hoNDArray<T> tmpX, tmpY, tmpR;
                tmpX.create(&dim2D, const_cast<T*>(pX+t*N2D));

                for ( int z=0; z<sz; z++ )
                {
                    tmpY.create(&dim2D, const_cast<T*>(pY+t*N2D*sz+z*N2D));
                    tmpR.create(&dim2D, pR+t*N2D*sz+z*N2D);
                    multiply(tmpX, tmpY, tmpR);
                }
            }
        }
        catch (...)
        {
            GADGET_ERROR_MSG("Errors in multiplyOver3rdDimension(const hoNDArray<float>& x3D, const hoNDArray<float>& y4D, hoNDArray<float>& r) ... ");
            return false;
        }
        return true;
    }

    template<typename T> 
    bool multiplyOver4thDimension(const hoNDArray<T>& x4D, const hoNDArray<T>& y5D, hoNDArray<T>& r)
    {
        try
        {
            boost::shared_ptr< std::vector<size_t> > dimX = x4D.get_dimensions();
            boost::shared_ptr< std::vector<size_t> > dimY = y5D.get_dimensions();

            size_t NDim = dimY->size();

            GADGET_CHECK_RETURN_FALSE(dimX->size()>=4);
            GADGET_CHECK_RETURN_FALSE((*dimX)[0]==(*dimY)[0]);
            GADGET_CHECK_RETURN_FALSE((*dimX)[1]==(*dimY)[1]);
            GADGET_CHECK_RETURN_FALSE((*dimX)[2]==(*dimY)[2]);

            if ( !r.dimensions_equal(dimY.get()) )
            {
                r.create(dimY);
            }

            size_t RO = (*dimX)[0];
            size_t E1 = (*dimX)[1];
            size_t CHA = (*dimX)[2];

            int t, N3D = RO*E1*CHA;

            size_t N = (*dimY)[3];
            size_t num = x4D.get_number_of_elements()/(RO*E1*CHA);

            const T* pX = x4D.begin();
            const T* pY = y5D.begin();
            T* pR = r.begin();

            std::vector<size_t> dim3D(3);
            dim3D[0] = RO;
            dim3D[1] = E1;
            dim3D[2] = CHA;

            #pragma omp parallel for default(none) private(t) shared(N3D, N, dim3D, pX, pY, pR, num)
            for ( t=0; t<(int)num; t++ )
            {
                hoNDArray<T> tmpX, tmpY, tmpR;
                tmpX.create(&dim3D, const_cast<T*>(pX+t*N3D));

                for ( int n=0; n<N; n++ )
                {
                    tmpY.create(&dim3D, const_cast<T*>(pY+t*N3D*N+n*N3D));
                    tmpR.create(&dim3D, pR+t*N3D*N+n*N3D);
                    multiply(tmpX, tmpY, tmpR);
                }
            }
        }
        catch (...)
        {
            GADGET_ERROR_MSG("Errors in multiplyOver4thDimension(const hoNDArray<float>& x4D, const hoNDArray<float>& y5D, hoNDArray<float>& r) ... ");
            return false;
        }
        return true;
    }

    template<typename T> 
    bool multiplyOver4thDimensionExcept(const hoNDArray<T>& x4D, const hoNDArray<T>& y5D, size_t n, hoNDArray<T>& r, bool copyY2R)
    {
        try
        {
            boost::shared_ptr< std::vector<size_t> > dimX = x4D.get_dimensions();
            boost::shared_ptr< std::vector<size_t> > dimY = y5D.get_dimensions();

            size_t NDim = dimY->size();

            GADGET_CHECK_RETURN_FALSE(dimX->size()>=4);
            GADGET_CHECK_RETURN_FALSE((*dimX)[0]==(*dimY)[0]);
            GADGET_CHECK_RETURN_FALSE((*dimX)[1]==(*dimY)[1]);
            GADGET_CHECK_RETURN_FALSE((*dimX)[2]==(*dimY)[2]);

            const T* pX = x4D.begin();
            const T* pY = y5D.begin();
            T* pR = r.begin();

            if ( (pR!=pY) && (!r.dimensions_equal(dimY.get())) )
            {
                r.create(dimY);
                pR = r.begin();
            }

            size_t RO = (*dimX)[0];
            size_t E1 = (*dimX)[1];
            size_t CHA = (*dimX)[2];

            int t, N3D = RO*E1*CHA;

            size_t N = (*dimY)[3];
            size_t num = x4D.get_number_of_elements()/(RO*E1*CHA);

            std::vector<size_t> dim3D(3);
            dim3D[0] = RO;
            dim3D[1] = E1;
            dim3D[2] = CHA;

            #pragma omp parallel for default(none) private(t) shared(N3D, N, dim3D, pX, pY, pR, num, n, copyY2R)
            for ( t=0; t<(int)num; t++ )
            {
                hoNDArray<T> tmpX, tmpY, tmpR;
                tmpX.create(&dim3D, const_cast<T*>(pX+t*N3D));

                for ( int z=0; z<N; z++ )
                {
                    if ( z != n )
                    {
                        tmpY.create(&dim3D, const_cast<T*>(pY+t*N3D*N+z*N3D));
                        tmpR.create(&dim3D, pR+t*N3D*N+z*N3D);
                        multiply(tmpX, tmpY, tmpR);
                    }
                    else
                    {
                        if ( pR != pY )
                        {
                            if ( copyY2R )
                            {
                                memcpy(pR+t*N3D*N+z*N3D, const_cast<T*>(pY+t*N3D*N+z*N3D), sizeof(T)*N3D);
                            }
                        }
                    }
                }
            }
        }
        catch (...)
        {
            GADGET_ERROR_MSG("Errors in multiplyOver4thDimensionExcept(const hoNDArray<float>& x4D, const hoNDArray<float>& y5D, size_t n, hoNDArray<float>& r, bool copyY2R) ... ");
            return false;
        }
        return true;
    }

    template<typename T> 
    bool multiplyOver5thDimension(const hoNDArray<T>& x, const hoNDArray<T>& y, hoNDArray<T>& r)
    {
        try
        {
            boost::shared_ptr< std::vector<size_t> > dimX = x.get_dimensions();
            boost::shared_ptr< std::vector<size_t> > dimY = y.get_dimensions();

            size_t NDim = dimY->size();

            GADGET_CHECK_RETURN_FALSE(dimX->size()>=5);
            GADGET_CHECK_RETURN_FALSE((*dimX)[0]==(*dimY)[0]);
            GADGET_CHECK_RETURN_FALSE((*dimX)[1]==(*dimY)[1]);
            GADGET_CHECK_RETURN_FALSE((*dimX)[2]==(*dimY)[2]);
            GADGET_CHECK_RETURN_FALSE((*dimX)[3]==(*dimY)[3]);

            if ( !r.dimensions_equal(dimY.get()) )
            {
                r.create(dimY);
            }

            size_t RO = (*dimX)[0];
            size_t E1 = (*dimX)[1];
            size_t E2 = (*dimX)[2];
            size_t CHA = (*dimX)[3];

            int t, N4D = RO*E1*E2*CHA;

            size_t N = (*dimY)[4];
            size_t num = x.get_number_of_elements()/N4D;

            const T* pX = x.begin();
            const T* pY = y.begin();
            T* pR = r.begin();

            std::vector<size_t> dim4D(4);
            dim4D[0] = RO;
            dim4D[1] = E1;
            dim4D[2] = E2;
            dim4D[3] = CHA;

            #pragma omp parallel for default(none) private(t) shared(N4D, N, dim4D, pX, pY, pR, num)
            for ( t=0; t<(int)num; t++ )
            {
                hoNDArray<T> tmpX, tmpY, tmpR;
                tmpX.create(&dim4D, const_cast<T*>(pX+t*N4D));

                for ( int n=0; n<N; n++ )
                {
                    tmpY.create(&dim4D, const_cast<T*>(pY+t*N4D*N+n*N4D));
                    tmpR.create(&dim4D, pR+t*N4D*N+n*N4D);
                    multiply(tmpX, tmpY, tmpR);
                }
            }
        }
        catch (...)
        {
            GADGET_ERROR_MSG("Errors in multiplyOver5thDimension(const hoNDArray<float>& x, const hoNDArray<float>& y, hoNDArray<float>& r) ... ");
            return false;
        }
        return true;
    }

    template<typename T> 
    bool multiplyOver5thDimensionExcept(const hoNDArray<T>& x, const hoNDArray<T>& y, size_t n, hoNDArray<T>& r, bool copyY2R)
    {
        try
        {
            boost::shared_ptr< std::vector<size_t> > dimX = x.get_dimensions();
            boost::shared_ptr< std::vector<size_t> > dimY = y.get_dimensions();

            size_t NDim = dimY->size();

            GADGET_CHECK_RETURN_FALSE(dimX->size()>=5);
            GADGET_CHECK_RETURN_FALSE((*dimX)[0]==(*dimY)[0]);
            GADGET_CHECK_RETURN_FALSE((*dimX)[1]==(*dimY)[1]);
            GADGET_CHECK_RETURN_FALSE((*dimX)[2]==(*dimY)[2]);
            GADGET_CHECK_RETURN_FALSE((*dimX)[3]==(*dimY)[3]);

            const T* pX = x.begin();
            const T* pY = y.begin();
            T* pR = r.begin();

            if ( (pR!=pY) && (!r.dimensions_equal(dimY.get())) )
            {
                r.create(dimY);
                pR = r.begin();
            }

            size_t RO = (*dimX)[0];
            size_t E1 = (*dimX)[1];
            size_t E2 = (*dimX)[2];
            size_t CHA = (*dimX)[3];

            int t, N4D = RO*E1*E2*CHA;

            size_t N = (*dimY)[4];
            size_t num = x.get_number_of_elements()/N4D;

            std::vector<size_t> dim4D(4);
            dim4D[0] = RO;
            dim4D[1] = E1;
            dim4D[2] = E2;
            dim4D[3] = CHA;

            #pragma omp parallel for default(none) private(t) shared(N4D, dim4D, pX, pY, pR, num, n, N, copyY2R)
            for ( t=0; t<(int)num; t++ )
            {
                hoNDArray<T> tmpX, tmpY, tmpR;
                tmpX.create(&dim4D, const_cast<T*>(pX+t*N4D));

                for ( int z=0; z<N; z++ )
                {
                    if ( z != n )
                    {
                        tmpY.create(&dim4D, const_cast<T*>(pY+t*N4D*N+z*N4D));
                        tmpR.create(&dim4D, pR+t*N4D*N+z*N4D);
                        multiply(tmpX, tmpY, tmpR);
                    }
                    else
                    {
                        if ( pR != pY )
                        {
                            if ( copyY2R )
                            {
                                memcpy(pR+t*N4D*N+z*N4D, const_cast<T*>(pY+t*N4D*N+z*N4D), sizeof(T)*N4D);
                            }
                        }
                    }
                }
            }
        }
        catch (...)
        {
            GADGET_ERROR_MSG("Errors in multiplyOver5thDimensionExcept(const hoNDArray<T>& x, const hoNDArray<T>& y, size_t n, hoNDArray<T>& r, bool copyY2R) ... ");
            return false;
        }
        return true;
    }

    template <typename T> 
    bool multipleAdd(const hoNDArray<T>& x, const hoNDArray<T>& y, hoNDArray<T>& r)
    {
        GADGET_DEBUG_CHECK_RETURN_FALSE(x.get_number_of_elements()<=y.get_number_of_elements());
        if ( r.get_number_of_elements()!=y.get_number_of_elements())
        {
            r = y;
        }

        int Nx = x.get_number_of_elements();
        int N = y.get_number_of_elements() / Nx;

        int n;

        if ( typeid(T)==typeid(float) )
        {
            #ifdef GCC_OLD_FLAG
                #pragma omp parallel for default(none) private(n) shared(Nx, N)
            #else
                #pragma omp parallel for default(none) private(n) shared(x, y, r, Nx, N)
            #endif
            for ( n=0; n<N; n++ )
            {
                vsAdd(x.get_number_of_elements(), reinterpret_cast<const float*>(x.begin()), reinterpret_cast<const float*>(y.begin()+n*Nx), reinterpret_cast<float*>(r.begin()+n*Nx));
            }
        }
        else if ( typeid(T)==typeid(double) )
        {
            #ifdef GCC_OLD_FLAG
                #pragma omp parallel for default(none) private(n) shared(Nx, N)
            #else
                #pragma omp parallel for default(none) private(n) shared(x, y, r, Nx, N)
            #endif
            for ( n=0; n<N; n++ )
            {
                vdAdd(x.get_number_of_elements(), reinterpret_cast<const double*>(x.begin()), reinterpret_cast<const double*>(y.begin()+n*Nx), reinterpret_cast<double*>(r.begin()+n*Nx));
            }
        }
        else if ( typeid(T)==typeid(GT_Complex8) )
        {
            #ifdef GCC_OLD_FLAG
                #pragma omp parallel for default(none) private(n) shared(Nx, N)
            #else
                #pragma omp parallel for default(none) private(n) shared(x, y, r, Nx, N)
            #endif
            for ( n=0; n<N; n++ )
            {
                vcAdd(x.get_number_of_elements(), reinterpret_cast<const MKL_Complex8*>(x.begin()), reinterpret_cast<const MKL_Complex8*>(y.begin()+n*Nx), reinterpret_cast<MKL_Complex8*>(r.begin()+n*Nx));
            }
        }
        else if ( typeid(T)==typeid(GT_Complex16) )
        {
            #ifdef GCC_OLD_FLAG
                #pragma omp parallel for default(none) private(n) shared(Nx, N)
            #else
                #pragma omp parallel for default(none) private(n) shared(x, y, r, Nx, N)
            #endif
            for ( n=0; n<N; n++ )
            {
                vzAdd(x.get_number_of_elements(), reinterpret_cast<const MKL_Complex16*>(x.begin()), reinterpret_cast<const MKL_Complex16*>(y.begin()+n*Nx), reinterpret_cast<MKL_Complex16*>(r.begin()+n*Nx));
            }
        }
        else
        {
            GADGET_ERROR_MSG("multipleAdd : unsupported type " << typeid(T).name());
            return false;
        }

        GADGET_CHECK_RETURN_FALSE(vmlGetErrStatus()==0);

        return true;
    }

    template <typename T> 
    bool multipleMultiply(const hoNDArray<T>& x, const hoNDArray<T>& y, hoNDArray<T>& r)
    {
        GADGET_DEBUG_CHECK_RETURN_FALSE(x.get_number_of_elements()<=y.get_number_of_elements());
        if ( r.get_number_of_elements()!=y.get_number_of_elements())
        {
            r = y;
        }

        int Nx = x.get_number_of_elements();
        int N = y.get_number_of_elements() / Nx;

        int n;

        if ( typeid(T)==typeid(float) )
        {
            #ifdef GCC_OLD_FLAG
                #pragma omp parallel for default(none) private(n) shared(Nx, N)
            #else
                #pragma omp parallel for default(none) private(n) shared(x, y, r, Nx, N)
            #endif
            for ( n=0; n<N; n++ )
            {
                vsMul(x.get_number_of_elements(), reinterpret_cast<const float*>(x.begin()), reinterpret_cast<const float*>(y.begin()+n*Nx), reinterpret_cast<float*>(r.begin()+n*Nx));
            }
        }
        else if ( typeid(T)==typeid(double) )
        {
            #ifdef GCC_OLD_FLAG
                #pragma omp parallel for default(none) private(n) shared(Nx, N)
            #else
                #pragma omp parallel for default(none) private(n) shared(x, y, r, Nx, N)
            #endif
            for ( n=0; n<N; n++ )
            {
                vdMul(x.get_number_of_elements(), reinterpret_cast<const double*>(x.begin()), reinterpret_cast<const double*>(y.begin()+n*Nx), reinterpret_cast<double*>(r.begin()+n*Nx));
            }
        }
        else if ( typeid(T)==typeid(GT_Complex8) )
        {
            #ifdef GCC_OLD_FLAG
                #pragma omp parallel for default(none) private(n) shared(Nx, N)
            #else
                #pragma omp parallel for default(none) private(n) shared(x, y, r, Nx, N)
            #endif
            for ( n=0; n<N; n++ )
            {
                vcMul(x.get_number_of_elements(), reinterpret_cast<const MKL_Complex8*>(x.begin()), reinterpret_cast<const MKL_Complex8*>(y.begin()+n*Nx), reinterpret_cast<MKL_Complex8*>(r.begin()+n*Nx));
            }
        }
        else if ( typeid(T)==typeid(GT_Complex16) )
        {
            #ifdef GCC_OLD_FLAG
                #pragma omp parallel for default(none) private(n) shared(Nx, N)
            #else
                #pragma omp parallel for default(none) private(n) shared(x, y, r, Nx, N)
            #endif
            for ( n=0; n<N; n++ )
            {
                vzMul(x.get_number_of_elements(), reinterpret_cast<const MKL_Complex16*>(x.begin()), reinterpret_cast<const MKL_Complex16*>(y.begin()+n*Nx), reinterpret_cast<MKL_Complex16*>(r.begin()+n*Nx));
            }
        }
        else
        {
            GADGET_ERROR_MSG("multipleMultiply : unsupported type " << typeid(T).name());
            return false;
        }

        GADGET_CHECK_RETURN_FALSE(vmlGetErrStatus()==0);

        return true;
    }

    template <typename T> 
    bool cropUpTo10DArray(const hoNDArray<T>& x, hoNDArray<T>& r, const std::vector<size_t>& startND, std::vector<size_t>& size)
    {
        GADGET_CHECK_RETURN_FALSE( startND.size() == size.size() );
        GADGET_CHECK_RETURN_FALSE( startND.size() <= 10 );

        r.create(&size);
        if ( r.get_number_of_elements() == x.get_number_of_elements() )
        {
            r = x;
            return true;
        }

        std::vector<size_t> start(10, 0);
        std::vector<size_t> end(10, 0);

        size_t ii;
        for ( ii=0; ii<startND.size(); ii++ )
        {
            start[ii] = startND[ii];
            end[ii] = start[ii] + size[ii] - 1;
            GADGET_CHECK_RETURN_FALSE(end[ii] < x.get_size(ii));
        }

        // [Ro E1 Cha Slice E2 Con Phase Rep Set Seg]
        size_t ro, e1, cha, n, s, con, phs, rep, set, seg;

        std::vector<size_t> srcInd(10), dstInd(10);

        for ( seg=start[9]; seg<=end[9]; seg++ )
        {
            srcInd[9] = seg; dstInd[9] = seg-start[9];

            for ( set=start[8]; set<=end[8]; set++ )
            {
                srcInd[8] = set; dstInd[8] = set-start[8];

                for ( rep=start[7]; rep<=end[7]; rep++ )
                {
                    srcInd[7] = rep; dstInd[7] = rep-start[7];

                    for ( phs=start[6]; phs<=end[6]; phs++ )
                    {
                        srcInd[6] = phs; dstInd[6] = phs-start[6];

                        for ( con=start[5]; con<=end[5]; con++ )
                        {
                            srcInd[5] = con; dstInd[5] = con-start[5];

                            for ( s=start[4]; s<=end[4]; s++ )
                            {
                                srcInd[4] = s; dstInd[4] = s-start[4];

                                for ( n=start[3]; n<=end[3]; n++ )
                                {
                                    srcInd[3] = n; dstInd[3] = n-start[3];

                                    for ( cha=start[2]; cha<=end[2]; cha++ )
                                    {
                                        srcInd[2] = cha; dstInd[2] = cha-start[2];

                                        for ( e1=start[1]; e1<=end[1]; e1++ )
                                        {
                                            srcInd[1] = e1; dstInd[1] = e1-start[1];

                                            srcInd[0] = start[0];
                                            dstInd[0] = 0;

                                            int offsetSrc = x.calculate_offset(srcInd);
                                            int offsetDst = r.calculate_offset(dstInd);

                                            memcpy(r.begin()+offsetDst, x.begin()+offsetSrc, sizeof(T)*(end[0]-start[0]+1));

                                            /*for ( ro=start[0]; ro<=end[0]; ro++ )
                                            {
                                                srcInd[0] = ro;
                                                dstInd[0] = ro-start[0];

                                                int offsetSrc = x.calculate_offset(srcInd);
                                                int offsetDst = r.calculate_offset(dstInd);

                                                r(offsetDst) = x(offsetSrc);
                                            }*/
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        return true;
    }

    template <typename T> 
    bool setSubArrayUpTo10DArray(const hoNDArray<T>& x, hoNDArray<T>& r, const std::vector<size_t>& startND, std::vector<size_t>& size)
    {
        GADGET_CHECK_RETURN_FALSE( startND.size() == size.size() );
        GADGET_CHECK_RETURN_FALSE( startND.size() <= 10 );

        if ( r.get_number_of_elements() == x.get_number_of_elements() )
        {
            r = x;
            return true;
        }

        std::vector<size_t> start(10, 0);
        std::vector<size_t> end(10, 0);

        size_t ii;
        for ( ii=0; ii<startND.size(); ii++ )
        {
            start[ii] = startND[ii];
            end[ii] = start[ii] + size[ii] - 1;
            GADGET_CHECK_RETURN_FALSE(end[ii] < r.get_size(ii));
        }

        // [Ro E1 Cha Slice E2 Con Phase Rep Set Seg]
        size_t ro, e1, cha, n, s, con, phs, rep, set, seg;

        std::vector<size_t> srcInd(10), dstInd(10);

        for ( seg=start[9]; seg<=end[9]; seg++ )
        {
            dstInd[9] = seg; srcInd[9] = seg-start[9];

            for ( set=start[8]; set<=end[8]; set++ )
            {
                dstInd[8] = set; srcInd[8] = set-start[8];

                for ( rep=start[7]; rep<=end[7]; rep++ )
                {
                    dstInd[7] = rep; srcInd[7] = rep-start[7];

                    for ( phs=start[6]; phs<=end[6]; phs++ )
                    {
                        dstInd[6] = phs; srcInd[6] = phs-start[6];

                        for ( con=start[5]; con<=end[5]; con++ )
                        {
                            dstInd[5] = con; srcInd[5] = con-start[5];

                            for ( s=start[4]; s<=end[4]; s++ )
                            {
                                dstInd[4] = s; srcInd[4] = s-start[4];

                                for ( n=start[3]; n<=end[3]; n++ )
                                {
                                    dstInd[3] = n; srcInd[3] = n-start[3];

                                    for ( cha=start[2]; cha<=end[2]; cha++ )
                                    {
                                        dstInd[2] = cha; srcInd[2] = cha-start[2];

                                        for ( e1=start[1]; e1<=end[1]; e1++ )
                                        {
                                            dstInd[1] = e1; srcInd[1] = e1-start[1];

                                            dstInd[0] = start[0];
                                            srcInd[0] = 0;

                                            int offsetSrc = x.calculate_offset(srcInd);
                                            int offsetDst = r.calculate_offset(dstInd);

                                            memcpy(r.begin()+offsetDst, x.begin()+offsetSrc, sizeof(T)*(end[0]-start[0]+1));
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        return true;
    }

    template<typename T> 
    bool stdOver3rdDimension(const hoNDArray<T>& x, hoNDArray<T>& std, bool NMinusOne)
    {
        try
        {
            GADGET_DEBUG_CHECK_RETURN_FALSE(x.get_number_of_dimensions() >= 3);

            size_t RO = x.get_size(0);
            size_t E1 = x.get_size(1);
            size_t CHA = x.get_size(2);

            int num = (int)x.get_number_of_elements() / (RO*E1*CHA);

            boost::shared_ptr< std::vector<size_t> > dim = x.get_dimensions();

            std::vector<size_t> dimStd(*dim);
            dimStd.erase(dimStd.begin()+2);
            std.create(&dimStd);

            std::vector<size_t> dim3D(3);
            dim3D[0] = RO;
            dim3D[1] = E1;
            dim3D[2] = CHA;

            T S(CHA);
            if ( NMinusOne )
            {
                S = T(CHA-1);
            }

            T v(0), v1(0);
            T S2 = T(1.0)/S;
            T S3 = T(1.0)/T(CHA);

            int n;

            #ifdef GCC_OLD_FLAG
                #pragma omp parallel for default(none) private(n) shared(num, RO, E1, CHA, S, S2, S3, v, v1)
            #else
                #pragma omp parallel for default(none) private(n) shared(num, RO, E1, CHA, x, std, S, S2, S3, v, v1)
            #endif
            for ( n=0; n<num; n++ )
            {
                hoNDArray<T> xTmp(RO, E1, CHA, const_cast<T*>(x.begin()+n*RO*E1*CHA));
                hoNDArray<T> mean(RO, E1);

                size_t ro, e1, cha;
                for ( cha=0; cha<CHA; cha++ )
                {
                    for ( e1=0; e1<E1; e1++ )
                    {
                        for ( ro=0; ro<RO; ro++ )
                        {
                            mean(ro+e1*RO) += xTmp(cha*RO*E1+e1*RO+ro)*S3;
                        }
                    }
                }

                for ( e1=0; e1<E1; e1++ )
                {
                    for ( ro=0; ro<RO; ro++ )
                    {
                        int ind = e1*RO+ro;

                        v = 0; v1 = 0;
                        for ( cha=0; cha<CHA; cha++ )
                        {
                            v1 = std::abs(xTmp(cha*RO*E1+ind)-mean(ind));
                            v += v1*v1;
                        }

                        v /= S;
                        std(ind+n*RO*E1) = std::sqrt(v);
                    }
                }
            }
        }
        catch(...)
        {
            GADGET_ERROR_MSG("Errors in stdOver3rdDimension(const hoNDArray<T>& x, hoNDArray<T>& std, bool NMinusOne) ... ");
            return false;
        }

        return true;
    }

    /*template<typename T> 
    bool permuteLastTwoDimensions(const hoNDArray<T>& x, hoNDArray<T>& r)
    {
        try
        {
            boost::shared_ptr< std::vector<size_t> > dimX = x.get_dimensions();

            size_t NDim = dimX->size();

            if ( NDim <= 2 )
            {
                r = x;
                return true;
            }

            size_t E1 = x.get_size(NDim-2);
            size_t E2 = x.get_size(NDim-1);

            std::vector<size_t> dimR(*dimX);
            dimR[NDim-2] = E2;
            dimR[NDim-1] = E1;

            r.create(&dimR);

            size_t N = x.get_number_of_elements()/E1/E2;

            const T* pX = x.begin();
            T* pR = r.begin();

            int e2;

            #pragma omp parallel for default(none) private(e2) shared(E2, E1, pR, pX, N)
            for ( e2=0; e2<(int)E2; e2++ )
            {
                for ( size_t e1=0; e1<E1; e1++ )
                {
                    memcpy(pR+e1*N*E2+e2*N, pX+e2*N*E1+e1*N, sizeof(T)*N);
                }
            }
        }
        catch (...)
        {
            GADGET_ERROR_MSG("Errors in permuteLastTwoDimensions(const hoNDArray<T>& x, hoNDArray<T>& r) ... ");
            return false;
        }
        return true;
    }*/

    template<typename T> 
    bool cropOver3rdDimension(const hoNDArray<T>& x, hoNDArray<T>& r, size_t start, size_t end)
    {
        try
        {
            boost::shared_ptr< std::vector<size_t> > dimX = x.get_dimensions();

            size_t NDim = dimX->size();

            if ( NDim <= 2 )
            {
                r = x;
                return true;
            }

            size_t RO = x.get_size(0);
            size_t E1 = x.get_size(1);
            size_t E2 = x.get_size(2);

            size_t E2_R = end-start+1;

            if ( E2 <= E2_R )
            {
                r = x;
                return true;
            }

            std::vector<size_t> dimR(*dimX);
            dimR[2] = E2_R;

            r.create(&dimR);

            size_t N2D = RO*E1;
            size_t N3D = RO*E1*E2;
            size_t N3D_R = RO*E1*E2_R;

            size_t N = x.get_number_of_elements()/N3D;

            const T* pX = x.begin();
            T* pR = r.begin();

            size_t n;
            for ( n=0; n<N; n++ )
            {
                int e2;
                #pragma omp parallel for default(none) private(e2) shared(N2D, N3D, N3D_R, pX, pR, RO, E1, E2, n, start, end)
                for ( e2=start; e2<=end; e2++ )
                {
                    memcpy(pR+n*N3D_R+(e2-start)*N2D, pX+n*N3D+e2*N2D, sizeof(T)*N2D);
                }
            }
        }
        catch (...)
        {
            GADGET_ERROR_MSG("Errors in cropOver3rdDimension(const hoNDArray<T>& x, hoNDArray<T>& r, size_t start, size_t end) ... ");
            return false;
        }
        return true;
    }

    template<typename T> bool setSubArrayOver3rdDimension(const hoNDArray<T>& x, hoNDArray<T>& r, size_t start, size_t end)
    {
        try
        {
            boost::shared_ptr< std::vector<size_t> > dimR = r.get_dimensions();

            size_t NDim = dimR->size();

            if ( NDim <= 2 )
            {
                r = x;
                return true;
            }

            size_t RO = r.get_size(0);
            size_t E1 = r.get_size(1);
            size_t E2 = r.get_size(2);

            size_t E2_X = end-start+1;
            GADGET_CHECK_RETURN_FALSE( E2_X == x.get_size(2) );

            if ( E2_X >= E2 )
            {
                r = x;
                return true;
            }

            size_t N2D = RO*E1;
            size_t N3D = RO*E1*E2;
            size_t N3D_X = RO*E1*E2_X;

            size_t N = r.get_number_of_elements()/N3D;

            const T* pX = x.begin();
            T* pR = r.begin();

            size_t n;
            for ( n=0; n<N; n++ )
            {
                int e2;
                #pragma omp parallel for default(none) private(e2) shared(N2D, N3D, N3D_X, pX, pR, RO, E1, E2, n, start, end)
                for ( e2=start; e2<=end; e2++ )
                {
                    memcpy(pR+n*N3D+e2*N2D, pX+n*N3D_X+(e2-start)*N2D, sizeof(T)*N2D);
                }
            }
        }
        catch (...)
        {
            GADGET_ERROR_MSG("Errors in setSubArrayOver3rdDimension(const hoNDArray<T>& x, hoNDArray<T>& r, size_t start, size_t end) ... ");
            return false;
        }
        return true;
    }

    template<typename T> 
    bool permuteE2To3rdDimension(const hoNDArray<T>& x, hoNDArray<T>& r)
    {
        try
        {
            boost::shared_ptr< std::vector<size_t> > dimX = x.get_dimensions();

            size_t NDim = dimX->size();

            if ( NDim <= 5 )
            {
                r = x;
                return true;
            }

            size_t RO = x.get_size(0);
            size_t E1 = x.get_size(1);
            size_t CHA = x.get_size(2);
            size_t SLC = x.get_size(3);
            size_t E2 = x.get_size(4);

            std::vector<size_t> dimR(*dimX);
            dimR[2] = E2;
            dimR[3] = CHA;
            dimR[4] = SLC;

            r.create(&dimR);

            size_t N2D = RO*E1;
            size_t N5D = RO*E1*CHA*E2*SLC;

            size_t N = x.get_number_of_elements()/N5D;

            const T* pX = x.begin();
            T* pR = r.begin();

            size_t n;
            for ( n=0; n<N; n++ )
            {
                int e2;
                #pragma omp parallel for default(none) private(e2) shared(N5D, N2D, pX, pR, CHA, SLC, E2, n)
                for ( e2=0; e2<E2; e2++ )
                {
                    for ( size_t slc=0; slc<SLC; slc++ )
                    {
                        for ( size_t cha=0; cha<CHA; cha++ )
                        {
                            memcpy(pR+n*N5D+slc*CHA*E2*N2D+cha*E2*N2D+e2*N2D, pX+n*N5D+e2*SLC*CHA*N2D+slc*CHA*N2D+cha*N2D, sizeof(T)*N2D);
                        }
                    }
                }
            }
        }
        catch (...)
        {
            GADGET_ERROR_MSG("Errors in permuteE2To3rdDimension(const hoNDArray<T>& x, hoNDArray<T>& r) ... ");
            return false;
        }
        return true;
    }

    template<typename T> 
    bool permuteE2To5thDimension(const hoNDArray<T>& x, hoNDArray<T>& r)
    {
        try
        {
            boost::shared_ptr< std::vector<size_t> > dimX = x.get_dimensions();

            size_t NDim = dimX->size();

            if ( NDim < 5 )
            {
                r = x;
                return true;
            }

            size_t RO = x.get_size(0);
            size_t E1 = x.get_size(1);
            size_t E2 = x.get_size(2);
            size_t CHA = x.get_size(3);
            size_t SLC = x.get_size(4);

            std::vector<size_t> dimR(*dimX);
            dimR[2] = CHA;
            dimR[3] = SLC;
            dimR[4] = E2;

            r.create(&dimR);

            size_t N2D = RO*E1;
            size_t N5D = RO*E1*CHA*E2*SLC;

            size_t N = x.get_number_of_elements()/N5D;

            const T* pX = x.begin();
            T* pR = r.begin();

            size_t n;
            for ( n=0; n<N; n++ )
            {
                int e2;
                #pragma omp parallel for default(none) private(e2) shared(N5D, N2D, pX, pR, CHA, SLC, E2, n)
                for ( e2=0; e2<E2; e2++ )
                {
                    for ( size_t slc=0; slc<SLC; slc++ )
                    {
                        for ( size_t cha=0; cha<CHA; cha++ )
                        {
                            memcpy(pR+n*N5D+e2*SLC*CHA*N2D+slc*CHA*N2D+cha*N2D, pX+n*N5D+slc*CHA*E2*N2D+cha*E2*N2D+e2*N2D, sizeof(T)*N2D);
                        }
                    }
                }
            }
        }
        catch (...)
        {
            GADGET_ERROR_MSG("Errors in permuteE2To5thDimension(const hoNDArray<T>& x, hoNDArray<T>& r) ... ");
            return false;
        }
        return true;
    }

    template<typename T> 
    bool permuteROTo3rdDimensionFor3DRecon(const hoNDArray<T>& x, hoNDArray<T>& r)
    {
        try
        {
            boost::shared_ptr< std::vector<size_t> > dimX = x.get_dimensions();

            size_t NDim = dimX->size();

            if ( NDim < 3 )
            {
                r = x;
                return true;
            }

            size_t RO = x.get_size(0);
            size_t E1 = x.get_size(1);
            size_t E2 = x.get_size(2);

            std::vector<size_t> dimR(*dimX);
            dimR[0] = E1;
            dimR[1] = E2;
            dimR[2] = RO;

            r.create(&dimR);

            size_t N3D = RO*E1*E2;

            size_t N = x.get_number_of_elements()/N3D;

            const T* pX = x.begin();
            T* pR = r.begin();

            long long n;

            #pragma omp parallel for default(none) private(n) shared(RO, E1, E2, N, pR, N3D, pX)
            for ( n=0; n<(long long)N; n++ )
            {
                T* pRn = pR + n*N3D;
                T* pXn = const_cast<T*>(pX) + n*N3D;

                for ( size_t e2=0; e2<E2; e2++ )
                {
                    for ( size_t e1=0; e1<E1; e1++ )
                    {
                        for ( size_t ro=0; ro<RO; ro++ )
                        {
                            pRn[e1+e2*E1+ro*E1*E2] = pXn[ro+e1*RO+e2*RO*E1];
                        }
                    }
                }
            }
        }
        catch (...)
        {
            GADGET_ERROR_MSG("Errors in permuteROTo3rdDimensionFor3DRecon(const hoNDArray<T>& x, hoNDArray<T>& r) ... ");
            return false;
        }
        return true;
    }

    template<typename T> 
    bool permuteROTo4thDimensionFor3DRecon(const hoNDArray<T>& x, hoNDArray<T>& r)
    {
        try
        {
            boost::shared_ptr< std::vector<size_t> > dimX = x.get_dimensions();

            size_t NDim = dimX->size();

            if ( NDim < 4 )
            {
                r = x;
                return true;
            }

            size_t RO = x.get_size(0);
            size_t E1 = x.get_size(1);
            size_t E2 = x.get_size(2);
            size_t CHA = x.get_size(3);

            std::vector<size_t> dimR(*dimX);
            dimR[0] = E1;
            dimR[1] = E2;
            dimR[2] = CHA;
            dimR[3] = RO;

            r.create(&dimR);

            size_t N4D = RO*E1*E2*CHA;

            size_t N = x.get_number_of_elements()/N4D;

            const T* pX = x.begin();
            T* pR = r.begin();

            long long n;
            for ( n=0; n<(long long)N; n++ )
            {
                T* pRn = pR + n*N4D;
                T* pXn = const_cast<T*>(pX) + n*N4D;

                long long cha;

                #pragma omp parallel for default(none) private(cha) shared(RO, E1, E2, CHA, pXn, pRn)
                for ( cha=0; cha<(long long)CHA; cha++ )
                {
                    for ( size_t e2=0; e2<E2; e2++ )
                    {
                        for ( size_t e1=0; e1<E1; e1++ )
                        {
                            for ( size_t ro=0; ro<RO; ro++ )
                            {
                                pRn[e1+e2*E1+cha*E1*E2+ro*E1*E2*CHA] = pXn[ro+e1*RO+e2*RO*E1+cha*RO*E1*E2];
                            }
                        }
                    }
                }
            }
        }
        catch (...)
        {
            GADGET_ERROR_MSG("Errors in permuteROTo4thDimensionFor3DRecon(const hoNDArray<T>& x, hoNDArray<T>& r) ... ");
            return false;
        }
        return true;
    }

    template<typename T> 
    bool permuteROTo1stDimensionFor3DRecon(const hoNDArray<T>& x, hoNDArray<T>& r)
    {
        try
        {
            boost::shared_ptr< std::vector<size_t> > dimX = x.get_dimensions();

            size_t NDim = dimX->size();

            if ( NDim < 4 )
            {
                r = x;
                return true;
            }

            size_t E1 = x.get_size(0);
            size_t E2 = x.get_size(1);
            size_t CHA = x.get_size(2);
            size_t RO = x.get_size(3);

            std::vector<size_t> dimR(*dimX);
            dimR[0] = RO;
            dimR[1] = E1;
            dimR[2] = E2;
            dimR[3] = CHA;

            r.create(&dimR);

            size_t N4D = RO*E1*E2*CHA;

            size_t N = x.get_number_of_elements()/N4D;

            const T* pX = x.begin();
            T* pR = r.begin();

            long long n;
            for ( n=0; n<(long long)N; n++ )
            {
                T* pRn = pR + n*N4D;
                T* pXn = const_cast<T*>(pX) + n*N4D;

                long long cha;

                #pragma omp parallel for default(none) private(cha) shared(RO, E1, E2, CHA, pXn, pRn)
                for ( cha=0; cha<(long long)CHA; cha++ )
                {
                    for ( size_t e2=0; e2<E2; e2++ )
                    {
                        for ( size_t e1=0; e1<E1; e1++ )
                        {
                            size_t indRn = e1*RO+e2*RO*E1+cha*RO*E1*E2;
                            size_t indXn = e1+e2*E1+cha*E1*E2;
                            for ( size_t ro=0; ro<RO; ro++ )
                            {
                                pRn[ro+indRn] = pXn[indXn+ro*E1*E2*CHA];
                            }
                        }
                    }
                }
            }
        }
        catch (...)
        {
            GADGET_ERROR_MSG("Errors in permuteROTo1stDimensionFor3DRecon(const hoNDArray<T>& x, hoNDArray<T>& r) ... ");
            return false;
        }
        return true;
    }

    template<typename T> 
    bool permute3rdDimensionTo1stDimension(const hoNDArray<T>& x, hoNDArray<T>& r)
    {
        try
        {
            boost::shared_ptr< std::vector<size_t> > dimX = x.get_dimensions();

            size_t NDim = dimX->size();

            if ( NDim < 3 )
            {
                r = x;
                return true;
            }

            size_t RO = x.get_size(0);
            size_t E1 = x.get_size(1);
            size_t E2 = x.get_size(2);

            std::vector<size_t> dimR(*dimX);
            dimR[0] = E2;
            dimR[1] = RO;
            dimR[2] = E1;

            r.create(&dimR);

            size_t N3D = RO*E1*E2;

            size_t N = x.get_number_of_elements()/N3D;

            const T* pX = x.begin();
            T* pR = r.begin();

            long long n, e2;
            for ( n=0; n<(long long)N; n++ )
            {
                T* pRn = pR + n*N3D;
                T* pXn = const_cast<T*>(pX) + n*N3D;

                #pragma omp parallel for default(none) private(e2) shared(RO, E1, E2, pXn, pRn)
                for ( e2=0; e2<(long long)E2; e2++ )
                {
                    for ( size_t e1=0; e1<E1; e1++ )
                    {
                        size_t indRn = e2+e1*E2*RO;
                        size_t indXn = e1*RO+e2*RO*E1;
                        for ( size_t ro=0; ro<RO; ro++ )
                        {
                            pRn[ro*E2+indRn] = pXn[ro+indXn];
                        }
                    }
                }
            }
        }
        catch (...)
        {
            GADGET_ERROR_MSG("Errors in permute3rdDimensionTo1stDimension(const hoNDArray<T>& x, hoNDArray<T>& r) ... ");
            return false;
        }
        return true;
    }

    template<typename T> 
    bool permuteROTo5thDimensionFor3DRecon(const hoNDArray<T>& x, hoNDArray<T>& r)
    {
        try
        {
            boost::shared_ptr< std::vector<size_t> > dimX = x.get_dimensions();

            size_t NDim = dimX->size();

            if ( NDim < 5 )
            {
                r = x;
                return true;
            }

            size_t RO = x.get_size(0);
            size_t E1 = x.get_size(1);
            size_t E2 = x.get_size(2);
            size_t srcCHA = x.get_size(3);
            size_t dstCHA = x.get_size(4);

            std::vector<size_t> dimR(*dimX);
            dimR[0] = E1;
            dimR[1] = E2;
            dimR[2] = srcCHA;
            dimR[3] = dstCHA;
            dimR[4] = RO;

            r.create(&dimR);

            size_t N5D = RO*E1*E2*srcCHA*dstCHA;

            size_t N = x.get_number_of_elements()/N5D;

            const T* pX = x.begin();
            T* pR = r.begin();

            long long n;
            for ( n=0; n<(long long)N; n++ )
            {
                T* pRn = pR + n*N5D;
                T* pXn = const_cast<T*>(pX) + n*N5D;

                long long dcha;

                #pragma omp parallel for default(none) private(dcha) shared(RO, E1, E2, srcCHA, dstCHA, pXn, pRn)
                for ( dcha=0; dcha<(long long)dstCHA; dcha++ )
                {
                    for ( size_t scha=0; scha<(int)srcCHA; scha++ )
                    {
                        for ( size_t e2=0; e2<E2; e2++ )
                        {
                            for ( size_t e1=0; e1<E1; e1++ )
                            {
                                size_t indRn = e1+e2*E1+scha*E1*E2+dcha*E1*E2*srcCHA;
                                size_t indXn = e1*RO+e2*RO*E1+scha*RO*E1*E2+dcha*RO*E1*E2*srcCHA;
                                for ( size_t ro=0; ro<RO; ro++ )
                                {
                                    pRn[indRn+ro*E1*E2*srcCHA*dstCHA] = pXn[ro+indXn];
                                }
                            }
                        }
                    }
                }
            }
        }
        catch (...)
        {
            GADGET_ERROR_MSG("Errors in permuteROTo5thDimensionFor3DRecon(const hoNDArray<T>& x, hoNDArray<T>& r) ... ");
            return false;
        }
        return true;
    }

    template<typename T> 
    bool imageDomainUnwrapping2D(const hoNDArray<T>& x, const hoNDArray<T>& kernel, hoNDArray<T>& buf, hoNDArray<T>& y)
    {
        try
        {
            T* pX = const_cast<T*>(x.begin());
            T* ker = const_cast<T*>(kernel.begin());
            T* pY = y.begin();

            size_t ro = x.get_size(0);
            size_t e1 = x.get_size(1);
            size_t srcCHA = x.get_size(2);
            size_t dstCHA = kernel.get_size(3);

            if ( buf.get_number_of_elements() < ro*e1*srcCHA )
            {
                buf.create(ro, e1, srcCHA);
            }
            T* pBuf = buf.begin();

            long long dCha;

            //#pragma omp parallel default(shared)
            {
                //#ifdef WIN32
                //    int tid = omp_get_thread_num();
                //    DWORD_PTR mask = (1 << tid);
                //    // GADGET_MSG("thread id : " << tid << " - mask : " << mask);
                //    SetThreadAffinityMask( GetCurrentThread(), mask );
                //#endif // WIN32

                //#pragma omp for

                if ( typeid(T)==typeid(GT_Complex8) )
                {
                    for ( dCha=0; dCha<dstCHA; dCha++ )
                    {
                        vcMul(ro*e1*srcCHA, reinterpret_cast<MKL_Complex8*>(pX), 
                            reinterpret_cast<MKL_Complex8*>(ker+dCha*ro*e1*srcCHA), 
                            reinterpret_cast<MKL_Complex8*>(pBuf));

                        memcpy(pY+dCha*ro*e1, pBuf, sizeof(T)*ro*e1);
                        for ( size_t sCha=1; sCha<srcCHA; sCha++ )
                        {
                            vcAdd(ro*e1, reinterpret_cast<MKL_Complex8*>(pY+dCha*ro*e1), 
                                reinterpret_cast<MKL_Complex8*>(pBuf+sCha*ro*e1), 
                                reinterpret_cast<MKL_Complex8*>(pY+dCha*ro*e1));
                        }
                    }
                }
                else if ( typeid(T)==typeid(GT_Complex16) )
                {
                    for ( dCha=0; dCha<dstCHA; dCha++ )
                    {
                        vzMul(ro*e1*srcCHA, reinterpret_cast<MKL_Complex16*>(pX), 
                            reinterpret_cast<MKL_Complex16*>(ker+dCha*ro*e1*srcCHA), 
                            reinterpret_cast<MKL_Complex16*>(pBuf));

                        memcpy(pY+dCha*ro*e1, pBuf, sizeof(T)*ro*e1);
                        for ( size_t sCha=1; sCha<srcCHA; sCha++ )
                        {
                            vzAdd(ro*e1, reinterpret_cast<MKL_Complex16*>(pY+dCha*ro*e1), 
                                reinterpret_cast<MKL_Complex16*>(pBuf+sCha*ro*e1), 
                                reinterpret_cast<MKL_Complex16*>(pY+dCha*ro*e1));
                        }
                    }
                }
            }
        }
        catch (...)
        {
            GADGET_ERROR_MSG("Errors in imageDomainUnwrapping2D(const hoNDArray<T>& x, const hoNDArray<T>& ker, hoNDArray<T>& buf, hoNDArray<T>& y) ... ");
            return false;
        }
        return true;
    }

    template<typename T> 
    bool imageDomainUnwrapping2DT(const hoNDArray<T>& x, const hoNDArray<T>& kernel, hoNDArray<T>& buf, hoNDArray<T>& y)
    {
        try
        {
            size_t ro = x.get_size(0);
            size_t e1 = x.get_size(1);
            size_t srcCHA = x.get_size(2);
            size_t N = x.get_size(3);

            size_t dstCHA = kernel.get_size(3);
            size_t kerN = kernel.get_size(4);

            if ( buf.get_number_of_elements() < ro*e1*srcCHA )
            {
                buf.create(ro, e1, srcCHA);
            }
            T* pBuf = buf.begin();

            long long n, dCha;

            //#pragma omp parallel default(shared)
            {
                //#ifdef WIN32
                //    int tid = omp_get_thread_num();
                //    DWORD_PTR mask = (1 << tid);
                //    // GADGET_MSG("thread id : " << tid << " - mask : " << mask);
                //    SetThreadAffinityMask( GetCurrentThread(), mask );
                //#endif // WIN32

                //#pragma omp for

                if ( typeid(T)==typeid(GT_Complex8) )
                {
                    const T* pXN = x.begin();
                    T* pYN = y.begin();
                    T* pBufN = buf.begin();
                    const T* pKerN = kernel.begin();

                    #ifdef USE_OMP
                        omp_set_nested(1);
                    #endif // USE_OMP

                    //#pragma omp parallel for default(none) private(n) shared(N, ro, e1, srcCHA, dstCHA, kerN, pXN, pYN, pBufN, pKerN)
                    //for ( n=0; n<N; n++ )
                    //{
                    //    const T* ker = pKerN + n*ro*e1*srcCHA*dstCHA;
                    //    if ( kerN <= n )
                    //    {
                    //        ker = pKerN + (kerN-1)*ro*e1*srcCHA*dstCHA;
                    //    }

                    //    const T* pX = pXN + n*ro*e1*srcCHA;
                    //    T* pY = pYN + n*ro*e1*dstCHA;
                    //    T* pBuf =pBufN + n*ro*e1*srcCHA;

                    //    for ( size_t dCha=0; dCha<dstCHA; dCha++ )
                    //    {
                    //        vcMul(ro*e1*srcCHA, reinterpret_cast<const MKL_Complex8*>(pX), 
                    //            reinterpret_cast<const MKL_Complex8*>(ker+dCha*ro*e1*srcCHA), 
                    //            reinterpret_cast<MKL_Complex8*>(pBuf));

                    //        memcpy(pY+dCha*ro*e1, pBuf, sizeof(T)*ro*e1);
                    //        for ( size_t sCha=1; sCha<srcCHA; sCha++ )
                    //        {
                    //            vcAdd(ro*e1, reinterpret_cast<MKL_Complex8*>(pY+dCha*ro*e1), 
                    //                reinterpret_cast<MKL_Complex8*>(pBuf+sCha*ro*e1), 
                    //                reinterpret_cast<MKL_Complex8*>(pY+dCha*ro*e1));
                    //        }
                    //    }
                    //}

                    // #pragma omp parallel for default(none) private(dCha, n) shared(N, ro, e1, srcCHA, dstCHA, kerN, pXN, pYN, pBufN, pKerN)
                    for ( dCha=0; dCha<(long long)dstCHA; dCha++ )
                    {
                        for ( n=0; n<N; n++  )
                        {
                            const T* ker = pKerN + n*ro*e1*srcCHA*dstCHA;
                            if ( kerN <= n )
                            {
                                ker = pKerN + (kerN-1)*ro*e1*srcCHA*dstCHA;
                            }

                            const T* pX = pXN + n*ro*e1*srcCHA;
                            T* pBuf =pBufN + n*ro*e1*srcCHA;

                            vcMul(ro*e1*srcCHA, reinterpret_cast<const MKL_Complex8*>(pX), 
                                reinterpret_cast<const MKL_Complex8*>(ker+dCha*ro*e1*srcCHA), 
                                reinterpret_cast<MKL_Complex8*>(pBuf));
                        //}

                        //for ( n=0; n<N; n++  )
                        //{
                            T* pY = pYN + n*ro*e1*dstCHA;
                            //T* pBuf =pBufN + n*ro*e1*srcCHA;

                            memcpy(pY+dCha*ro*e1, pBuf, sizeof(T)*ro*e1);
                            for ( size_t sCha=1; sCha<srcCHA; sCha++ )
                            {
                                vcAdd(ro*e1, reinterpret_cast<MKL_Complex8*>(pY+dCha*ro*e1), 
                                    reinterpret_cast<MKL_Complex8*>(pBuf+sCha*ro*e1), 
                                    reinterpret_cast<MKL_Complex8*>(pY+dCha*ro*e1));
                            }
                        }
                    }
                }
                else if ( typeid(T)==typeid(GT_Complex16) )
                {
                    for ( n=0; n<N; n++ )
                    {
                        const T* ker = kernel.begin() + n*ro*e1*srcCHA*dstCHA;
                        if ( kerN <= n )
                        {
                            ker = kernel.begin() + (kerN-1)*ro*e1*srcCHA*dstCHA;
                        }

                        const T* pX = x.begin() + n*ro*e1*srcCHA;
                        T* pY = y.begin() + n*ro*e1*dstCHA;

                        for ( size_t dCha=0; dCha<dstCHA; dCha++ )
                        {
                            vzMul(ro*e1*srcCHA, reinterpret_cast<const MKL_Complex16*>(pX), 
                                reinterpret_cast<const MKL_Complex16*>(ker+dCha*ro*e1*srcCHA), 
                                reinterpret_cast<MKL_Complex16*>(pBuf));

                            memcpy(pY+dCha*ro*e1, pBuf, sizeof(T)*ro*e1);
                            for ( size_t sCha=1; sCha<srcCHA; sCha++ )
                            {
                                vzAdd(ro*e1, reinterpret_cast<MKL_Complex16*>(pY+dCha*ro*e1), 
                                    reinterpret_cast<MKL_Complex16*>(pBuf+sCha*ro*e1), 
                                    reinterpret_cast<MKL_Complex16*>(pY+dCha*ro*e1));
                            }
                        }
                    }
                }
            }
        }
        catch (...)
        {
            GADGET_ERROR_MSG("Errors in imageDomainUnwrapping2DT(const hoNDArray<T>& x, const hoNDArray<T>& ker, hoNDArray<T>& buf, hoNDArray<T>& y) ... ");
            return false;
        }
        return true;
    }

    template EXPORTCPUCOREMATH bool sumOverLastDimension(const hoNDArray<float>& x, hoNDArray<float>& r);
    template EXPORTCPUCOREMATH bool sumOverLastDimension(const hoNDArray<double>& x, hoNDArray<double>& r);
    template EXPORTCPUCOREMATH bool sumOverLastDimension(const hoNDArray<GT_Complex8>& x, hoNDArray<GT_Complex8>& r);
    template EXPORTCPUCOREMATH bool sumOverLastDimension(const hoNDArray<GT_Complex16>& x, hoNDArray<GT_Complex16>& r);

    template EXPORTCPUCOREMATH bool sumOverSecondLastDimension(const hoNDArray<float>& x, hoNDArray<float>& r);
    template EXPORTCPUCOREMATH bool sumOverSecondLastDimension(const hoNDArray<double>& x, hoNDArray<double>& r);
    template EXPORTCPUCOREMATH bool sumOverSecondLastDimension(const hoNDArray<GT_Complex8>& x, hoNDArray<GT_Complex8>& r);
    template EXPORTCPUCOREMATH bool sumOverSecondLastDimension(const hoNDArray<GT_Complex16>& x, hoNDArray<GT_Complex16>& r);

    template EXPORTCPUCOREMATH bool multiplyOverLastDimension(const hoNDArray<float>& x, const hoNDArray<float>& y, hoNDArray<float>& r);
    template EXPORTCPUCOREMATH bool multiplyOverLastDimension(const hoNDArray<double>& x, const hoNDArray<double>& y, hoNDArray<double>& r);
    template EXPORTCPUCOREMATH bool multiplyOverLastDimension(const hoNDArray<GT_Complex8>& x, const hoNDArray<GT_Complex8>& y, hoNDArray<GT_Complex8>& r);
    template EXPORTCPUCOREMATH bool multiplyOverLastDimension(const hoNDArray<GT_Complex16>& x, const hoNDArray<GT_Complex16>& y, hoNDArray<GT_Complex16>& r);

    template EXPORTCPUCOREMATH bool divideOverLastDimension(const hoNDArray<float>& x, const hoNDArray<float>& y, hoNDArray<float>& r);
    template EXPORTCPUCOREMATH bool divideOverLastDimension(const hoNDArray<double>& x, const hoNDArray<double>& y, hoNDArray<double>& r);
    template EXPORTCPUCOREMATH bool divideOverLastDimension(const hoNDArray<GT_Complex8>& x, const hoNDArray<GT_Complex8>& y, hoNDArray<GT_Complex8>& r);
    template EXPORTCPUCOREMATH bool divideOverLastDimension(const hoNDArray<GT_Complex16>& x, const hoNDArray<GT_Complex16>& y, hoNDArray<GT_Complex16>& r);

    template EXPORTCPUCOREMATH bool sumOver1stDimension(const hoNDArray<float>& x, hoNDArray<float>& r);
    template EXPORTCPUCOREMATH bool sumOver1stDimension(const hoNDArray<double>& x, hoNDArray<double>& r);
    template EXPORTCPUCOREMATH bool sumOver1stDimension(const hoNDArray<GT_Complex8>& x, hoNDArray<GT_Complex8>& r);
    template EXPORTCPUCOREMATH bool sumOver1stDimension(const hoNDArray<GT_Complex16>& x, hoNDArray<GT_Complex16>& r);

    template EXPORTCPUCOREMATH bool sumOver2ndDimension(const hoNDArray<float>& x, hoNDArray<float>& r);
    template EXPORTCPUCOREMATH bool sumOver2ndDimension(const hoNDArray<double>& x, hoNDArray<double>& r);
    template EXPORTCPUCOREMATH bool sumOver2ndDimension(const hoNDArray<GT_Complex8>& x, hoNDArray<GT_Complex8>& r);
    template EXPORTCPUCOREMATH bool sumOver2ndDimension(const hoNDArray<GT_Complex16>& x, hoNDArray<GT_Complex16>& r);

    template EXPORTCPUCOREMATH bool sumOver3rdDimension(const hoNDArray<float>& x, hoNDArray<float>& r);
    template EXPORTCPUCOREMATH bool sumOver3rdDimension(const hoNDArray<double>& x, hoNDArray<double>& r);
    template EXPORTCPUCOREMATH bool sumOver3rdDimension(const hoNDArray<GT_Complex8>& x, hoNDArray<GT_Complex8>& r);
    template EXPORTCPUCOREMATH bool sumOver3rdDimension(const hoNDArray<GT_Complex16>& x, hoNDArray<GT_Complex16>& r);

    template EXPORTCPUCOREMATH bool sumOver4thDimension(const hoNDArray<float>& x, hoNDArray<float>& r);
    template EXPORTCPUCOREMATH bool sumOver4thDimension(const hoNDArray<double>& x, hoNDArray<double>& r);
    template EXPORTCPUCOREMATH bool sumOver4thDimension(const hoNDArray<GT_Complex8>& x, hoNDArray<GT_Complex8>& r);
    template EXPORTCPUCOREMATH bool sumOver4thDimension(const hoNDArray<GT_Complex16>& x, hoNDArray<GT_Complex16>& r);

    template EXPORTCPUCOREMATH bool sumOver5thDimension(const hoNDArray<float>& x, hoNDArray<float>& r);
    template EXPORTCPUCOREMATH bool sumOver5thDimension(const hoNDArray<double>& x, hoNDArray<double>& r);
    template EXPORTCPUCOREMATH bool sumOver5thDimension(const hoNDArray<GT_Complex8>& x, hoNDArray<GT_Complex8>& r);
    template EXPORTCPUCOREMATH bool sumOver5thDimension(const hoNDArray<GT_Complex16>& x, hoNDArray<GT_Complex16>& r);

    template EXPORTCPUCOREMATH bool multiplyOver3rdDimension(const hoNDArray<float>& x3D, const hoNDArray<float>& y4D, hoNDArray<float>& r);
    template EXPORTCPUCOREMATH bool multiplyOver3rdDimension(const hoNDArray<double>& x3D, const hoNDArray<double>& y4D, hoNDArray<double>& r);
    template EXPORTCPUCOREMATH bool multiplyOver3rdDimension(const hoNDArray<GT_Complex8>& x3D, const hoNDArray<GT_Complex8>& y4D, hoNDArray<GT_Complex8>& r);
    template EXPORTCPUCOREMATH bool multiplyOver3rdDimension(const hoNDArray<GT_Complex16>& x3D, const hoNDArray<GT_Complex16>& y4D, hoNDArray<GT_Complex16>& r);

    template EXPORTCPUCOREMATH bool multiplyOver4thDimension(const hoNDArray<float>& x4D, const hoNDArray<float>& y5D, hoNDArray<float>& r);
    template EXPORTCPUCOREMATH bool multiplyOver4thDimension(const hoNDArray<double>& x4D, const hoNDArray<double>& y5D, hoNDArray<double>& r);
    template EXPORTCPUCOREMATH bool multiplyOver4thDimension(const hoNDArray<GT_Complex8>& x4D, const hoNDArray<GT_Complex8>& y5D, hoNDArray<GT_Complex8>& r);
    template EXPORTCPUCOREMATH bool multiplyOver4thDimension(const hoNDArray<GT_Complex16>& x4D, const hoNDArray<GT_Complex16>& y5D, hoNDArray<GT_Complex16>& r);

    template EXPORTCPUCOREMATH bool multiplyOver4thDimensionExcept(const hoNDArray<float>& x4D, const hoNDArray<float>& y5D, size_t n, hoNDArray<float>& r, bool copyY2R);
    template EXPORTCPUCOREMATH bool multiplyOver4thDimensionExcept(const hoNDArray<double>& x4D, const hoNDArray<double>& y5D, size_t n, hoNDArray<double>& r, bool copyY2R);
    template EXPORTCPUCOREMATH bool multiplyOver4thDimensionExcept(const hoNDArray<GT_Complex8>& x4D, const hoNDArray<GT_Complex8>& y5D, size_t n, hoNDArray<GT_Complex8>& r, bool copyY2R);
    template EXPORTCPUCOREMATH bool multiplyOver4thDimensionExcept(const hoNDArray<GT_Complex16>& x4D, const hoNDArray<GT_Complex16>& y5D, size_t n, hoNDArray<GT_Complex16>& r, bool copyY2R);

    template EXPORTCPUCOREMATH bool multiplyOver5thDimension(const hoNDArray<float>& x, const hoNDArray<float>& y, hoNDArray<float>& r);
    template EXPORTCPUCOREMATH bool multiplyOver5thDimension(const hoNDArray<double>& x, const hoNDArray<double>& y, hoNDArray<double>& r);
    template EXPORTCPUCOREMATH bool multiplyOver5thDimension(const hoNDArray<GT_Complex8>& x, const hoNDArray<GT_Complex8>& y, hoNDArray<GT_Complex8>& r);
    template EXPORTCPUCOREMATH bool multiplyOver5thDimension(const hoNDArray<GT_Complex16>& x, const hoNDArray<GT_Complex16>& y, hoNDArray<GT_Complex16>& r);

    template EXPORTCPUCOREMATH bool multiplyOver5thDimensionExcept(const hoNDArray<float>& x, const hoNDArray<float>& y, size_t n, hoNDArray<float>& r, bool copyY2R);
    template EXPORTCPUCOREMATH bool multiplyOver5thDimensionExcept(const hoNDArray<double>& x, const hoNDArray<double>& y, size_t n, hoNDArray<double>& r, bool copyY2R);
    template EXPORTCPUCOREMATH bool multiplyOver5thDimensionExcept(const hoNDArray<GT_Complex8>& x, const hoNDArray<GT_Complex8>& y, size_t n, hoNDArray<GT_Complex8>& r, bool copyY2R);
    template EXPORTCPUCOREMATH bool multiplyOver5thDimensionExcept(const hoNDArray<GT_Complex16>& x, const hoNDArray<GT_Complex16>& y, size_t n, hoNDArray<GT_Complex16>& r, bool copyY2R);

    template EXPORTCPUCOREMATH bool multipleAdd(const hoNDArray<float>& x, const hoNDArray<float>& y, hoNDArray<float>& r);
    template EXPORTCPUCOREMATH bool multipleAdd(const hoNDArray<double>& x, const hoNDArray<double>& y, hoNDArray<double>& r);
    template EXPORTCPUCOREMATH bool multipleAdd(const hoNDArray<GT_Complex8>& x, const hoNDArray<GT_Complex8>& y, hoNDArray<GT_Complex8>& r);
    template EXPORTCPUCOREMATH bool multipleAdd(const hoNDArray<GT_Complex16>& x, const hoNDArray<GT_Complex16>& y, hoNDArray<GT_Complex16>& r);

    template EXPORTCPUCOREMATH bool multipleMultiply(const hoNDArray<float>& x, const hoNDArray<float>& y, hoNDArray<float>& r);
    template EXPORTCPUCOREMATH bool multipleMultiply(const hoNDArray<double>& x, const hoNDArray<double>& y, hoNDArray<double>& r);
    template EXPORTCPUCOREMATH bool multipleMultiply(const hoNDArray<GT_Complex8>& x, const hoNDArray<GT_Complex8>& y, hoNDArray<GT_Complex8>& r);
    template EXPORTCPUCOREMATH bool multipleMultiply(const hoNDArray<GT_Complex16>& x, const hoNDArray<GT_Complex16>& y, hoNDArray<GT_Complex16>& r);

    template EXPORTCPUCOREMATH bool cropUpTo10DArray(const hoNDArray<short>& x, hoNDArray<short>& r, const std::vector<size_t>& start, std::vector<size_t>& size);
    template EXPORTCPUCOREMATH bool cropUpTo10DArray(const hoNDArray<unsigned short>& x, hoNDArray<unsigned short>& r, const std::vector<size_t>& start, std::vector<size_t>& size);
    template EXPORTCPUCOREMATH bool cropUpTo10DArray(const hoNDArray<float>& x, hoNDArray<float>& r, const std::vector<size_t>& start, std::vector<size_t>& size);
    template EXPORTCPUCOREMATH bool cropUpTo10DArray(const hoNDArray<double>& x, hoNDArray<double>& r, const std::vector<size_t>& start, std::vector<size_t>& size);
    template EXPORTCPUCOREMATH bool cropUpTo10DArray(const hoNDArray<GT_Complex8>& x, hoNDArray<GT_Complex8>& r, const std::vector<size_t>& start, std::vector<size_t>& size);
    template EXPORTCPUCOREMATH bool cropUpTo10DArray(const hoNDArray<GT_Complex16>& x, hoNDArray<GT_Complex16>& r, const std::vector<size_t>& start, std::vector<size_t>& size);

    template EXPORTCPUCOREMATH bool setSubArrayUpTo10DArray(const hoNDArray<short>& x, hoNDArray<short>& r, const std::vector<size_t>& start, std::vector<size_t>& size);
    template EXPORTCPUCOREMATH bool setSubArrayUpTo10DArray(const hoNDArray<unsigned short>& x, hoNDArray<unsigned short>& r, const std::vector<size_t>& start, std::vector<size_t>& size);
    template EXPORTCPUCOREMATH bool setSubArrayUpTo10DArray(const hoNDArray<float>& x, hoNDArray<float>& r, const std::vector<size_t>& start, std::vector<size_t>& size);
    template EXPORTCPUCOREMATH bool setSubArrayUpTo10DArray(const hoNDArray<double>& x, hoNDArray<double>& r, const std::vector<size_t>& start, std::vector<size_t>& size);
    template EXPORTCPUCOREMATH bool setSubArrayUpTo10DArray(const hoNDArray<GT_Complex8>& x, hoNDArray<GT_Complex8>& r, const std::vector<size_t>& start, std::vector<size_t>& size);
    template EXPORTCPUCOREMATH bool setSubArrayUpTo10DArray(const hoNDArray<GT_Complex16>& x, hoNDArray<GT_Complex16>& r, const std::vector<size_t>& start, std::vector<size_t>& size);

    template EXPORTCPUCOREMATH bool cropOver3rdDimension(const hoNDArray<short>& x, hoNDArray<short>& r, size_t start, size_t end);
    template EXPORTCPUCOREMATH bool cropOver3rdDimension(const hoNDArray<unsigned short>& x, hoNDArray<unsigned short>& r, size_t start, size_t end);
    template EXPORTCPUCOREMATH bool cropOver3rdDimension(const hoNDArray<float>& x, hoNDArray<float>& r, size_t start, size_t end);
    template EXPORTCPUCOREMATH bool cropOver3rdDimension(const hoNDArray<double>& x, hoNDArray<double>& r, size_t start, size_t end);
    template EXPORTCPUCOREMATH bool cropOver3rdDimension(const hoNDArray<GT_Complex8>& x, hoNDArray<GT_Complex8>& r, size_t start, size_t end);
    template EXPORTCPUCOREMATH bool cropOver3rdDimension(const hoNDArray<GT_Complex16>& x, hoNDArray<GT_Complex16>& r, size_t start, size_t end);

    template EXPORTCPUCOREMATH bool setSubArrayOver3rdDimension(const hoNDArray<short>& x, hoNDArray<short>& r, size_t start, size_t end);
    template EXPORTCPUCOREMATH bool setSubArrayOver3rdDimension(const hoNDArray<unsigned short>& x, hoNDArray<unsigned short>& r, size_t start, size_t end);
    template EXPORTCPUCOREMATH bool setSubArrayOver3rdDimension(const hoNDArray<float>& x, hoNDArray<float>& r, size_t start, size_t end);
    template EXPORTCPUCOREMATH bool setSubArrayOver3rdDimension(const hoNDArray<double>& x, hoNDArray<double>& r, size_t start, size_t end);
    template EXPORTCPUCOREMATH bool setSubArrayOver3rdDimension(const hoNDArray<GT_Complex8>& x, hoNDArray<GT_Complex8>& r, size_t start, size_t end);
    template EXPORTCPUCOREMATH bool setSubArrayOver3rdDimension(const hoNDArray<GT_Complex16>& x, hoNDArray<GT_Complex16>& r, size_t start, size_t end);

    template EXPORTCPUCOREMATH bool stdOver3rdDimension(const hoNDArray<float>& x, hoNDArray<float>& std, bool NMinusOne);
    template EXPORTCPUCOREMATH bool stdOver3rdDimension(const hoNDArray<double>& x, hoNDArray<double>& std, bool NMinusOne);
    template EXPORTCPUCOREMATH bool stdOver3rdDimension(const hoNDArray<GT_Complex8>& x, hoNDArray<GT_Complex8>& std, bool NMinusOne);
    template EXPORTCPUCOREMATH bool stdOver3rdDimension(const hoNDArray<GT_Complex16>& x, hoNDArray<GT_Complex16>& std, bool NMinusOne);

    //template EXPORTCPUCOREMATH bool permuteLastTwoDimensions(const hoNDArray<float>& x, hoNDArray<float>& r);
    //template EXPORTCPUCOREMATH bool permuteLastTwoDimensions(const hoNDArray<double>& x, hoNDArray<double>& r);
    //template EXPORTCPUCOREMATH bool permuteLastTwoDimensions(const hoNDArray<GT_Complex8>& x, hoNDArray<GT_Complex8>& r);
    //template EXPORTCPUCOREMATH bool permuteLastTwoDimensions(const hoNDArray<GT_Complex16>& x, hoNDArray<GT_Complex16>& r);

    template EXPORTCPUCOREMATH bool permuteE2To3rdDimension(const hoNDArray<float>& x, hoNDArray<float>& r);
    template EXPORTCPUCOREMATH bool permuteE2To3rdDimension(const hoNDArray<double>& x, hoNDArray<double>& r);
    template EXPORTCPUCOREMATH bool permuteE2To3rdDimension(const hoNDArray<GT_Complex8>& x, hoNDArray<GT_Complex8>& r);
    template EXPORTCPUCOREMATH bool permuteE2To3rdDimension(const hoNDArray<GT_Complex16>& x, hoNDArray<GT_Complex16>& r);

    template EXPORTCPUCOREMATH bool permuteE2To5thDimension(const hoNDArray<float>& x, hoNDArray<float>& r);
    template EXPORTCPUCOREMATH bool permuteE2To5thDimension(const hoNDArray<double>& x, hoNDArray<double>& r);
    template EXPORTCPUCOREMATH bool permuteE2To5thDimension(const hoNDArray<GT_Complex8>& x, hoNDArray<GT_Complex8>& r);
    template EXPORTCPUCOREMATH bool permuteE2To5thDimension(const hoNDArray<GT_Complex16>& x, hoNDArray<GT_Complex16>& r);

    template EXPORTCPUCOREMATH bool permuteROTo3rdDimensionFor3DRecon(const hoNDArray<float>& x, hoNDArray<float>& r);
    template EXPORTCPUCOREMATH bool permuteROTo3rdDimensionFor3DRecon(const hoNDArray<double>& x, hoNDArray<double>& r);
    template EXPORTCPUCOREMATH bool permuteROTo3rdDimensionFor3DRecon(const hoNDArray<GT_Complex8>& x, hoNDArray<GT_Complex8>& r);
    template EXPORTCPUCOREMATH bool permuteROTo3rdDimensionFor3DRecon(const hoNDArray<GT_Complex16>& x, hoNDArray<GT_Complex16>& r);

    template EXPORTCPUCOREMATH bool permuteROTo4thDimensionFor3DRecon(const hoNDArray<float>& x, hoNDArray<float>& r);
    template EXPORTCPUCOREMATH bool permuteROTo4thDimensionFor3DRecon(const hoNDArray<double>& x, hoNDArray<double>& r);
    template EXPORTCPUCOREMATH bool permuteROTo4thDimensionFor3DRecon(const hoNDArray<GT_Complex8>& x, hoNDArray<GT_Complex8>& r);
    template EXPORTCPUCOREMATH bool permuteROTo4thDimensionFor3DRecon(const hoNDArray<GT_Complex16>& x, hoNDArray<GT_Complex16>& r);

    template EXPORTCPUCOREMATH bool permuteROTo1stDimensionFor3DRecon(const hoNDArray<float>& x, hoNDArray<float>& r);
    template EXPORTCPUCOREMATH bool permuteROTo1stDimensionFor3DRecon(const hoNDArray<double>& x, hoNDArray<double>& r);
    template EXPORTCPUCOREMATH bool permuteROTo1stDimensionFor3DRecon(const hoNDArray<GT_Complex8>& x, hoNDArray<GT_Complex8>& r);
    template EXPORTCPUCOREMATH bool permuteROTo1stDimensionFor3DRecon(const hoNDArray<GT_Complex16>& x, hoNDArray<GT_Complex16>& r);

    template EXPORTCPUCOREMATH bool permute3rdDimensionTo1stDimension(const hoNDArray<float>& x, hoNDArray<float>& r);
    template EXPORTCPUCOREMATH bool permute3rdDimensionTo1stDimension(const hoNDArray<double>& x, hoNDArray<double>& r);
    template EXPORTCPUCOREMATH bool permute3rdDimensionTo1stDimension(const hoNDArray<GT_Complex8>& x, hoNDArray<GT_Complex8>& r);
    template EXPORTCPUCOREMATH bool permute3rdDimensionTo1stDimension(const hoNDArray<GT_Complex16>& x, hoNDArray<GT_Complex16>& r);

    template EXPORTCPUCOREMATH bool permuteROTo5thDimensionFor3DRecon(const hoNDArray<float>& x, hoNDArray<float>& r);
    template EXPORTCPUCOREMATH bool permuteROTo5thDimensionFor3DRecon(const hoNDArray<double>& x, hoNDArray<double>& r);
    template EXPORTCPUCOREMATH bool permuteROTo5thDimensionFor3DRecon(const hoNDArray<GT_Complex8>& x, hoNDArray<GT_Complex8>& r);
    template EXPORTCPUCOREMATH bool permuteROTo5thDimensionFor3DRecon(const hoNDArray<GT_Complex16>& x, hoNDArray<GT_Complex16>& r);

    template EXPORTCPUCOREMATH bool imageDomainUnwrapping2D(const hoNDArray<GT_Complex8>& x, const hoNDArray<GT_Complex8>& ker, hoNDArray<GT_Complex8>& buf, hoNDArray<GT_Complex8>& y);
    template EXPORTCPUCOREMATH bool imageDomainUnwrapping2D(const hoNDArray<GT_Complex16>& x, const hoNDArray<GT_Complex16>& ker, hoNDArray<GT_Complex16>& buf, hoNDArray<GT_Complex16>& y);

    template EXPORTCPUCOREMATH bool imageDomainUnwrapping2DT(const hoNDArray<GT_Complex8>& x, const hoNDArray<GT_Complex8>& ker, hoNDArray<GT_Complex8>& buf, hoNDArray<GT_Complex8>& y);
    template EXPORTCPUCOREMATH bool imageDomainUnwrapping2DT(const hoNDArray<GT_Complex16>& x, const hoNDArray<GT_Complex16>& ker, hoNDArray<GT_Complex16>& buf, hoNDArray<GT_Complex16>& y);

    #endif // USE_MKL

    //
    // Instantiation
    //

    template EXPORTCPUCOREMATH void clear<short>( hoNDArray<short>& );
    template EXPORTCPUCOREMATH void clear<unsigned short>( hoNDArray<unsigned short>& );
    template EXPORTCPUCOREMATH void clear<int>( hoNDArray<int>& );
    template EXPORTCPUCOREMATH void clear<size_t>( hoNDArray<size_t>& );

    template EXPORTCPUCOREMATH boost::shared_ptr< hoNDArray<float> > abs<float>( hoNDArray<float>* );
    template EXPORTCPUCOREMATH void abs_inplace<float>( hoNDArray<float>* );
    template EXPORTCPUCOREMATH boost::shared_ptr< hoNDArray<float> > abs_square<float>( hoNDArray<float>* );
    template EXPORTCPUCOREMATH boost::shared_ptr< hoNDArray<float> > sqrt<float>( hoNDArray<float>* );
    template EXPORTCPUCOREMATH void sqrt_inplace<float>( hoNDArray<float>* );
    template EXPORTCPUCOREMATH boost::shared_ptr< hoNDArray<float> > square<float>( hoNDArray<float>* );
    template EXPORTCPUCOREMATH void square_inplace<float>( hoNDArray<float>* );
    template EXPORTCPUCOREMATH boost::shared_ptr< hoNDArray<float> > reciprocal<float>( hoNDArray<float>* );
    template EXPORTCPUCOREMATH void reciprocal_inplace<float>( hoNDArray<float>* );
    template EXPORTCPUCOREMATH boost::shared_ptr< hoNDArray<float> > reciprocal_sqrt<float>( hoNDArray<float>* );
    template EXPORTCPUCOREMATH void reciprocal_sqrt_inplace<float>( hoNDArray<float>* );
    template EXPORTCPUCOREMATH boost::shared_ptr< hoNDArray<float> > sgn<float>( hoNDArray<float>* );
    template EXPORTCPUCOREMATH void sgn_inplace<float>( hoNDArray<float>* );
    template EXPORTCPUCOREMATH void clear<float>( hoNDArray<float>* );
    template EXPORTCPUCOREMATH void clear<float>( hoNDArray<float>& );
    template EXPORTCPUCOREMATH void fill<float>( hoNDArray<float>*, float );
    template EXPORTCPUCOREMATH void clamp<float>( hoNDArray<float>*, float, float );
    template EXPORTCPUCOREMATH void clamp_min<float>( hoNDArray<float>*, float );
    template EXPORTCPUCOREMATH void clamp_max<float>( hoNDArray<float>*, float );
    template EXPORTCPUCOREMATH void normalize<float>( hoNDArray<float>*, float );
    template EXPORTCPUCOREMATH void shrink1<float>( hoNDArray<float>*, float, hoNDArray<float>* );
    template EXPORTCPUCOREMATH void pshrink<float>( hoNDArray<float>*, float,float, hoNDArray<float>* );
    template EXPORTCPUCOREMATH void shrinkd<float> ( hoNDArray<float>*, hoNDArray<float>*, float, hoNDArray<float>* );
    template EXPORTCPUCOREMATH void pshrinkd<float> ( hoNDArray<float>*, hoNDArray<float>*, float, float, hoNDArray<float>* );

    template EXPORTCPUCOREMATH boost::shared_ptr< hoNDArray<double> > abs<double>( hoNDArray<double>* );
    template EXPORTCPUCOREMATH void abs_inplace<double>( hoNDArray<double>* );
    template EXPORTCPUCOREMATH boost::shared_ptr< hoNDArray<double> > abs_square<double>( hoNDArray<double>* );
    template EXPORTCPUCOREMATH boost::shared_ptr< hoNDArray<double> > sqrt<double>( hoNDArray<double>* );
    template EXPORTCPUCOREMATH void sqrt_inplace<double>( hoNDArray<double>* );
    template EXPORTCPUCOREMATH boost::shared_ptr< hoNDArray<double> > square<double>( hoNDArray<double>* );
    template EXPORTCPUCOREMATH void square_inplace<double>( hoNDArray<double>* );
    template EXPORTCPUCOREMATH boost::shared_ptr< hoNDArray<double> > reciprocal<double>( hoNDArray<double>* );
    template EXPORTCPUCOREMATH void reciprocal_inplace<double>( hoNDArray<double>* );
    template EXPORTCPUCOREMATH boost::shared_ptr< hoNDArray<double> > reciprocal_sqrt<double>( hoNDArray<double>* );
    template EXPORTCPUCOREMATH void reciprocal_sqrt_inplace<double>( hoNDArray<double>* );
    template EXPORTCPUCOREMATH boost::shared_ptr< hoNDArray<double> > sgn<double>( hoNDArray<double>* );
    template EXPORTCPUCOREMATH void sgn_inplace<double>( hoNDArray<double>* );
    template EXPORTCPUCOREMATH void clear<double>( hoNDArray<double>* );
    template EXPORTCPUCOREMATH void clear<double>( hoNDArray<double>& );
    template EXPORTCPUCOREMATH void fill<double>( hoNDArray<double>*, double );
    template EXPORTCPUCOREMATH void clamp<double>( hoNDArray<double>*, double, double );
    template EXPORTCPUCOREMATH void clamp_min<double>( hoNDArray<double>*, double );
    template EXPORTCPUCOREMATH void clamp_max<double>( hoNDArray<double>*, double );
    template EXPORTCPUCOREMATH void normalize<double>( hoNDArray<double>*, double );
    template EXPORTCPUCOREMATH void shrink1<double>( hoNDArray<double>*, double, hoNDArray<double>* );
    template EXPORTCPUCOREMATH void pshrink<double>( hoNDArray<double>*, double,double, hoNDArray<double>* );
    template EXPORTCPUCOREMATH void shrinkd<double> ( hoNDArray<double>*, hoNDArray<double>*, double, hoNDArray<double>* );
    template EXPORTCPUCOREMATH void pshrinkd<double> ( hoNDArray<double>*, hoNDArray<double>*, double, double, hoNDArray<double>* );

    template EXPORTCPUCOREMATH boost::shared_ptr< hoNDArray<float> > abs< std::complex<float> >( hoNDArray< std::complex<float> >* );
    template EXPORTCPUCOREMATH boost::shared_ptr< hoNDArray<float> > abs_square< std::complex<float> >( hoNDArray< std::complex<float> >* );
    template EXPORTCPUCOREMATH boost::shared_ptr< hoNDArray< std::complex<float> > > sqrt< std::complex<float> >( hoNDArray< std::complex<float> >* );
    template EXPORTCPUCOREMATH void sqrt_inplace< std::complex<float> >( hoNDArray< std::complex<float> >* );
    template EXPORTCPUCOREMATH boost::shared_ptr< hoNDArray< std::complex<float> > > square< std::complex<float> >( hoNDArray< std::complex<float> >* );
    template EXPORTCPUCOREMATH void square_inplace< std::complex<float> >( hoNDArray< std::complex<float> >* );
    template EXPORTCPUCOREMATH boost::shared_ptr< hoNDArray< std::complex<float> > > reciprocal< std::complex<float> >( hoNDArray< std::complex<float> >* );
    template EXPORTCPUCOREMATH void reciprocal_inplace< std::complex<float> >( hoNDArray< std::complex<float> >* );
    template EXPORTCPUCOREMATH boost::shared_ptr< hoNDArray< std::complex<float> > > reciprocal_sqrt< std::complex<float> >( hoNDArray< std::complex<float> >* );
    template EXPORTCPUCOREMATH void reciprocal_sqrt_inplace< std::complex<float> >( hoNDArray< std::complex<float> >* );
    template EXPORTCPUCOREMATH void clear< std::complex<float> >( hoNDArray< std::complex<float> >* );
    template EXPORTCPUCOREMATH void clear< std::complex<float> >( hoNDArray< std::complex<float> >& );
    template EXPORTCPUCOREMATH void fill< std::complex<float> >( hoNDArray< std::complex<float> >*, std::complex<float> );
    template EXPORTCPUCOREMATH void clamp< std::complex<float> >( hoNDArray< std::complex<float> >*, float, float );
    template EXPORTCPUCOREMATH void clamp_min< std::complex<float> >( hoNDArray< std::complex<float> >*, float );
    template EXPORTCPUCOREMATH void clamp_max<std::complex<float> >( hoNDArray< std::complex<float> >*, float );
    template EXPORTCPUCOREMATH void normalize< std::complex<float> >( hoNDArray< std::complex<float> >*, float );
    template EXPORTCPUCOREMATH void shrink1< std::complex<float> >( hoNDArray< std::complex<float> >*, float, hoNDArray< std::complex<float> >* );
    template EXPORTCPUCOREMATH void pshrink< std::complex<float> >( hoNDArray< std::complex<float> >*, float,float, hoNDArray< std::complex<float> >* );
    template EXPORTCPUCOREMATH void shrinkd< std::complex<float> > ( hoNDArray< std::complex<float> >*, hoNDArray<float>*, float, hoNDArray< std::complex<float> >* );
    template EXPORTCPUCOREMATH void pshrinkd< std::complex<float> > ( hoNDArray< std::complex<float> >*, hoNDArray<float>*, float, float, hoNDArray< std::complex<float> >* );

    template EXPORTCPUCOREMATH boost::shared_ptr< hoNDArray<double> > abs< std::complex<double> >( hoNDArray< std::complex<double> >* );
    template EXPORTCPUCOREMATH boost::shared_ptr< hoNDArray<double> > abs_square< std::complex<double> >( hoNDArray< std::complex<double> >* );
    template EXPORTCPUCOREMATH boost::shared_ptr< hoNDArray< std::complex<double> > > sqrt< std::complex<double> >( hoNDArray< std::complex<double> >* );
    template EXPORTCPUCOREMATH void sqrt_inplace< std::complex<double> >( hoNDArray< std::complex<double> >* );
    template EXPORTCPUCOREMATH boost::shared_ptr< hoNDArray< std::complex<double> > > square< std::complex<double> >( hoNDArray< std::complex<double> >* );
    template EXPORTCPUCOREMATH void square_inplace< std::complex<double> >( hoNDArray< std::complex<double> >* );
    template EXPORTCPUCOREMATH boost::shared_ptr< hoNDArray< std::complex<double> > > reciprocal< std::complex<double> >( hoNDArray< std::complex<double> >* );
    template EXPORTCPUCOREMATH void reciprocal_inplace< std::complex<double> >( hoNDArray< std::complex<double> >* );
    template EXPORTCPUCOREMATH boost::shared_ptr< hoNDArray< std::complex<double> > > reciprocal_sqrt< std::complex<double> >( hoNDArray< std::complex<double> >* );
    template EXPORTCPUCOREMATH void reciprocal_sqrt_inplace< std::complex<double> >( hoNDArray< std::complex<double> >* );
    template EXPORTCPUCOREMATH void clear< std::complex<double> >( hoNDArray< std::complex<double> >* );
    template EXPORTCPUCOREMATH void clear< std::complex<double> >( hoNDArray< std::complex<double> >& );
    template EXPORTCPUCOREMATH void fill< std::complex<double> >( hoNDArray< std::complex<double> >*, std::complex<double> );
    template EXPORTCPUCOREMATH void clamp< std::complex<double> >( hoNDArray< std::complex<double> >*, double, double );
    template EXPORTCPUCOREMATH void clamp_min< std::complex<double> >( hoNDArray< std::complex<double> >*, double );
    template EXPORTCPUCOREMATH void clamp_max<std::complex<double> >( hoNDArray< std::complex<double> >*, double );
    template EXPORTCPUCOREMATH void normalize< std::complex<double> >( hoNDArray< std::complex<double> >*, double );
    template EXPORTCPUCOREMATH void shrink1< std::complex<double> >( hoNDArray< std::complex<double> >*, double, hoNDArray< std::complex<double> >* );
    template EXPORTCPUCOREMATH void pshrink< std::complex<double> >( hoNDArray< std::complex<double> >*, double,double, hoNDArray< std::complex<double> >* );
    template EXPORTCPUCOREMATH void shrinkd< std::complex<double> > ( hoNDArray< std::complex<double> >*, hoNDArray<double>*, double, hoNDArray< std::complex<double> >* );
    template EXPORTCPUCOREMATH void pshrinkd< std::complex<double> > ( hoNDArray< std::complex<double> >*, hoNDArray<double>*, double, double, hoNDArray< std::complex<double> >* );

    template EXPORTCPUCOREMATH boost::shared_ptr< hoNDArray<float> > abs< complext<float> >( hoNDArray< complext<float> >* );
    template EXPORTCPUCOREMATH boost::shared_ptr< hoNDArray<float> > abs_square< complext<float> >( hoNDArray< complext<float> >* );
    template EXPORTCPUCOREMATH boost::shared_ptr< hoNDArray< complext<float> > > sqrt< complext<float> >( hoNDArray< complext<float> >* );
    template EXPORTCPUCOREMATH void sqrt_inplace< complext<float> >( hoNDArray< complext<float> >* );
    template EXPORTCPUCOREMATH boost::shared_ptr< hoNDArray< complext<float> > > square< complext<float> >( hoNDArray< complext<float> >* );
    template EXPORTCPUCOREMATH void square_inplace< complext<float> >( hoNDArray< complext<float> >* );
    template EXPORTCPUCOREMATH boost::shared_ptr< hoNDArray< complext<float> > > reciprocal< complext<float> >( hoNDArray< complext<float> >* );
    template EXPORTCPUCOREMATH void reciprocal_inplace< complext<float> >( hoNDArray< complext<float> >* );
    template EXPORTCPUCOREMATH boost::shared_ptr< hoNDArray< complext<float> > > reciprocal_sqrt< complext<float> >( hoNDArray< complext<float> >* );
    template EXPORTCPUCOREMATH void reciprocal_sqrt_inplace< complext<float> >( hoNDArray< complext<float> >* );
    template EXPORTCPUCOREMATH void clear< complext<float> >( hoNDArray< complext<float> >* );
    template EXPORTCPUCOREMATH void clear< complext<float> >( hoNDArray< complext<float> >& );
    template EXPORTCPUCOREMATH void fill< complext<float> >( hoNDArray< complext<float> >*, complext<float> );
    template EXPORTCPUCOREMATH void clamp< complext<float> >( hoNDArray< complext<float> >*, float, float );
    template EXPORTCPUCOREMATH void clamp_min< complext<float> >( hoNDArray< complext<float> >*, float );
    template EXPORTCPUCOREMATH void clamp_max<complext<float> >( hoNDArray< complext<float> >*, float );
    template EXPORTCPUCOREMATH void normalize< complext<float> >( hoNDArray< complext<float> >*, float );
    template EXPORTCPUCOREMATH void shrink1< complext<float> >( hoNDArray< complext<float> >*, float, hoNDArray< complext<float> >* );
    template EXPORTCPUCOREMATH void pshrink< complext<float> >( hoNDArray< complext<float> >*, float,float, hoNDArray< complext<float> >* );
    template EXPORTCPUCOREMATH void shrinkd< complext<float> > ( hoNDArray< complext<float> >*, hoNDArray<float>*, float, hoNDArray< complext<float> >* );
    template EXPORTCPUCOREMATH void pshrinkd< complext<float> > ( hoNDArray< complext<float> >*, hoNDArray<float>*, float, float, hoNDArray< complext<float> >* );

    template EXPORTCPUCOREMATH boost::shared_ptr< hoNDArray<double> > abs< complext<double> >( hoNDArray< complext<double> >* );
    template EXPORTCPUCOREMATH boost::shared_ptr< hoNDArray<double> > abs_square< complext<double> >( hoNDArray< complext<double> >* );
    template EXPORTCPUCOREMATH boost::shared_ptr< hoNDArray< complext<double> > > sqrt< complext<double> >( hoNDArray< complext<double> >* );
    template EXPORTCPUCOREMATH void sqrt_inplace< complext<double> >( hoNDArray< complext<double> >* );
    template EXPORTCPUCOREMATH boost::shared_ptr< hoNDArray< complext<double> > > square< complext<double> >( hoNDArray< complext<double> >* );
    template EXPORTCPUCOREMATH void square_inplace< complext<double> >( hoNDArray< complext<double> >* );
    template EXPORTCPUCOREMATH boost::shared_ptr< hoNDArray< complext<double> > > reciprocal< complext<double> >( hoNDArray< complext<double> >* );
    template EXPORTCPUCOREMATH void reciprocal_inplace< complext<double> >( hoNDArray< complext<double> >* );
    template EXPORTCPUCOREMATH boost::shared_ptr< hoNDArray< complext<double> > > reciprocal_sqrt< complext<double> >( hoNDArray< complext<double> >* );
    template EXPORTCPUCOREMATH void reciprocal_sqrt_inplace< complext<double> >( hoNDArray< complext<double> >* );
    template EXPORTCPUCOREMATH void clear< complext<double> >( hoNDArray< complext<double> >* );
    template EXPORTCPUCOREMATH void clear< complext<double> >( hoNDArray< complext<double> >& );
    template EXPORTCPUCOREMATH void fill< complext<double> >( hoNDArray< complext<double> >*, complext<double> );
    template EXPORTCPUCOREMATH void clamp< complext<double> >( hoNDArray< complext<double> >*, double, double );
    template EXPORTCPUCOREMATH void clamp_min< complext<double> >( hoNDArray< complext<double> >*, double );
    template EXPORTCPUCOREMATH void clamp_max<complext<double> >( hoNDArray< complext<double> >*, double );
    template EXPORTCPUCOREMATH void normalize< complext<double> >( hoNDArray< complext<double> >*, double );
    template EXPORTCPUCOREMATH void shrink1< complext<double> >( hoNDArray< complext<double> >*, double, hoNDArray< complext<double> >* );
    template EXPORTCPUCOREMATH void pshrink< complext<double> >( hoNDArray< complext<double> >*, double,double, hoNDArray< complext<double> >* );
    template EXPORTCPUCOREMATH void shrinkd< complext<double> > ( hoNDArray< complext<double> >*, hoNDArray<double>*, double, hoNDArray< complext<double> >* );
    template EXPORTCPUCOREMATH void pshrinkd< complext<double> > ( hoNDArray< complext<double> >*, hoNDArray<double>*, double, double, hoNDArray< complext<double> >* );

    template EXPORTCPUCOREMATH boost::shared_ptr< hoNDArray< std::complex<float> > > real_to_complex< std::complex<float> >( hoNDArray<float>* );
    template EXPORTCPUCOREMATH boost::shared_ptr< hoNDArray< std::complex<float> > > real_imag_to_complex< std::complex<float> >( hoNDArray<float>*, hoNDArray<float>* );

    template EXPORTCPUCOREMATH bool real_imag_to_complex(const hoNDArray<float>& real, const hoNDArray<float>& imag, hoNDArray< std::complex<float> >& cplx);
    template EXPORTCPUCOREMATH bool real_imag_to_complex(const hoNDArray<float>& real, const hoNDArray<float>& imag, hoNDArray< float_complext >& cplx);

    template EXPORTCPUCOREMATH bool complex_to_real_imag(const hoNDArray< std::complex<float> >& cplx, hoNDArray<float>& real, hoNDArray<float>& imag);
    //template EXPORTCPUCOREMATH bool complex_to_real_imag(const hoNDArray< float_complext >& cplx, hoNDArray<float>& real, hoNDArray<float>& imag);

    template EXPORTCPUCOREMATH bool complex_to_real(const hoNDArray< std::complex<float> >& cplx, hoNDArray<float>& real);
    template EXPORTCPUCOREMATH bool complex_to_imag(const hoNDArray< std::complex<float> >& cplx, hoNDArray<float>& imag);

    template EXPORTCPUCOREMATH boost::shared_ptr< hoNDArray<float_complext> > real_to_complex<float_complext>( hoNDArray<float>* );
    template EXPORTCPUCOREMATH boost::shared_ptr< hoNDArray<float_complext> > real_imag_to_complex<float_complext>( hoNDArray<float>*, hoNDArray<float>* );
    template EXPORTCPUCOREMATH boost::shared_ptr< hoNDArray<float> > real<float>( hoNDArray<float>* );
    template EXPORTCPUCOREMATH boost::shared_ptr< hoNDArray<float> > real<std::complex<float> >( hoNDArray< std::complex<float> >* );
    template EXPORTCPUCOREMATH boost::shared_ptr< hoNDArray<float> > real<float_complext>( hoNDArray<float_complext>* );
    template EXPORTCPUCOREMATH boost::shared_ptr< hoNDArray<float> > imag<float>( hoNDArray<float>* );
    template EXPORTCPUCOREMATH boost::shared_ptr< hoNDArray<float> > imag<std::complex<float> >( hoNDArray< std::complex<float> >* );
    template EXPORTCPUCOREMATH boost::shared_ptr< hoNDArray<float> > imag<float_complext>( hoNDArray<float_complext>* );

    template EXPORTCPUCOREMATH boost::shared_ptr< hoNDArray<float> > conj<float>( hoNDArray<float>* );
    template EXPORTCPUCOREMATH boost::shared_ptr< hoNDArray<std::complex<float> > > conj<std::complex<float> >( hoNDArray<std::complex<float> >* );
    template EXPORTCPUCOREMATH boost::shared_ptr< hoNDArray<float_complext> > conj<float_complext>( hoNDArray<float_complext>* );


    template EXPORTCPUCOREMATH boost::shared_ptr< hoNDArray< std::complex<double> > > real_to_complex< std::complex<double> >( hoNDArray<double>* );
    template EXPORTCPUCOREMATH boost::shared_ptr< hoNDArray< std::complex<double> > > real_imag_to_complex< std::complex<double> >( hoNDArray<double>*, hoNDArray<double>* );

    template EXPORTCPUCOREMATH bool real_imag_to_complex(const hoNDArray<double>& real, const hoNDArray<double>& imag, hoNDArray< std::complex<double> >& cplx);
    template EXPORTCPUCOREMATH bool real_imag_to_complex(const hoNDArray<double>& real, const hoNDArray<double>& imag, hoNDArray< double_complext >& cplx);

    template EXPORTCPUCOREMATH bool complex_to_real_imag(const hoNDArray< std::complex<double> >& cplx, hoNDArray<double>& real, hoNDArray<double>& imag);
    //template EXPORTCPUCOREMATH bool complex_to_real_imag(const hoNDArray< double >& cplx, hoNDArray<double>& real, hoNDArray<double>& imag);
    //template EXPORTCPUCOREMATH bool complex_to_real_imag(const hoNDArray< float >& cplx, hoNDArray<float>& real, hoNDArray<float>& imag);
    // template EXPORTCPUCOREMATH bool complex_to_real_imag(const hoNDArray< double_complext >& cplx, hoNDArray<double>& real, hoNDArray<double>& imag);

    template EXPORTCPUCOREMATH bool complex_to_real(const hoNDArray< std::complex<double> >& cplx, hoNDArray<double>& real);
    template EXPORTCPUCOREMATH bool complex_to_imag(const hoNDArray< std::complex<double> >& cplx, hoNDArray<double>& imag);

    template EXPORTCPUCOREMATH boost::shared_ptr< hoNDArray<double_complext> > real_to_complex<double_complext>( hoNDArray<double>* );
    template EXPORTCPUCOREMATH boost::shared_ptr< hoNDArray<double_complext> > real_imag_to_complex<double_complext>( hoNDArray<double>*, hoNDArray<double>* );
    template EXPORTCPUCOREMATH boost::shared_ptr< hoNDArray<double> > real<double>( hoNDArray<double>* );
    template EXPORTCPUCOREMATH boost::shared_ptr< hoNDArray<double> > real<std::complex<double> >( hoNDArray< std::complex<double> >* );
    template EXPORTCPUCOREMATH boost::shared_ptr< hoNDArray<double> > real<double_complext>( hoNDArray<double_complext>* );
    template EXPORTCPUCOREMATH boost::shared_ptr< hoNDArray<double> > imag<std::complex<double> >( hoNDArray< std::complex<double> >* );
    template EXPORTCPUCOREMATH boost::shared_ptr< hoNDArray<double> > imag<double>( hoNDArray<double>* );
    template EXPORTCPUCOREMATH boost::shared_ptr< hoNDArray<double> > imag<double_complext>( hoNDArray<double_complext>* );
    template EXPORTCPUCOREMATH boost::shared_ptr< hoNDArray<double> > conj<double>( hoNDArray<double>* );
    template EXPORTCPUCOREMATH boost::shared_ptr< hoNDArray<std::complex<double> > > conj<std::complex<double> >( hoNDArray<std::complex<double> >* );
    template EXPORTCPUCOREMATH boost::shared_ptr< hoNDArray<double_complext> > conj<double_complext>( hoNDArray<double_complext>* );
}
