#include "hoNDArray_elemwise.h"
#include "hoNDArray_operators.h"
#include "hoNDArray_blas.h"
#include "complext.h"
#include "hoArmadillo.h"

#include <complex>

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
        for( int i = 0; i < res->get_number_of_elements(); i++ ){
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
        for( int i = 0; i < x->get_number_of_elements(); i++ ) 
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
            throw std::runtime_error("Gadgetron::real_to_std_complex(): Invalid input array");

        boost::shared_ptr< hoNDArray<T> > result(new hoNDArray<T>());
        result->create(x->get_dimensions());
        arma::Col<typename stdType<T>::Type> aRes = as_arma_col(result.get());
        aRes = arma::Col<typename stdType<T>::Type>(as_arma_col(x), arma::Col<typename realType<T>::Type>(x->get_number_of_elements()).zeros());
        return result;
    }

    template<class T> void clear( hoNDArray<T> *x )
    {
        if( x == 0x0 )
            throw std::runtime_error("Gadgetron::clear(): Invalid input array");

        memset( x->get_data_ptr(), 0, x->get_number_of_elements()*sizeof(T));
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
        hoNDA_clamp( T _min, T _max ) : min(_min), max(_max) {}
        T operator()(const T &x) const 
        {
            if( x < min ) return min;
            else if ( x > max) return max;
            else return x;
        }
        T min, max;
    };

    template<typename T> struct hoNDA_clamp< std::complex<T> > //: public thrust::unary_function< std::complex<T>, std::complex<T> >
    {
        hoNDA_clamp( T _min, T _max ) : min(_min), max(_max) {}
        std::complex<T> operator()(const std::complex<T> &x) const 
        {
            if( real(x) < min ) return std::complex<T>(min);
            else if ( real(x) > max) return std::complex<T>(max);
            else return std::complex<T>(real(x));
        }
        T min, max;
    };

    template<typename T> struct hoNDA_clamp< complext<T> > //: public thrust::unary_function< complext<T>, complext<T> >
    {
        hoNDA_clamp( T _min, T _max ) : min(_min), max(_max) {}
        complext<T> operator()(const complext<T> &x) const 
        {
            if( real(x) < min ) return complext<T>(min);
            else if ( real(x) > max) return complext<T>(max);
            else return complext<T>(real(x));
        }
        T min, max;
    };

    template<class T> void clamp( hoNDArray<T> *x, typename realType<T>::Type min, typename realType<T>::Type max )
    { 
        if( x == 0x0 )
            throw std::runtime_error("Gadgetron::clamp(): Invalid input array");

        hoNDA_clamp<T> functor(min, max);
        std::transform(x->begin(),x->end(),x->begin(),functor);
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

        unsigned int max_idx = amax(x);
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
        for( int i = 0; i < x->get_number_of_elements(); i++ ) {
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
        for( int i = 0; i < x->get_number_of_elements(); i++ ) {
            T prev = x->get_data_ptr()[i];
            typename realType<T>::Type absPrev = abs(prev);
            T sgnPrev = (absPrev <= typename realType<T>::Type(0)) ? T(0) : prev/absPrev;
            outPtr[i] = sgnPrev*std::max(absPrev-gamma*std::pow(absPrev,p-1), typename realType<T>::Type(0));
        }
    }

    template<class T> void shrinkd ( hoNDArray<T> *_x, hoNDArray<typename realType<T>::Type> *_s, typename realType<T>::Type gamma, hoNDArray<T> *out )
    {
        if( _x == 0x0  || _s == 0 )
            throw std::runtime_error("Gadgetron::pshrinkd(): Invalid input array");

        T *outPtr = (out==0x0) ? _x->get_data_ptr() : out->get_data_ptr();

#ifdef USE_OMP
#pragma omp parallel for
#endif
        for( int i = 0; i < _x->get_number_of_elements(); i++ ) {
            T x = _x->get_data_ptr()[i];
            typename realType<T>::Type s = _s->get_data_ptr()[i];
            outPtr[i] = x/s*std::max(s-gamma,typename realType<T>::Type(0));
        } 
    }


    template<class T> void pshrinkd ( hoNDArray<T> *_x, hoNDArray<typename realType<T>::Type> *_s, typename realType<T>::Type gamma,typename realType<T>::Type p, hoNDArray<T> *out )
    {
        if( _x == 0x0  || _s == 0 )
            throw std::runtime_error("Gadgetron::shrinkd(): Invalid input array");

        T *outPtr = (out==0x0) ? _x->get_data_ptr() : out->get_data_ptr();

#ifdef USE_OMP
#pragma omp parallel for
#endif
        for( int i = 0; i < _x->get_number_of_elements(); i++ ) {
            T x = _x->get_data_ptr()[i];
            typename realType<T>::Type s = _s->get_data_ptr()[i];
            outPtr[i] = x/s*std::max(s-gamma*std::pow(s,p-1),typename realType<T>::Type(0));
        }
    }
    //
    // Instantiation
    //

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
    template EXPORTCPUCOREMATH void fill<float>( hoNDArray<float>*, float );
    template EXPORTCPUCOREMATH void clamp<float>( hoNDArray<float>*, float, float );
    template EXPORTCPUCOREMATH void clamp_min<float>( hoNDArray<float>*, float );
    template EXPORTCPUCOREMATH void clamp_max<float>( hoNDArray<float>*, float );
    template EXPORTCPUCOREMATH void normalize<float>( hoNDArray<float>*, float );
    template EXPORTCPUCOREMATH void shrink1<float>( hoNDArray<float>*, float, hoNDArray<float>* );
    template EXPORTCPUCOREMATH void pshrink<float>( hoNDArray<float>*, float,float, hoNDArray<float>* );
    template EXPORTCPUCOREMATH void shrinkd<float> ( hoNDArray<float>*, hoNDArray<float>*, float, hoNDArray<float>* );
    template EXPORTCPUCOREMATH void pshrinkd<float> ( hoNDArray<float>*, hoNDArray<float>*, float,float, hoNDArray<float>* );

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
    template EXPORTCPUCOREMATH void fill<double>( hoNDArray<double>*, double );
    template EXPORTCPUCOREMATH void clamp<double>( hoNDArray<double>*, double, double );
    template EXPORTCPUCOREMATH void clamp_min<double>( hoNDArray<double>*, double );
    template EXPORTCPUCOREMATH void clamp_max<double>( hoNDArray<double>*, double );
    template EXPORTCPUCOREMATH void normalize<double>( hoNDArray<double>*, double );
    template EXPORTCPUCOREMATH void shrink1<double>( hoNDArray<double>*, double, hoNDArray<double>* );
    template EXPORTCPUCOREMATH void pshrink<double>( hoNDArray<double>*, double,double, hoNDArray<double>* );
    template EXPORTCPUCOREMATH void shrinkd<double> ( hoNDArray<double>*, hoNDArray<double>*, double, hoNDArray<double>* );
    template EXPORTCPUCOREMATH void pshrinkd<double> ( hoNDArray<double>*, hoNDArray<double>*, double,double, hoNDArray<double>* );

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
    template EXPORTCPUCOREMATH void fill< std::complex<float> >( hoNDArray< std::complex<float> >*, std::complex<float> );
    template EXPORTCPUCOREMATH void clamp< std::complex<float> >( hoNDArray< std::complex<float> >*, float, float );
    template EXPORTCPUCOREMATH void clamp_min< std::complex<float> >( hoNDArray< std::complex<float> >*, float );
    template EXPORTCPUCOREMATH void clamp_max<std::complex<float> >( hoNDArray< std::complex<float> >*, float );
    template EXPORTCPUCOREMATH void normalize< std::complex<float> >( hoNDArray< std::complex<float> >*, float );
    template EXPORTCPUCOREMATH void shrink1< std::complex<float> >( hoNDArray< std::complex<float> >*, float, hoNDArray< std::complex<float> >* );
    template EXPORTCPUCOREMATH void pshrink< std::complex<float> >( hoNDArray< std::complex<float> >*, float,float, hoNDArray< std::complex<float> >* );
    template EXPORTCPUCOREMATH void shrinkd< std::complex<float> > ( hoNDArray< std::complex<float> >*, hoNDArray<float>*, float, hoNDArray< std::complex<float> >* );
    template EXPORTCPUCOREMATH void pshrinkd< std::complex<float> > ( hoNDArray< std::complex<float> >*, hoNDArray<float>*, float,float, hoNDArray< std::complex<float> >* );

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
    template EXPORTCPUCOREMATH void fill< std::complex<double> >( hoNDArray< std::complex<double> >*, std::complex<double> );
    template EXPORTCPUCOREMATH void clamp< std::complex<double> >( hoNDArray< std::complex<double> >*, double, double );
    template EXPORTCPUCOREMATH void clamp_min< std::complex<double> >( hoNDArray< std::complex<double> >*, double );
    template EXPORTCPUCOREMATH void clamp_max<std::complex<double> >( hoNDArray< std::complex<double> >*, double );
    template EXPORTCPUCOREMATH void normalize< std::complex<double> >( hoNDArray< std::complex<double> >*, double );
    template EXPORTCPUCOREMATH void shrink1< std::complex<double> >( hoNDArray< std::complex<double> >*, double, hoNDArray< std::complex<double> >* );
    template EXPORTCPUCOREMATH void pshrink< std::complex<double> >( hoNDArray< std::complex<double> >*, double,double, hoNDArray< std::complex<double> >* );
    template EXPORTCPUCOREMATH void shrinkd< std::complex<double> > ( hoNDArray< std::complex<double> >*, hoNDArray<double>*, double, hoNDArray< std::complex<double> >* );
    template EXPORTCPUCOREMATH void pshrinkd< std::complex<double> > ( hoNDArray< std::complex<double> >*, hoNDArray<double>*, double,double, hoNDArray< std::complex<double> >* );

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
    template EXPORTCPUCOREMATH void fill< complext<float> >( hoNDArray< complext<float> >*, complext<float> );
    template EXPORTCPUCOREMATH void clamp< complext<float> >( hoNDArray< complext<float> >*, float, float );
    template EXPORTCPUCOREMATH void clamp_min< complext<float> >( hoNDArray< complext<float> >*, float );
    template EXPORTCPUCOREMATH void clamp_max<complext<float> >( hoNDArray< complext<float> >*, float );
    template EXPORTCPUCOREMATH void normalize< complext<float> >( hoNDArray< complext<float> >*, float );
    template EXPORTCPUCOREMATH void shrink1< complext<float> >( hoNDArray< complext<float> >*, float, hoNDArray< complext<float> >* );
    template EXPORTCPUCOREMATH void pshrink< complext<float> >( hoNDArray< complext<float> >*, float,float, hoNDArray< complext<float> >* );
    template EXPORTCPUCOREMATH void shrinkd< complext<float> > ( hoNDArray< complext<float> >*, hoNDArray<float>*, float, hoNDArray< complext<float> >* );

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
    template EXPORTCPUCOREMATH void fill< complext<double> >( hoNDArray< complext<double> >*, complext<double> );
    template EXPORTCPUCOREMATH void clamp< complext<double> >( hoNDArray< complext<double> >*, double, double );
    template EXPORTCPUCOREMATH void clamp_min< complext<double> >( hoNDArray< complext<double> >*, double );
    template EXPORTCPUCOREMATH void clamp_max<complext<double> >( hoNDArray< complext<double> >*, double );
    template EXPORTCPUCOREMATH void normalize< complext<double> >( hoNDArray< complext<double> >*, double );
    template EXPORTCPUCOREMATH void shrink1< complext<double> >( hoNDArray< complext<double> >*, double, hoNDArray< complext<double> >* );
    template EXPORTCPUCOREMATH void pshrink< complext<double> >( hoNDArray< complext<double> >*, double,double, hoNDArray< complext<double> >* );
    template EXPORTCPUCOREMATH void shrinkd< complext<double> > ( hoNDArray< complext<double> >*, hoNDArray<double>*, double, hoNDArray< complext<double> >* );

    template EXPORTCPUCOREMATH boost::shared_ptr< hoNDArray< std::complex<float> > > real_to_complex< std::complex<float> >( hoNDArray<float>* );
    template EXPORTCPUCOREMATH boost::shared_ptr< hoNDArray<float_complext> > real_to_complex<float_complext>( hoNDArray<float>* );
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
    template EXPORTCPUCOREMATH boost::shared_ptr< hoNDArray<double_complext> > real_to_complex<double_complext>( hoNDArray<double>* );
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
