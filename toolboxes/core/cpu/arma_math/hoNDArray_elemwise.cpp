#include "hoNDArray_elemwise.h"
#include "hoNDArray_operators.h"
#include "hoNDArray_blas.h"
#include "complext.h"

#include <complex>

using namespace Gadgetron;

template<class T> boost::shared_ptr< hoNDArray<typename realType<T>::Type> > 
Gadgetron::abs( hoNDArray<T> *x )
{ 
  if( x == 0x0 )
    BOOST_THROW_EXCEPTION(runtime_error("Gadgetron::abs(): Invalid input array"));
   
  boost::shared_ptr< hoNDArray<typename realType<T>::Type> > result(new hoNDArray<typename realType<T>::Type>());
  result->create(x->get_dimensions());
  arma::Col<typename realType<T>::Type> aRes = as_arma_col(result.get());
  aRes = arma::abs(as_arma_col(x));
  return result;
}

template<class T> void 
Gadgetron::abs_inplace( hoNDArray<T> *x )
{ 
  if( x == 0x0 )
    BOOST_THROW_EXCEPTION(runtime_error("Gadgetron::abs_inplace(): Invalid input array"));
   
  arma::Col<typename realType<T>::Type> aRes = as_arma_col(x);
  aRes = arma::abs(aRes);
}  

template<class T> boost::shared_ptr< hoNDArray<typename realType<T>::Type> > 
Gadgetron::abs_square( hoNDArray<T> *x )
{ 
  if( x == 0x0 )
    BOOST_THROW_EXCEPTION(runtime_error("Gadgetron::abs_square(): Invalid input array"));
   
  boost::shared_ptr< hoNDArray<typename realType<T>::Type> > result(new hoNDArray<typename realType<T>::Type>());
  result->create(x->get_dimensions());
  arma::Col<typename realType<T>::Type> aRes = as_arma_col(result.get());
  aRes = arma::square(abs(as_arma_col(x)));
  return result;
}

template<class T> boost::shared_ptr< hoNDArray<T> > 
Gadgetron::sqrt( hoNDArray<T> *x )
{ 
  if( x == 0x0 )
    BOOST_THROW_EXCEPTION(runtime_error("Gadgetron::sqrt(): Invalid input array"));
   
  boost::shared_ptr< hoNDArray<T> > result(new hoNDArray<T>());
  result->create(x->get_dimensions());
  arma::Col<typename stdType<T>::Type> aRes = as_arma_col(result.get());
  aRes = arma::sqrt(as_arma_col(x));
  return result;
}

template<class T> void 
Gadgetron::sqrt_inplace( hoNDArray<T> *x )
{ 
  if( x == 0x0 )
    BOOST_THROW_EXCEPTION(runtime_error("Gadgetron::sqrt_inplace(): Invalid input array"));
   
  arma::Col<typename stdType<T>::Type> aRes = as_arma_col(x);
  aRes = arma::sqrt(aRes);
}
 
template<class T> boost::shared_ptr< hoNDArray<T> > Gadgetron::square( hoNDArray<T> *x )
{ 
  if( x == 0x0 )
    BOOST_THROW_EXCEPTION(runtime_error("Gadgetron::square(): Invalid input array"));
   
  boost::shared_ptr< hoNDArray<T> > result(new hoNDArray<T>());
  result->create(x->get_dimensions());
  arma::Col<typename stdType<T>::Type> aRes = as_arma_col(result.get());
  aRes = arma::square(as_arma_col(x));
  return result;
}

template<class T> void 
Gadgetron::square_inplace( hoNDArray<T> *x )
{ 
  if( x == 0x0 )
    BOOST_THROW_EXCEPTION(runtime_error("Gadgetron::square_inplace(): Invalid input array"));
   
  arma::Col<typename stdType<T>::Type> aRes = as_arma_col(x);
  aRes = arma::square(aRes);
}  
  
template<class T> boost::shared_ptr< hoNDArray<T> > 
Gadgetron::reciprocal( hoNDArray<T> *x )
{ 
  if( x == 0x0 )
    BOOST_THROW_EXCEPTION(runtime_error("Gadgetron::reciprocal(): Invalid input array"));
  
  arma::Col<typename stdType<T>::Type> ones(x->get_number_of_elements());
  ones.ones();
  boost::shared_ptr< hoNDArray<T> > result(new hoNDArray<T>());
  result->create(x->get_dimensions());
  arma::Col<typename stdType<T>::Type> aRes = as_arma_col(result.get());
  aRes = ones/as_arma_col(x);
  return result;
}

template<class T> void
Gadgetron::reciprocal_inplace( hoNDArray<T> *x )
{ 
  if( x == 0x0 )
    BOOST_THROW_EXCEPTION(runtime_error("Gadgetron::reciprocal_inplace(): Invalid input array"));
  
  arma::Col<typename stdType<T>::Type> aRes = as_arma_col(x);
  arma::Col<typename stdType<T>::Type> ones(x->get_number_of_elements());
  ones.ones();
  aRes = ones/aRes;
}
 
template<class T> boost::shared_ptr< hoNDArray<T> > 
Gadgetron::reciprocal_sqrt( hoNDArray<T> *x )
{ 
  if( x == 0x0 )
    BOOST_THROW_EXCEPTION(runtime_error("Gadgetron::reciprocal_sqrt(): Invalid input array"));
  
  arma::Col<typename stdType<T>::Type> ones(x->get_number_of_elements());
  ones.ones();   
  boost::shared_ptr< hoNDArray<T> > result(new hoNDArray<T>());
  result->create(x->get_dimensions());
  arma::Col<typename stdType<T>::Type> aRes = as_arma_col(result.get());
  aRes = ones/arma::sqrt(as_arma_col(x));
  return result;
}
 
template<class T> void Gadgetron::reciprocal_sqrt_inplace( hoNDArray<T> *x )
{ 
  if( x == 0x0 )
    BOOST_THROW_EXCEPTION(runtime_error("Gadgetron::reciprocal_sqrt_inplace(): Invalid input array"));
  
  arma::Col<typename stdType<T>::Type> ones(x->get_number_of_elements());
  ones.ones();
  arma::Col<typename stdType<T>::Type> aRes = as_arma_col(x);
  aRes = ones/arma::sqrt(aRes);
}
 
template<class T> boost::shared_ptr< hoNDArray<T> > 
Gadgetron::sgn( hoNDArray<T> *x )
{
  if( x == 0x0 )
    BOOST_THROW_EXCEPTION(runtime_error("Gadgetron::sgn(): Invalid input array"));
  
  boost::shared_ptr< hoNDArray<T> > res( new hoNDArray<T>() );
  res->create(x->get_dimensions());   
  for( int i = 0; i < res->get_number_of_elements(); i++ ){
    res->get_data_ptr()[i] = sgn(x->get_data_ptr()[i]);
  }
  return res;
}
 
template<class T> void Gadgetron::sgn_inplace( hoNDArray<T> *x )
{
  if( x == 0x0 )
    BOOST_THROW_EXCEPTION(runtime_error("Gadgetron::sgn_inplace(): Invalid input array"));
   
  for( int i = 0; i < x->get_number_of_elements(); i++ ) 
    x->get_data_ptr()[i] = sgn(x->get_data_ptr()[i]);
}

template<class T> boost::shared_ptr< hoNDArray<T> > 
Gadgetron::real( hoNDArray< std::complex<T> > *x )
{
  if( x == 0x0 )
    BOOST_THROW_EXCEPTION(runtime_error("Gadgetron::real(): Invalid input array"));
  
  boost::shared_ptr< hoNDArray<T> > result(new hoNDArray<T>());
  result->create(x->get_dimensions());
  arma::Col<T> aRes = as_arma_col(result.get());
  aRes = arma::real(as_arma_col(x));
  return result;
}

template<class T> boost::shared_ptr< hoNDArray<T> > 
Gadgetron::real( hoNDArray< complext<T> > *x )
{
  if( x == 0x0 )
    BOOST_THROW_EXCEPTION(runtime_error("Gadgetron::real(): Invalid input array"));
  
  boost::shared_ptr< hoNDArray<T> > result(new hoNDArray<T>());
  result->create(x->get_dimensions());
  arma::Col<T> aRes = as_arma_col(result.get());
  aRes = arma::real(as_arma_col(x));
  return result;
}

template<class T> boost::shared_ptr< hoNDArray<T> > 
Gadgetron::imag( hoNDArray< std::complex<T> > *x )
{
  if( x == 0x0 )
    BOOST_THROW_EXCEPTION(runtime_error("Gadgetron::imag(): Invalid input array"));
  
  boost::shared_ptr< hoNDArray<T> > result(new hoNDArray<T>());
  result->create(x->get_dimensions());
  arma::Col<T> aRes = as_arma_col(result.get());
  aRes = arma::imag(as_arma_col(x));
  return result;
}

template<class T> boost::shared_ptr< hoNDArray<T> > 
Gadgetron::imag( hoNDArray< complext<T> > *x )
{
  if( x == 0x0 )
    BOOST_THROW_EXCEPTION(runtime_error("Gadgetron::imag(): Invalid input array"));
  
  boost::shared_ptr< hoNDArray<T> > result(new hoNDArray<T>());
  result->create(x->get_dimensions());
  arma::Col<T> aRes = as_arma_col(result.get());
  aRes = arma::imag(as_arma_col(x));
  return result;
}

template<class T> boost::shared_ptr< hoNDArray< std::complex<T> > > 
Gadgetron::real_to_std_complex( hoNDArray<T> *x )
{
  if( x == 0x0 )
    BOOST_THROW_EXCEPTION(runtime_error("Gadgetron::real_to_std_complex(): Invalid input array"));
  
  boost::shared_ptr< hoNDArray< std::complex<T> > > result(new hoNDArray< std::complex<T> >());
  result->create(x->get_dimensions());
  arma::Col< std::complex<T> > aRes = as_arma_col(result.get());
  aRes = arma::Col< std::complex<T> >(as_arma_col(x), arma::Col<T>(x->get_number_of_elements()).zeros());
  return result;
}

template<class T> boost::shared_ptr< hoNDArray< complext<T> > > 
Gadgetron::real_to_complext( hoNDArray<T> *x )
{
  if( x == 0x0 )
    BOOST_THROW_EXCEPTION(runtime_error("Gadgetron::real_to_std_complex(): Invalid input array"));
  
  boost::shared_ptr< hoNDArray< complext<T> > > result(new hoNDArray< complext<T> >());
  result->create(x->get_dimensions());
  arma::Col< std::complex<T> > aRes = as_arma_col(result.get());
  aRes = arma::Col< std::complex<T> >(as_arma_col(x), arma::Col<T>(x->get_number_of_elements()).zeros());
  return result;
}

template<class T> void Gadgetron::clear( hoNDArray<T> *x )
{
  if( x == 0x0 )
    BOOST_THROW_EXCEPTION(runtime_error("Gadgetron::clear(): Invalid input array"));
   
  memset( x->get_data_ptr(), 0, x->get_number_of_elements()*sizeof(T));
}

template<class T> void Gadgetron::fill( hoNDArray<T> *x, T val )
{
  if( x == 0x0 )
    BOOST_THROW_EXCEPTION(runtime_error("Gadgetron::fill(): Invalid input array"));
  
  arma::Col<typename stdType<T>::Type> aRes = as_arma_col(x);
  aRes.fill(*((typename stdType<T>::Type*)&val));
}

template<class T> void 
Gadgetron::clamp( hoNDArray<T> *x, T min, T max )
{
  if( x == 0x0 )
    BOOST_THROW_EXCEPTION(runtime_error("Gadgetron::clamp_inplace(): Invalid input array"));
   
  for( int i = 0; i < x->get_number_of_elements(); i++ ){
    T tmp = x->get_data_ptr()[i];
    if( tmp < min ) 
      x->get_data_ptr()[i] = min;
    else if ( tmp > max) 
      x->get_data_ptr()[i] = max;
  }
}

template<class T> void 
Gadgetron::clamp_min( hoNDArray<T> *x, T min )
{
  if( x == 0x0 )
    BOOST_THROW_EXCEPTION(runtime_error("Gadgetron::clamp_min_inplace(): Invalid input array"));
   
  boost::shared_ptr< hoNDArray<T> > res( new hoNDArray<T>() );
  res->create( x->get_dimensions() );
   
  for( int i = 0; i < res->get_number_of_elements(); i++ ){
    T tmp = x->get_data_ptr()[i];
    if( tmp < min) 
      x->get_data_ptr()[i] = min;
  }
}

template<class T> void 
Gadgetron::clamp_max( hoNDArray<T> *x, T max )
{
  if( x == 0x0 )
    BOOST_THROW_EXCEPTION(runtime_error("Gadgetron::clamp_max_inplace(): Invalid input array"));
   
  boost::shared_ptr< hoNDArray<T> > res( new hoNDArray<T>() );
  res->create( x->get_dimensions() );
   
  for( int i = 0; i < res->get_number_of_elements(); i++ ){
    T tmp = x->get_data_ptr()[i];
    if( tmp > max) 
      x->get_data_ptr()[i] = max;
  }
}

template<class T> void 
Gadgetron::normalize( hoNDArray<T> *x, typename realType<T>::Type val )
{
  if( x == 0x0 )
    BOOST_THROW_EXCEPTION(runtime_error("Gadgetron::normalize(): Invalid input array"));
  
  unsigned int max_idx = amax(x);
  T max_val_before = x->get_data_ptr()[max_idx];
  typename realType<T>::Type scale = val/abs(max_val_before);
  *x *= scale;
}

template<class T> void 
Gadgetron::shrink1( hoNDArray<T> *x, typename realType<T>::Type gamma )
{
  if( x == 0x0 )
    BOOST_THROW_EXCEPTION(runtime_error("Gadgetron::shrink1(): Invalid input array"));
 
  for( int i = 0; i < x->get_number_of_elements(); i++ ) {
    T prev = x->get_data_ptr()[i];
    typename realType<T>::Type absPrev = abs(prev);
    T sgnPrev = (absPrev <= typename realType<T>::Type(0)) ? T(0) : prev/absPrev;
    x->get_data_ptr()[i] = sgnPrev*std::max(absPrev-gamma, typename realType<T>::Type(0));
  } 
}

template<class T> void 
Gadgetron::shrinkd ( hoNDArray<T> *_x, hoNDArray<typename realType<T>::Type> *_s, typename realType<T>::Type gamma )
{
  if( _x == 0x0  || _s == 0 )
    BOOST_THROW_EXCEPTION(runtime_error("Gadgetron::shrinkd(): Invalid input array"));
  
  for( int i = 0; i < _x->get_number_of_elements(); i++ ) {
    T x = _x->get_data_ptr()[i];
    typename realType<T>::Type s = _s->get_data_ptr()[i];
    _x->get_data_ptr()[i] = x/s*std::max(s-gamma,typename realType<T>::Type(0));
  } 
}


//
// Instantiation
//

template EXPORTCPUCOREMATH boost::shared_ptr< hoNDArray<float> > Gadgetron::abs<float>( hoNDArray<float>* );
template EXPORTCPUCOREMATH void Gadgetron::abs_inplace<float>( hoNDArray<float>* );
template EXPORTCPUCOREMATH boost::shared_ptr< hoNDArray<float> > Gadgetron::abs_square<float>( hoNDArray<float>* );
template EXPORTCPUCOREMATH boost::shared_ptr< hoNDArray<float> > Gadgetron::sqrt<float>( hoNDArray<float>* );
template EXPORTCPUCOREMATH void Gadgetron::sqrt_inplace<float>( hoNDArray<float>* );
template EXPORTCPUCOREMATH boost::shared_ptr< hoNDArray<float> > Gadgetron::square<float>( hoNDArray<float>* );
template EXPORTCPUCOREMATH void Gadgetron::square_inplace<float>( hoNDArray<float>* );
template EXPORTCPUCOREMATH boost::shared_ptr< hoNDArray<float> > Gadgetron::reciprocal<float>( hoNDArray<float>* );
template EXPORTCPUCOREMATH void Gadgetron::reciprocal_inplace<float>( hoNDArray<float>* );
template EXPORTCPUCOREMATH boost::shared_ptr< hoNDArray<float> > Gadgetron::reciprocal_sqrt<float>( hoNDArray<float>* );
template EXPORTCPUCOREMATH void Gadgetron::reciprocal_sqrt_inplace<float>( hoNDArray<float>* );
template EXPORTCPUCOREMATH boost::shared_ptr< hoNDArray<float> > Gadgetron::sgn<float>( hoNDArray<float>* );
template EXPORTCPUCOREMATH void Gadgetron::sgn_inplace<float>( hoNDArray<float>* );
template EXPORTCPUCOREMATH void Gadgetron::clear<float>( hoNDArray<float>* );
template EXPORTCPUCOREMATH void Gadgetron::fill<float>( hoNDArray<float>*, float );
template EXPORTCPUCOREMATH void Gadgetron::clamp<float>( hoNDArray<float>*, float, float );
template EXPORTCPUCOREMATH void Gadgetron::clamp_min<float>( hoNDArray<float>*, float );
template EXPORTCPUCOREMATH void Gadgetron::clamp_max<float>( hoNDArray<float>*, float );
template EXPORTCPUCOREMATH void Gadgetron::normalize<float>( hoNDArray<float>*, float );
template EXPORTCPUCOREMATH void Gadgetron::shrink1<float>( hoNDArray<float>*, float );
template EXPORTCPUCOREMATH void Gadgetron::shrinkd<float> ( hoNDArray<float>*, hoNDArray<float>*, float );

template EXPORTCPUCOREMATH boost::shared_ptr< hoNDArray<double> > Gadgetron::abs<double>( hoNDArray<double>* );
template EXPORTCPUCOREMATH void Gadgetron::abs_inplace<double>( hoNDArray<double>* );
template EXPORTCPUCOREMATH boost::shared_ptr< hoNDArray<double> > Gadgetron::abs_square<double>( hoNDArray<double>* );
template EXPORTCPUCOREMATH boost::shared_ptr< hoNDArray<double> > Gadgetron::sqrt<double>( hoNDArray<double>* );
template EXPORTCPUCOREMATH void Gadgetron::sqrt_inplace<double>( hoNDArray<double>* );
template EXPORTCPUCOREMATH boost::shared_ptr< hoNDArray<double> > Gadgetron::square<double>( hoNDArray<double>* );
template EXPORTCPUCOREMATH void Gadgetron::square_inplace<double>( hoNDArray<double>* );
template EXPORTCPUCOREMATH boost::shared_ptr< hoNDArray<double> > Gadgetron::reciprocal<double>( hoNDArray<double>* );
template EXPORTCPUCOREMATH void Gadgetron::reciprocal_inplace<double>( hoNDArray<double>* );
template EXPORTCPUCOREMATH boost::shared_ptr< hoNDArray<double> > Gadgetron::reciprocal_sqrt<double>( hoNDArray<double>* );
template EXPORTCPUCOREMATH void Gadgetron::reciprocal_sqrt_inplace<double>( hoNDArray<double>* );
template EXPORTCPUCOREMATH boost::shared_ptr< hoNDArray<double> > Gadgetron::sgn<double>( hoNDArray<double>* );
template EXPORTCPUCOREMATH void Gadgetron::sgn_inplace<double>( hoNDArray<double>* );
template EXPORTCPUCOREMATH void Gadgetron::clear<double>( hoNDArray<double>* );
template EXPORTCPUCOREMATH void Gadgetron::fill<double>( hoNDArray<double>*, double );
template EXPORTCPUCOREMATH void Gadgetron::clamp<double>( hoNDArray<double>*, double, double );
template EXPORTCPUCOREMATH void Gadgetron::clamp_min<double>( hoNDArray<double>*, double );
template EXPORTCPUCOREMATH void Gadgetron::clamp_max<double>( hoNDArray<double>*, double );
template EXPORTCPUCOREMATH void Gadgetron::normalize<double>( hoNDArray<double>*, double );
template EXPORTCPUCOREMATH void Gadgetron::shrink1<double>( hoNDArray<double>*, double );
template EXPORTCPUCOREMATH void Gadgetron::shrinkd<double> ( hoNDArray<double>*, hoNDArray<double>*, double );

template EXPORTCPUCOREMATH boost::shared_ptr< hoNDArray<float> > Gadgetron::abs< std::complex<float> >( hoNDArray< std::complex<float> >* );
template EXPORTCPUCOREMATH boost::shared_ptr< hoNDArray<float> > Gadgetron::abs_square< std::complex<float> >( hoNDArray< std::complex<float> >* );
template EXPORTCPUCOREMATH boost::shared_ptr< hoNDArray< std::complex<float> > > Gadgetron::sqrt< std::complex<float> >( hoNDArray< std::complex<float> >* );
template EXPORTCPUCOREMATH void Gadgetron::sqrt_inplace< std::complex<float> >( hoNDArray< std::complex<float> >* );
template EXPORTCPUCOREMATH boost::shared_ptr< hoNDArray< std::complex<float> > > Gadgetron::square< std::complex<float> >( hoNDArray< std::complex<float> >* );
template EXPORTCPUCOREMATH void Gadgetron::square_inplace< std::complex<float> >( hoNDArray< std::complex<float> >* );
template EXPORTCPUCOREMATH boost::shared_ptr< hoNDArray< std::complex<float> > > Gadgetron::reciprocal< std::complex<float> >( hoNDArray< std::complex<float> >* );
template EXPORTCPUCOREMATH void Gadgetron::reciprocal_inplace< std::complex<float> >( hoNDArray< std::complex<float> >* );
template EXPORTCPUCOREMATH boost::shared_ptr< hoNDArray< std::complex<float> > > Gadgetron::reciprocal_sqrt< std::complex<float> >( hoNDArray< std::complex<float> >* );
template EXPORTCPUCOREMATH void Gadgetron::reciprocal_sqrt_inplace< std::complex<float> >( hoNDArray< std::complex<float> >* );
template EXPORTCPUCOREMATH void Gadgetron::clear< std::complex<float> >( hoNDArray< std::complex<float> >* );
template EXPORTCPUCOREMATH void Gadgetron::fill< std::complex<float> >( hoNDArray< std::complex<float> >*, std::complex<float> );
template EXPORTCPUCOREMATH void Gadgetron::normalize< std::complex<float> >( hoNDArray< std::complex<float> >*, float );
template EXPORTCPUCOREMATH void Gadgetron::shrink1< std::complex<float> >( hoNDArray< std::complex<float> >*, float );
template EXPORTCPUCOREMATH void Gadgetron::shrinkd< std::complex<float> > ( hoNDArray< std::complex<float> >*, hoNDArray<float>*, float );

template EXPORTCPUCOREMATH boost::shared_ptr< hoNDArray<double> > Gadgetron::abs< std::complex<double> >( hoNDArray< std::complex<double> >* );
template EXPORTCPUCOREMATH boost::shared_ptr< hoNDArray<double> > Gadgetron::abs_square< std::complex<double> >( hoNDArray< std::complex<double> >* );
template EXPORTCPUCOREMATH boost::shared_ptr< hoNDArray< std::complex<double> > > Gadgetron::sqrt< std::complex<double> >( hoNDArray< std::complex<double> >* );
template EXPORTCPUCOREMATH void Gadgetron::sqrt_inplace< std::complex<double> >( hoNDArray< std::complex<double> >* );
template EXPORTCPUCOREMATH boost::shared_ptr< hoNDArray< std::complex<double> > > Gadgetron::square< std::complex<double> >( hoNDArray< std::complex<double> >* );
template EXPORTCPUCOREMATH void Gadgetron::square_inplace< std::complex<double> >( hoNDArray< std::complex<double> >* );
template EXPORTCPUCOREMATH boost::shared_ptr< hoNDArray< std::complex<double> > > Gadgetron::reciprocal< std::complex<double> >( hoNDArray< std::complex<double> >* );
template EXPORTCPUCOREMATH void Gadgetron::reciprocal_inplace< std::complex<double> >( hoNDArray< std::complex<double> >* );
template EXPORTCPUCOREMATH boost::shared_ptr< hoNDArray< std::complex<double> > > Gadgetron::reciprocal_sqrt< std::complex<double> >( hoNDArray< std::complex<double> >* );
template EXPORTCPUCOREMATH void Gadgetron::reciprocal_sqrt_inplace< std::complex<double> >( hoNDArray< std::complex<double> >* );
template EXPORTCPUCOREMATH void Gadgetron::clear< std::complex<double> >( hoNDArray< std::complex<double> >* );
template EXPORTCPUCOREMATH void Gadgetron::fill< std::complex<double> >( hoNDArray< std::complex<double> >*, std::complex<double> );
template EXPORTCPUCOREMATH void Gadgetron::normalize< std::complex<double> >( hoNDArray< std::complex<double> >*, double );
template EXPORTCPUCOREMATH void Gadgetron::shrink1< std::complex<double> >( hoNDArray< std::complex<double> >*, double );
template EXPORTCPUCOREMATH void Gadgetron::shrinkd< std::complex<double> > ( hoNDArray< std::complex<double> >*, hoNDArray<double>*, double );

template EXPORTCPUCOREMATH boost::shared_ptr< hoNDArray<float> > Gadgetron::abs< complext<float> >( hoNDArray< complext<float> >* );
template EXPORTCPUCOREMATH boost::shared_ptr< hoNDArray<float> > Gadgetron::abs_square< complext<float> >( hoNDArray< complext<float> >* );
template EXPORTCPUCOREMATH boost::shared_ptr< hoNDArray< complext<float> > > Gadgetron::sqrt< complext<float> >( hoNDArray< complext<float> >* );
template EXPORTCPUCOREMATH void Gadgetron::sqrt_inplace< complext<float> >( hoNDArray< complext<float> >* );
template EXPORTCPUCOREMATH boost::shared_ptr< hoNDArray< complext<float> > > Gadgetron::square< complext<float> >( hoNDArray< complext<float> >* );
template EXPORTCPUCOREMATH void Gadgetron::square_inplace< complext<float> >( hoNDArray< complext<float> >* );
template EXPORTCPUCOREMATH boost::shared_ptr< hoNDArray< complext<float> > > Gadgetron::reciprocal< complext<float> >( hoNDArray< complext<float> >* );
template EXPORTCPUCOREMATH void Gadgetron::reciprocal_inplace< complext<float> >( hoNDArray< complext<float> >* );
template EXPORTCPUCOREMATH boost::shared_ptr< hoNDArray< complext<float> > > Gadgetron::reciprocal_sqrt< complext<float> >( hoNDArray< complext<float> >* );
template EXPORTCPUCOREMATH void Gadgetron::reciprocal_sqrt_inplace< complext<float> >( hoNDArray< complext<float> >* );
template EXPORTCPUCOREMATH void Gadgetron::clear< complext<float> >( hoNDArray< complext<float> >* );
template EXPORTCPUCOREMATH void Gadgetron::fill< complext<float> >( hoNDArray< complext<float> >*, complext<float> );
template EXPORTCPUCOREMATH void Gadgetron::normalize< complext<float> >( hoNDArray< complext<float> >*, float );
template EXPORTCPUCOREMATH void Gadgetron::shrink1< complext<float> >( hoNDArray< complext<float> >*, float );
template EXPORTCPUCOREMATH void Gadgetron::shrinkd< complext<float> > ( hoNDArray< complext<float> >*, hoNDArray<float>*, float );

template EXPORTCPUCOREMATH boost::shared_ptr< hoNDArray<double> > Gadgetron::abs< complext<double> >( hoNDArray< complext<double> >* );
template EXPORTCPUCOREMATH boost::shared_ptr< hoNDArray<double> > Gadgetron::abs_square< complext<double> >( hoNDArray< complext<double> >* );
template EXPORTCPUCOREMATH boost::shared_ptr< hoNDArray< complext<double> > > Gadgetron::sqrt< complext<double> >( hoNDArray< complext<double> >* );
template EXPORTCPUCOREMATH void Gadgetron::sqrt_inplace< complext<double> >( hoNDArray< complext<double> >* );
template EXPORTCPUCOREMATH boost::shared_ptr< hoNDArray< complext<double> > > Gadgetron::square< complext<double> >( hoNDArray< complext<double> >* );
template EXPORTCPUCOREMATH void Gadgetron::square_inplace< complext<double> >( hoNDArray< complext<double> >* );
template EXPORTCPUCOREMATH boost::shared_ptr< hoNDArray< complext<double> > > Gadgetron::reciprocal< complext<double> >( hoNDArray< complext<double> >* );
template EXPORTCPUCOREMATH void Gadgetron::reciprocal_inplace< complext<double> >( hoNDArray< complext<double> >* );
template EXPORTCPUCOREMATH boost::shared_ptr< hoNDArray< complext<double> > > Gadgetron::reciprocal_sqrt< complext<double> >( hoNDArray< complext<double> >* );
template EXPORTCPUCOREMATH void Gadgetron::reciprocal_sqrt_inplace< complext<double> >( hoNDArray< complext<double> >* );
template EXPORTCPUCOREMATH void Gadgetron::clear< complext<double> >( hoNDArray< complext<double> >* );
template EXPORTCPUCOREMATH void Gadgetron::fill< complext<double> >( hoNDArray< complext<double> >*, complext<double> );
template EXPORTCPUCOREMATH void Gadgetron::normalize< complext<double> >( hoNDArray< complext<double> >*, double );
template EXPORTCPUCOREMATH void Gadgetron::shrink1< complext<double> >( hoNDArray< complext<double> >*, double );
template EXPORTCPUCOREMATH void Gadgetron::shrinkd< complext<double> > ( hoNDArray< complext<double> >*, hoNDArray<double>*, double );

template EXPORTCPUCOREMATH boost::shared_ptr< hoNDArray< std::complex<float> > > Gadgetron::real_to_std_complex<float>( hoNDArray<float>* );
template EXPORTCPUCOREMATH boost::shared_ptr< hoNDArray< complext<float> > > Gadgetron::real_to_complext<float>( hoNDArray<float>* );
template EXPORTCPUCOREMATH boost::shared_ptr< hoNDArray<float> > Gadgetron::real<float>( hoNDArray< std::complex<float> >* );
template EXPORTCPUCOREMATH boost::shared_ptr< hoNDArray<float> > Gadgetron::real<float>( hoNDArray< complext<float> >* );
template EXPORTCPUCOREMATH boost::shared_ptr< hoNDArray<float> > Gadgetron::imag<float>( hoNDArray< std::complex<float> >* );
template EXPORTCPUCOREMATH boost::shared_ptr< hoNDArray<float> > Gadgetron::imag<float>( hoNDArray< complext<float> >* );

template EXPORTCPUCOREMATH boost::shared_ptr< hoNDArray< std::complex<double> > > Gadgetron::real_to_std_complex<double>( hoNDArray<double>* );
template EXPORTCPUCOREMATH boost::shared_ptr< hoNDArray< complext<double> > > Gadgetron::real_to_complext<double>( hoNDArray<double>* );
template EXPORTCPUCOREMATH boost::shared_ptr< hoNDArray<double> > Gadgetron::real<double>( hoNDArray< std::complex<double> >* );
template EXPORTCPUCOREMATH boost::shared_ptr< hoNDArray<double> > Gadgetron::real<double>( hoNDArray< complext<double> >* );
template EXPORTCPUCOREMATH boost::shared_ptr< hoNDArray<double> > Gadgetron::imag<double>( hoNDArray< std::complex<double> >* );
template EXPORTCPUCOREMATH boost::shared_ptr< hoNDArray<double> > Gadgetron::imag<double>( hoNDArray< complext<double> >* );


