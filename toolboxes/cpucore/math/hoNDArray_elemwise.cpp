#include "hoNDArray_elemwise.h"
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

template<class T> boost::shared_ptr< hoNDArray<T> > 
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
Gadgetron::clamp( hoNDArray<T> *x, T min, T max )
{
  if( x == 0x0 )
    BOOST_THROW_EXCEPTION(runtime_error("Gadgetron::clamp(): Invalid input array"));
   
  boost::shared_ptr< hoNDArray<T> > res( new hoNDArray<T>() );
  res->create(x->get_dimensions());
   
  for( int i = 0; i < res->get_number_of_elements(); i++ ){
    T tmp = x->get_data_ptr()[i];
    if( tmp < min ) 
      res->get_data_ptr()[i] = min;
    else if ( tmp > max) 
      res->get_data_ptr()[i] = max;
    else
      res->get_data_ptr()[i] = tmp;
  }
  return res;
}

template<class T> void 
Gadgetron::clamp_inplace( hoNDArray<T> *x, T min, T max )
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

template<class T> boost::shared_ptr< hoNDArray<T> > 
Gadgetron::clamp_min( hoNDArray<T> *x, T min )
{
  if( x == 0x0 )
    BOOST_THROW_EXCEPTION(runtime_error("Gadgetron::clamp_min(): Invalid input array"));
    
  boost::shared_ptr< hoNDArray<T> > res( new hoNDArray<T>() );
  res->create( x->get_dimensions() );
    
  for( int i = 0; i < res->get_number_of_elements(); i++ ){
    T tmp = x->get_data_ptr()[i];
    if( tmp < min) 
      x->get_data_ptr()[i] = min;
    else
      x->get_data_ptr()[i] = tmp;
  }
  return res;
}

template<class T> void 
Gadgetron::clamp_min_inplace( hoNDArray<T> *x, T min )
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

template<class T> boost::shared_ptr< hoNDArray<T> > 
Gadgetron::clamp_max( hoNDArray<T> *x, T max )
{
  if( x == 0x0 )
    BOOST_THROW_EXCEPTION(runtime_error("Gadgetron::clamp_max(): Invalid input array"));
    
  boost::shared_ptr< hoNDArray<T> > res( new hoNDArray<T>() );
  res->create( x->get_dimensions() );
    
  for( int i = 0; i < res->get_number_of_elements(); i++ ){
    T tmp = x->get_data_ptr()[i];
    if( tmp > max) 
      x->get_data_ptr()[i] = max;
    else
      x->get_data_ptr()[i] = tmp;
  }
  return res;
}

template<class T> void 
Gadgetron::clamp_max_inplace( hoNDArray<T> *x, T max )
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

//
// Instantiation
//

template EXPORTCPUCORE boost::shared_ptr< hoNDArray<float> > Gadgetron::abs<float>( hoNDArray<float>* );
template EXPORTCPUCORE void Gadgetron::abs_inplace<float>( hoNDArray<float>* );
template EXPORTCPUCORE boost::shared_ptr< hoNDArray<float> > Gadgetron::sqrt<float>( hoNDArray<float>* );
template EXPORTCPUCORE void Gadgetron::sqrt_inplace<float>( hoNDArray<float>* );
template EXPORTCPUCORE boost::shared_ptr< hoNDArray<float> > Gadgetron::square<float>( hoNDArray<float>* );
template EXPORTCPUCORE void Gadgetron::square_inplace<float>( hoNDArray<float>* );
template EXPORTCPUCORE boost::shared_ptr< hoNDArray<float> > Gadgetron::reciprocal<float>( hoNDArray<float>* );
template EXPORTCPUCORE boost::shared_ptr< hoNDArray<float> > Gadgetron::reciprocal_inplace<float>( hoNDArray<float>* );
template EXPORTCPUCORE boost::shared_ptr< hoNDArray<float> > Gadgetron::reciprocal_sqrt<float>( hoNDArray<float>* );
template EXPORTCPUCORE void Gadgetron::reciprocal_sqrt_inplace<float>( hoNDArray<float>* );
template EXPORTCPUCORE boost::shared_ptr< hoNDArray<float> > Gadgetron::sgn<float>( hoNDArray<float>* );
template EXPORTCPUCORE void Gadgetron::sgn_inplace<float>( hoNDArray<float>* );
template EXPORTCPUCORE boost::shared_ptr< hoNDArray<float> > Gadgetron::clamp<float>( hoNDArray<float>*, float, float );
template EXPORTCPUCORE void Gadgetron::clamp_inplace<float>( hoNDArray<float>*, float, float );
template EXPORTCPUCORE boost::shared_ptr< hoNDArray<float> > Gadgetron::clamp_min<float>( hoNDArray<float>*, float );
template EXPORTCPUCORE void Gadgetron::clamp_min_inplace<float>( hoNDArray<float>*, float );
template EXPORTCPUCORE boost::shared_ptr< hoNDArray<float> > Gadgetron::clamp_max<float>( hoNDArray<float>*, float );
template EXPORTCPUCORE void Gadgetron::clamp_max_inplace<float>( hoNDArray<float>*, float );

template EXPORTCPUCORE boost::shared_ptr< hoNDArray<double> > Gadgetron::abs<double>( hoNDArray<double>* );
template EXPORTCPUCORE void Gadgetron::abs_inplace<double>( hoNDArray<double>* );
template EXPORTCPUCORE boost::shared_ptr< hoNDArray<double> > Gadgetron::sqrt<double>( hoNDArray<double>* );
template EXPORTCPUCORE void Gadgetron::sqrt_inplace<double>( hoNDArray<double>* );
template EXPORTCPUCORE boost::shared_ptr< hoNDArray<double> > Gadgetron::square<double>( hoNDArray<double>* );
template EXPORTCPUCORE void Gadgetron::square_inplace<double>( hoNDArray<double>* );
template EXPORTCPUCORE boost::shared_ptr< hoNDArray<double> > Gadgetron::reciprocal<double>( hoNDArray<double>* );
template EXPORTCPUCORE boost::shared_ptr< hoNDArray<double> > Gadgetron::reciprocal_inplace<double>( hoNDArray<double>* );
template EXPORTCPUCORE boost::shared_ptr< hoNDArray<double> > Gadgetron::reciprocal_sqrt<double>( hoNDArray<double>* );
template EXPORTCPUCORE void Gadgetron::reciprocal_sqrt_inplace<double>( hoNDArray<double>* );
template EXPORTCPUCORE boost::shared_ptr< hoNDArray<double> > Gadgetron::sgn<double>( hoNDArray<double>* );
template EXPORTCPUCORE void Gadgetron::sgn_inplace<double>( hoNDArray<double>* );
template EXPORTCPUCORE boost::shared_ptr< hoNDArray<double> > Gadgetron::clamp<double>( hoNDArray<double>*, double, double );
template EXPORTCPUCORE void Gadgetron::clamp_inplace<double>( hoNDArray<double>*, double, double );
template EXPORTCPUCORE boost::shared_ptr< hoNDArray<double> > Gadgetron::clamp_min<double>( hoNDArray<double>*, double );
template EXPORTCPUCORE void Gadgetron::clamp_min_inplace<double>( hoNDArray<double>*, double );
template EXPORTCPUCORE boost::shared_ptr< hoNDArray<double> > Gadgetron::clamp_max<double>( hoNDArray<double>*, double );
template EXPORTCPUCORE void Gadgetron::clamp_max_inplace<double>( hoNDArray<double>*, double );

template EXPORTCPUCORE boost::shared_ptr< hoNDArray<float> > Gadgetron::abs< std::complex<float> >( hoNDArray< std::complex<float> >* );
template EXPORTCPUCORE boost::shared_ptr< hoNDArray< std::complex<float> > > Gadgetron::sqrt< std::complex<float> >( hoNDArray< std::complex<float> >* );
template EXPORTCPUCORE void Gadgetron::sqrt_inplace< std::complex<float> >( hoNDArray< std::complex<float> >* );
template EXPORTCPUCORE boost::shared_ptr< hoNDArray< std::complex<float> > > Gadgetron::square< std::complex<float> >( hoNDArray< std::complex<float> >* );
template EXPORTCPUCORE void Gadgetron::square_inplace< std::complex<float> >( hoNDArray< std::complex<float> >* );
template EXPORTCPUCORE boost::shared_ptr< hoNDArray< std::complex<float> > > Gadgetron::reciprocal< std::complex<float> >( hoNDArray< std::complex<float> >* );
template EXPORTCPUCORE boost::shared_ptr< hoNDArray< std::complex<float> > > Gadgetron::reciprocal_inplace< std::complex<float> >( hoNDArray< std::complex<float> >* );
template EXPORTCPUCORE boost::shared_ptr< hoNDArray< std::complex<float> > > Gadgetron::reciprocal_sqrt< std::complex<float> >( hoNDArray< std::complex<float> >* );
template EXPORTCPUCORE void Gadgetron::reciprocal_sqrt_inplace< std::complex<float> >( hoNDArray< std::complex<float> >* );

template EXPORTCPUCORE boost::shared_ptr< hoNDArray<double> > Gadgetron::abs< std::complex<double> >( hoNDArray< std::complex<double> >* );
template EXPORTCPUCORE boost::shared_ptr< hoNDArray< std::complex<double> > > Gadgetron::sqrt< std::complex<double> >( hoNDArray< std::complex<double> >* );
template EXPORTCPUCORE void Gadgetron::sqrt_inplace< std::complex<double> >( hoNDArray< std::complex<double> >* );
template EXPORTCPUCORE boost::shared_ptr< hoNDArray< std::complex<double> > > Gadgetron::square< std::complex<double> >( hoNDArray< std::complex<double> >* );
template EXPORTCPUCORE void Gadgetron::square_inplace< std::complex<double> >( hoNDArray< std::complex<double> >* );
template EXPORTCPUCORE boost::shared_ptr< hoNDArray< std::complex<double> > > Gadgetron::reciprocal< std::complex<double> >( hoNDArray< std::complex<double> >* );
template EXPORTCPUCORE boost::shared_ptr< hoNDArray< std::complex<double> > > Gadgetron::reciprocal_inplace< std::complex<double> >( hoNDArray< std::complex<double> >* );
template EXPORTCPUCORE boost::shared_ptr< hoNDArray< std::complex<double> > > Gadgetron::reciprocal_sqrt< std::complex<double> >( hoNDArray< std::complex<double> >* );
template EXPORTCPUCORE void Gadgetron::reciprocal_sqrt_inplace< std::complex<double> >( hoNDArray< std::complex<double> >* );

template EXPORTCPUCORE boost::shared_ptr< hoNDArray<float> > Gadgetron::abs< complext<float> >( hoNDArray< complext<float> >* );
template EXPORTCPUCORE boost::shared_ptr< hoNDArray< complext<float> > > Gadgetron::sqrt< complext<float> >( hoNDArray< complext<float> >* );
template EXPORTCPUCORE void Gadgetron::sqrt_inplace< complext<float> >( hoNDArray< complext<float> >* );
template EXPORTCPUCORE boost::shared_ptr< hoNDArray< complext<float> > > Gadgetron::square< complext<float> >( hoNDArray< complext<float> >* );
template EXPORTCPUCORE void Gadgetron::square_inplace< complext<float> >( hoNDArray< complext<float> >* );
template EXPORTCPUCORE boost::shared_ptr< hoNDArray< complext<float> > > Gadgetron::reciprocal< complext<float> >( hoNDArray< complext<float> >* );
template EXPORTCPUCORE boost::shared_ptr< hoNDArray< complext<float> > > Gadgetron::reciprocal_inplace< complext<float> >( hoNDArray< complext<float> >* );
template EXPORTCPUCORE boost::shared_ptr< hoNDArray< complext<float> > > Gadgetron::reciprocal_sqrt< complext<float> >( hoNDArray< complext<float> >* );
template EXPORTCPUCORE void Gadgetron::reciprocal_sqrt_inplace< complext<float> >( hoNDArray< complext<float> >* );

template EXPORTCPUCORE boost::shared_ptr< hoNDArray<double> > Gadgetron::abs< complext<double> >( hoNDArray< complext<double> >* );
template EXPORTCPUCORE boost::shared_ptr< hoNDArray< complext<double> > > Gadgetron::sqrt< complext<double> >( hoNDArray< complext<double> >* );
template EXPORTCPUCORE void Gadgetron::sqrt_inplace< complext<double> >( hoNDArray< complext<double> >* );
template EXPORTCPUCORE boost::shared_ptr< hoNDArray< complext<double> > > Gadgetron::square< complext<double> >( hoNDArray< complext<double> >* );
template EXPORTCPUCORE void Gadgetron::square_inplace< complext<double> >( hoNDArray< complext<double> >* );
template EXPORTCPUCORE boost::shared_ptr< hoNDArray< complext<double> > > Gadgetron::reciprocal< complext<double> >( hoNDArray< complext<double> >* );
template EXPORTCPUCORE boost::shared_ptr< hoNDArray< complext<double> > > Gadgetron::reciprocal_inplace< complext<double> >( hoNDArray< complext<double> >* );
template EXPORTCPUCORE boost::shared_ptr< hoNDArray< complext<double> > > Gadgetron::reciprocal_sqrt< complext<double> >( hoNDArray< complext<double> >* );
template EXPORTCPUCORE void Gadgetron::reciprocal_sqrt_inplace< complext<double> >( hoNDArray< complext<double> >* );
