#include "hoNDArray_operators.h"

using namespace Gadgetron;

// Private utility to verify array dimensions. 
// It "replaces" NDArray::dimensions_equal() to support batch mode.
//
template<class T,class S> static bool compatible_dimensions( hoNDArray<T> &x, hoNDArray<S> &y )
{
  bool retVal = true;
  for (int i = 0; i < y.get_number_of_dimensions(); i++){
    retVal &= (x.get_size(i) == y.get_size(i));
  }
  return retVal;
}


template<class T> hoNDArray<T>& Gadgetron::operator+= (hoNDArray<T> &x, const hoNDArray<T> &y)
{
  if( compatible_dimensions(x,y) ){
    unsigned int num_batches = x.get_number_of_elements()/y.get_number_of_elements();
    for( unsigned int batch=0; batch<num_batches; batch++ ){	
      hoNDArray tmp;
      tmp.create( y->get_dimensions(), x.get_data_ptr()+batch*y.get_number_of_elements() );
      arma::Col<typename stdType<T>::Type> aRes = as_arma_col(&tmp);
      aRes += as_arma_col(&y);
      return x;
    }
  } else {
    BOOST_THROW_EXCEPTION(runtime_error("hoNDArray::operator+=: Incompatible array dimensions"));
  }
}

template<class T> hoNDArray<T>& Gadgetron::operator+= (hoNDArray<T> &x, const T &y)
{
  arma::Col<typename stdType<T>::Type> aRes = as_arma_col(&x);
  aRes += y;
  return x;  
}

template<class T> hoNDArray<T>& Gadgetron::operator-= (hoNDArray<T> &x, const hoNDArray<T> &y)
{
  if( compatible_dimensions(x,y) ){
    unsigned int num_batches = x.get_number_of_elements()/y.get_number_of_elements();
    for( unsigned int batch=0; batch<num_batches; batch++ ){	
      hoNDArray tmp;
      tmp.create( y->get_dimensions(), x.get_data_ptr()+batch*y.get_number_of_elements() );
      arma::Col<typename stdType<T>::Type> aRes = as_arma_col(&tmp);
      aRes -= as_arma_col(&y);
      return x;
    }
  } else {
    BOOST_THROW_EXCEPTION(runtime_error("hoNDArray::operator-=: Incompatible array dimensions"));
  }
}

template<class T> hoNDArray<T>& Gadgetron::operator-= (hoNDArray<T> &x, const T &y)
{
  arma::Col<typename stdType<T>::Type> aRes = as_arma_col(&x);
  aRes -= y;
  return x;  
}

template<class T> hoNDArray<T>& Gadgetron::operator*= (hoNDArray<T> &x, const hoNDArray<T> &y)
{
  if( compatible_dimensions(x,y) ){
    unsigned int num_batches = x.get_number_of_elements()/y.get_number_of_elements();
    for( unsigned int batch=0; batch<num_batches; batch++ ){	
      hoNDArray tmp;
      tmp.create( y->get_dimensions(), x.get_data_ptr()+batch*y.get_number_of_elements() );
      arma::Col<typename stdType<T>::Type> aRes = as_arma_col(&tmp);
      aRes *= as_arma_col(&y);
      return x;
    }
  } else {
    BOOST_THROW_EXCEPTION(runtime_error("hoNDArray::operator*=: Incompatible array dimensions"));
  }
}

template<class T> hoNDArray<T>& Gadgetron::operator*= (hoNDArray<T> &x, const T &y)
{
  arma::Col<typename stdType<T>::Type> aRes = as_arma_col(&x);
  aRes *= y;
  return x;  
}

template<class T> hoNDArray<T>& Gadgetron::operator/= (hoNDArray<T> &x, const hoNDArray<T> &y)
{
  if( compatible_dimensions(x,y) ){
    unsigned int num_batches = x.get_number_of_elements()/y.get_number_of_elements();
    for( unsigned int batch=0; batch<num_batches; batch++ ){	
      hoNDArray tmp;
      tmp.create( y->get_dimensions(), x.get_data_ptr()+batch*y.get_number_of_elements() );
      arma::Col<typename stdType<T>::Type> aRes = as_arma_col(&tmp);
      aRes /= as_arma_col(&y);
      return x;
    }
  } else {
    BOOST_THROW_EXCEPTION(runtime_error("hoNDArray::operator/=: Incompatible array dimensions"));
  }
}

template<class T> hoNDArray<T>& Gadgetron::operator/= (hoNDArray<T> &x, const T &y)
{
  arma::Col<typename stdType<T>::Type> aRes = as_arma_col(&x);
  aRes /= y;
  return x;  
}
