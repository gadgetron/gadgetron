#include "cuNDArray_elemwise.h"
#include "cuNDArray_operators.h"
#include "cuNDArray_blas.h"
#include "complext.h"

#include <complex>
#include <thrust/functional.h>

using namespace Gadgetron;
using namespace std;

template<typename T> struct cuNDA_abs : public thrust::unary_function<T,typename realType<T>::Type>
{
  __device__ typename Gadgetron::realType<T>::Type operator()(const T &x) const {return abs(x);}
};

template<class T> boost::shared_ptr< cuNDArray<typename realType<T>::Type> > 
Gadgetron::abs( cuNDArray<T> *x )
{ 
  if( x == 0x0 )
    BOOST_THROW_EXCEPTION(runtime_error("Gadgetron::abs(): Invalid input array"));
   
  boost::shared_ptr< cuNDArray<typename realType<T>::Type> > result(new cuNDArray<typename realType<T>::Type>());
  result->create(x->get_dimensions());
  thrust::device_ptr<typename realType<T>::Type> resPtr = result->get_device_ptr();
  thrust::device_ptr<T> xPtr = x->get_device_ptr();
  thrust::transform(xPtr,xPtr+x->get_number_of_elements(),resPtr,cuNDA_abs<T>());
  return result;
}

template<class T> void 
Gadgetron::abs_inplace( cuNDArray<T> *x )
{ 
  if( x == 0x0 )
    BOOST_THROW_EXCEPTION(runtime_error("Gadgetron::abs_inplace(): Invalid input array"));
   
  thrust::device_ptr<T> xPtr = x->get_device_ptr();
  thrust::transform(xPtr,xPtr+x->get_number_of_elements(),xPtr,cuNDA_abs<T>());
}  
  
template<typename T> struct cuNDA_abs_square : public thrust::unary_function<T,typename realType<T>::Type>
{
  __device__ typename Gadgetron::realType<T>::Type operator()(const T &x) const 
  { 
    typename realType<T>::Type tmp = abs(x);
    return tmp*tmp;
  }
};

template<class T> boost::shared_ptr< cuNDArray<typename realType<T>::Type> > 
Gadgetron::abs_square( cuNDArray<T> *x )
{ 
  if( x == 0x0 )
    BOOST_THROW_EXCEPTION(runtime_error("Gadgetron::abs_square(): Invalid input array"));
   
  boost::shared_ptr< cuNDArray<typename realType<T>::Type> > result(new cuNDArray<typename realType<T>::Type>());
  result->create(x->get_dimensions());
  thrust::device_ptr<typename realType<T>::Type> resPtr = result->get_device_ptr();
  thrust::device_ptr<T> xPtr = x->get_device_ptr();
  thrust::transform(xPtr,xPtr+x->get_number_of_elements(),resPtr,cuNDA_abs_square<T>());
  return result;
}

template<typename T> struct cuNDA_sqrt : public thrust::unary_function<T,T>
{
  __device__ T operator()(const T &x) const {return sqrt(x);}
};

template<class T> boost::shared_ptr< cuNDArray<T> > 
Gadgetron::sqrt( cuNDArray<T> *x )
{ 
  if( x == 0x0 )
    BOOST_THROW_EXCEPTION(runtime_error("Gadgetron::sqrt(): Invalid input array"));
   
  boost::shared_ptr< cuNDArray<T> > result(new cuNDArray<T>());
  result->create(x->get_dimensions());
  thrust::device_ptr<T> resPtr = result->get_device_ptr();
  thrust::device_ptr<T> xPtr = x->get_device_ptr();
  thrust::transform(xPtr,xPtr+x->get_number_of_elements(),resPtr,cuNDA_sqrt<T>());
  return result;
}

template<class T> void 
Gadgetron::sqrt_inplace( cuNDArray<T> *x )
{ 
  if( x == 0x0 )
    BOOST_THROW_EXCEPTION(runtime_error("Gadgetron::sqrt_inplace(): Invalid input array"));
   
  thrust::device_ptr<T> xPtr = x->get_device_ptr();
  thrust::transform(xPtr,xPtr+x->get_number_of_elements(),xPtr,cuNDA_sqrt<T>());
}
 
template<typename T> struct cuNDA_square : public thrust::unary_function<T,T>
{
  __device__ T operator()(const T &x) const {return x*x;}
};

template<class T> boost::shared_ptr< cuNDArray<T> > Gadgetron::square( cuNDArray<T> *x )
{ 
  if( x == 0x0 )
    BOOST_THROW_EXCEPTION(runtime_error("Gadgetron::square(): Invalid input array"));
   
  boost::shared_ptr< cuNDArray<T> > result(new cuNDArray<T>());
  result->create(x->get_dimensions());
  thrust::device_ptr<T> resPtr = result->get_device_ptr();
  thrust::device_ptr<T> xPtr = x->get_device_ptr();
  thrust::transform(xPtr,xPtr+x->get_number_of_elements(),resPtr,cuNDA_square<T>());
  return result;
}

template<class T> void 
Gadgetron::square_inplace( cuNDArray<T> *x )
{ 
  if( x == 0x0 )
    BOOST_THROW_EXCEPTION(runtime_error("Gadgetron::square_inplace(): Invalid input array"));
   
  thrust::device_ptr<T> xPtr = x->get_device_ptr();
  thrust::transform(xPtr,xPtr+x->get_number_of_elements(),xPtr,cuNDA_square<T>());
}  

template<typename T> struct cuNDA_reciprocal : public thrust::unary_function<T,T>
{
  __device__ T operator()(const T &x) const {return T(1)/x;}
};

template<class T> boost::shared_ptr< cuNDArray<T> > Gadgetron::reciprocal( cuNDArray<T> *x )
{ 
  if( x == 0x0 )
    BOOST_THROW_EXCEPTION(runtime_error("Gadgetron::reciprocal(): Invalid input array"));
   
  boost::shared_ptr< cuNDArray<T> > result(new cuNDArray<T>());
  result->create(x->get_dimensions());
  thrust::device_ptr<T> resPtr = result->get_device_ptr();
  thrust::device_ptr<T> xPtr = x->get_device_ptr();
  thrust::transform(xPtr,xPtr+x->get_number_of_elements(),resPtr,cuNDA_reciprocal<T>());
  return result;
}

template<class T> void 
Gadgetron::reciprocal_inplace( cuNDArray<T> *x )
{ 
  if( x == 0x0 )
    BOOST_THROW_EXCEPTION(runtime_error("Gadgetron::reciprocal_inplace(): Invalid input array"));
   
  thrust::device_ptr<T> xPtr = x->get_device_ptr();
  thrust::transform(xPtr,xPtr+x->get_number_of_elements(),xPtr,cuNDA_reciprocal<T>());
}  
 
template<typename T> struct cuNDA_reciprocal_sqrt : public thrust::unary_function<T,T>
{
  __device__ T operator()(const T &x) const {return T(1)/sqrt(x);}
};

template<class T> boost::shared_ptr< cuNDArray<T> > Gadgetron::reciprocal_sqrt( cuNDArray<T> *x )
{ 
  if( x == 0x0 )
    BOOST_THROW_EXCEPTION(runtime_error("Gadgetron::reciprocal_sqrt(): Invalid input array"));
   
  boost::shared_ptr< cuNDArray<T> > result(new cuNDArray<T>());
  result->create(x->get_dimensions());
  thrust::device_ptr<T> resPtr = result->get_device_ptr();
  thrust::device_ptr<T> xPtr = x->get_device_ptr();
  thrust::transform(xPtr,xPtr+x->get_number_of_elements(),resPtr,cuNDA_reciprocal_sqrt<T>());
  return result;
}

template<class T> void 
Gadgetron::reciprocal_sqrt_inplace( cuNDArray<T> *x )
{ 
  if( x == 0x0 )
    BOOST_THROW_EXCEPTION(runtime_error("Gadgetron::reciprocal_sqrt_inplace(): Invalid input array"));
   
  thrust::device_ptr<T> xPtr = x->get_device_ptr();
  thrust::transform(xPtr,xPtr+x->get_number_of_elements(),xPtr,cuNDA_reciprocal_sqrt<T>());
}  

template<typename T> struct cuNDA_sgn : public thrust::unary_function<T,T>
{
  __device__ T operator()(const T &x) const {return sgn(x);}
};

template<class T> boost::shared_ptr< cuNDArray<T> > Gadgetron::sgn( cuNDArray<T> *x )
{ 
  if( x == 0x0 )
    BOOST_THROW_EXCEPTION(runtime_error("Gadgetron::sgn(): Invalid input array"));
   
  boost::shared_ptr< cuNDArray<T> > result(new cuNDArray<T>());
  result->create(x->get_dimensions());
  thrust::device_ptr<T> resPtr = result->get_device_ptr();
  thrust::device_ptr<T> xPtr = x->get_device_ptr();
  thrust::transform(xPtr,xPtr+x->get_number_of_elements(),resPtr,cuNDA_sgn<T>());
  return result;
}

template<class T> void 
Gadgetron::sgn_inplace( cuNDArray<T> *x )
{ 
  if( x == 0x0 )
    BOOST_THROW_EXCEPTION(runtime_error("Gadgetron::sgn_inplace(): Invalid input array"));
   
  thrust::device_ptr<T> xPtr = x->get_device_ptr();
  thrust::transform(xPtr,xPtr+x->get_number_of_elements(),xPtr,cuNDA_sgn<T>());
}  
 
template<typename T> struct cuNDA_real : public thrust::unary_function<T,typename realType<T>::Type>
{
  __device__ typename realType<T>::Type operator()(const T &x) const {return real(x);}
};

/*
template<class T> boost::shared_ptr< cuNDArray<T> > 
Gadgetron::real( cuNDArray< std::complex<T> > *x )
{ 
  if( x == 0x0 )
    BOOST_THROW_EXCEPTION(runtime_error("Gadgetron::real(): Invalid input array"));
   
  boost::shared_ptr< cuNDArray<T> > result(new cuNDArray<T>());
  result->create(x->get_dimensions());
  thrust::device_ptr<T> resPtr = result->get_device_ptr();
  thrust::device_ptr< std::complex<T> > xPtr = x->get_device_ptr();
  thrust::transform(xPtr,xPtr+x->get_number_of_elements(),resPtr,cuNDA_real< std::complex<T> >());
  return result;
  }*/

template<class T> boost::shared_ptr< cuNDArray<T> > 
Gadgetron::real( cuNDArray< complext<T> > *x )
{ 
  if( x == 0x0 )
    BOOST_THROW_EXCEPTION(runtime_error("Gadgetron::real(): Invalid input array"));
   
  boost::shared_ptr< cuNDArray<T> > result(new cuNDArray<T>());
  result->create(x->get_dimensions());
  thrust::device_ptr<T> resPtr = result->get_device_ptr();
  thrust::device_ptr< complext<T> > xPtr = x->get_device_ptr();
  thrust::transform(xPtr,xPtr+x->get_number_of_elements(),resPtr,cuNDA_real< complext<T> >());
  return result;
}

template <typename T> struct cuNDA_imag : public thrust::unary_function<T,typename realType<T>::Type>
{
  __device__ typename realType<T>::Type operator()(const T &x) const {return imag(x);}
};

/*
template<class T> boost::shared_ptr< cuNDArray<T> > 
Gadgetron::imag( cuNDArray< std::complex<T> > *x )
{ 
  if( x == 0x0 )
    BOOST_THROW_EXCEPTION(runtime_error("Gadgetron::imag(): Invalid input array"));
   
  boost::shared_ptr< cuNDArray<T> > result(new cuNDArray<T>());
  result->create(x->get_dimensions());
  thrust::device_ptr<T> resPtr = result->get_device_ptr();
  thrust::device_ptr< std::complex<T> > xPtr = x->get_device_ptr();
  thrust::transform(xPtr,xPtr+x->get_number_of_elements(),resPtr,cuNDA_imag< std::complex<T> >());
  return result;
  }*/

template<class T> boost::shared_ptr< cuNDArray<T> > 
Gadgetron::imag( cuNDArray< complext<T> > *x )
{ 
  if( x == 0x0 )
    BOOST_THROW_EXCEPTION(runtime_error("Gadgetron::imag(): Invalid input array"));
   
  boost::shared_ptr< cuNDArray<T> > result(new cuNDArray<T>());
  result->create(x->get_dimensions());
  thrust::device_ptr<T> resPtr = result->get_device_ptr();
  thrust::device_ptr< complext<T> > xPtr = x->get_device_ptr();
  thrust::transform(xPtr,xPtr+x->get_number_of_elements(),resPtr,cuNDA_imag< complext<T> >());
  return result;
}


template <typename T> struct cuNDA_real_to_complex : public thrust::unary_function<typename realType<T>::Type,T>
{
  __device__ T operator()(const typename realType<T>::Type &x) const {return T(x);}
};

/*
template<class T> boost::shared_ptr< cuNDArray< std::complex<T> > > 
Gadgetron::real_to_std_complex( cuNDArray<T> *x )
{
  if( x == 0x0 )
    BOOST_THROW_EXCEPTION(runtime_error("Gadgetron::real_to_std_complex(): Invalid input array"));
   
  boost::shared_ptr< cuNDArray< std::complex<T> > > result(new cuNDArray< std::complex<T> >());
  result->create(x->get_dimensions());
  thrust::device_ptr<std::complex<T> > resPtr = result->get_device_ptr();
  thrust::device_ptr<T> xPtr = x->get_device_ptr();
  thrust::transform(xPtr,xPtr+x->get_number_of_elements(),resPtr,cuNDA_real_to_complex<std::complex<T> >());
  return result;
  }*/

template<class T> boost::shared_ptr< cuNDArray< complext<T> > > 
Gadgetron::real_to_complext( cuNDArray<T> *x )
{
  if( x == 0x0 )
    BOOST_THROW_EXCEPTION(runtime_error("Gadgetron::real_to_complext(): Invalid input array"));
   
  boost::shared_ptr< cuNDArray< complext<T> > > result(new cuNDArray< complext<T> >());
  result->create(x->get_dimensions());
  thrust::device_ptr<complext<T> > resPtr = result->get_device_ptr();
  thrust::device_ptr<T> xPtr = x->get_device_ptr();
  thrust::transform(xPtr,xPtr+x->get_number_of_elements(),resPtr,cuNDA_real_to_complex<complext<T> >());
  return result;
}

template<class T> void Gadgetron::clear( cuNDArray<T> *x )
{
  if( x == 0x0 )
    BOOST_THROW_EXCEPTION(runtime_error("Gadgetron::clear(): Invalid input array"));
  
  cudaMemset(x->get_data_ptr(),0,sizeof(T)*x->get_number_of_elements());
}

template<class T> void 
Gadgetron::fill( cuNDArray<T> *x, T val )
{ 
  if( x == 0x0 )
    BOOST_THROW_EXCEPTION(runtime_error("Gadgetron::fill_inplace(): Invalid input array"));
  
  thrust::device_ptr<T> devPtr = x->get_device_ptr();
  thrust::fill(devPtr,devPtr+x->get_number_of_elements(),val);
}  

template<typename T> struct cuNDA_clamp : public thrust::unary_function<T,T>
{
  cuNDA_clamp( T _min, T _max ) : min(_min), max(_max) {}
  __device__ T operator()(const T &x) const 
  {
    if( x < min ) return min;
    else if ( x > max) return max;
    else return x;
  }
  T min, max;
};

template<class T> void 
Gadgetron::clamp( cuNDArray<T> *x, T min, T max )
{ 
  if( x == 0x0 )
    BOOST_THROW_EXCEPTION(runtime_error("Gadgetron::clamp(): Invalid input array"));
   
  thrust::device_ptr<T> xPtr = x->get_device_ptr();
  thrust::transform(xPtr,xPtr+x->get_number_of_elements(),xPtr,cuNDA_clamp<T>(min, max));
}  

template<typename T> struct cuNDA_clamp_min : public thrust::unary_function<T,T>
{
  cuNDA_clamp_min( T _min ) : min(_min) {}
  __device__ T operator()(const T &x) const 
  {
    if( x < min ) return min;
    else return x;
  }
  T min;
};

template<class T> void 
Gadgetron::clamp_min( cuNDArray<T> *x, T min )
{ 
  if( x == 0x0 )
    BOOST_THROW_EXCEPTION(runtime_error("Gadgetron::clamp_min(): Invalid input array"));
   
  thrust::device_ptr<T> xPtr = x->get_device_ptr();
  thrust::transform(xPtr,xPtr+x->get_number_of_elements(),xPtr,cuNDA_clamp_min<T>(min));
}  

template<typename T> struct cuNDA_clamp_max : public thrust::unary_function<T,T>
{
  cuNDA_clamp_max( T _max ) : max(_max) {}
  __device__ T operator()(const T &x) const 
  {
    if( x > max ) return max;
    else return x;
  }
  T max;
};

template<class T> void 
Gadgetron::clamp_max( cuNDArray<T> *x, T max )
{ 
  if( x == 0x0 )
    BOOST_THROW_EXCEPTION(runtime_error("Gadgetron::clamp_max(): Invalid input array"));
   
  thrust::device_ptr<T> xPtr = x->get_device_ptr();
  thrust::transform(xPtr,xPtr+x->get_number_of_elements(),xPtr,cuNDA_clamp_max<T>(max));
}  

template<class T> void 
Gadgetron::normalize( cuNDArray<T> *x, typename realType<T>::Type val )
{
  if( x == 0x0 )
    BOOST_THROW_EXCEPTION(runtime_error("Gadgetron::normalize(): Invalid input array"));
  
  unsigned int max_idx = amax(x);
  T max_val_before;
  CUDA_CALL(cudaMemcpy(&max_val_before, &x->get_data_ptr()[max_idx], sizeof(T), cudaMemcpyDeviceToHost));
  typename realType<T>::Type scale = val/abs(max_val_before);
  *x *= scale;
}


template<typename T> struct cuNDA_shrink1 : public thrust::unary_function<T,T>
{
  cuNDA_shrink1( typename realType<T>::Type _gamma ) : gamma(_gamma) {}
  __device__ T operator()(const T &x) const {
    typename realType<T>::Type absX = abs(x);
    T sgnX = (absX <= typename realType<T>::Type(0)) ? T(0) : x/absX;
    return sgnX*max(absX-gamma, typename realType<T>::Type(0));
  }
  typename realType<T>::Type gamma;
};

template<class T> void 
Gadgetron::shrink1( cuNDArray<T> *x, typename realType<T>::Type gamma )
{ 
  if( x == 0x0 )
    BOOST_THROW_EXCEPTION(runtime_error("Gadgetron::shrink1(): Invalid input array"));
  
  thrust::device_ptr<T> xPtr = x->get_device_ptr();
  thrust::transform(xPtr,xPtr+x->get_number_of_elements(),xPtr,cuNDA_shrink1<T>(gamma));
}  

template<typename T> struct cuNDA_shrinkd : public thrust::binary_function<T,typename realType<T>::Type,T>
{
  cuNDA_shrinkd( typename realType<T>::Type _gamma ) : gamma(_gamma) {}
  __device__ T operator()(const T &x, const typename realType<T>::Type &s) const {
    return x/s*max(s-gamma,typename realType<T>::Type(0));
  }
  typename realType<T>::Type gamma;
};

template<class T> void 
Gadgetron::shrinkd( cuNDArray<T> *x, cuNDArray<typename realType<T>::Type> *s, typename realType<T>::Type gamma )
{ 
  if( x == 0x0 )
    BOOST_THROW_EXCEPTION(runtime_error("Gadgetron::shrinkd(): Invalid input array"));
  
  thrust::device_ptr<T> xPtr = x->get_device_ptr();
  thrust::device_ptr<typename realType<T>::Type> sPtr = s->get_device_ptr();
  thrust::transform(xPtr,xPtr+x->get_number_of_elements(),sPtr,xPtr,cuNDA_shrinkd<T>(gamma));
}  


//
// Instantiation
//

template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float> > Gadgetron::abs<float>( cuNDArray<float>* );
template EXPORTGPUCORE void Gadgetron::abs_inplace<float>( cuNDArray<float>* );
template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float> > Gadgetron::abs_square<float>( cuNDArray<float>* );
template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float> > Gadgetron::sqrt<float>( cuNDArray<float>* );
template EXPORTGPUCORE void Gadgetron::sqrt_inplace<float>( cuNDArray<float>* );
template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float> > Gadgetron::square<float>( cuNDArray<float>* );
template EXPORTGPUCORE void Gadgetron::square_inplace<float>( cuNDArray<float>* );
template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float> > Gadgetron::reciprocal<float>( cuNDArray<float>* );
template EXPORTGPUCORE void Gadgetron::reciprocal_inplace<float>( cuNDArray<float>* );
template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float> > Gadgetron::reciprocal_sqrt<float>( cuNDArray<float>* );
template EXPORTGPUCORE void Gadgetron::reciprocal_sqrt_inplace<float>( cuNDArray<float>* );
template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float> > Gadgetron::sgn<float>( cuNDArray<float>* );
template EXPORTGPUCORE void Gadgetron::sgn_inplace<float>( cuNDArray<float>* );
template EXPORTGPUCORE void Gadgetron::clear<float>( cuNDArray<float>* );
template EXPORTGPUCORE void Gadgetron::fill<float>( cuNDArray<float>*, float );
template EXPORTGPUCORE void Gadgetron::clamp<float>( cuNDArray<float>*, float, float );
template EXPORTGPUCORE void Gadgetron::clamp_min<float>( cuNDArray<float>*, float );
template EXPORTGPUCORE void Gadgetron::clamp_max<float>( cuNDArray<float>*, float );
template EXPORTGPUCORE void Gadgetron::normalize<float>( cuNDArray<float>*, float );
template EXPORTGPUCORE void Gadgetron::shrink1<float>( cuNDArray<float>*, float );
template EXPORTGPUCORE void Gadgetron::shrinkd<float> ( cuNDArray<float>*, cuNDArray<float>*, float );

template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double> > Gadgetron::abs<double>( cuNDArray<double>* );
template EXPORTGPUCORE void Gadgetron::abs_inplace<double>( cuNDArray<double>* );
template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double> > Gadgetron::abs_square<double>( cuNDArray<double>* );
template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double> > Gadgetron::sqrt<double>( cuNDArray<double>* );
template EXPORTGPUCORE void Gadgetron::sqrt_inplace<double>( cuNDArray<double>* );
template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double> > Gadgetron::square<double>( cuNDArray<double>* );
template EXPORTGPUCORE void Gadgetron::square_inplace<double>( cuNDArray<double>* );
template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double> > Gadgetron::reciprocal<double>( cuNDArray<double>* );
template EXPORTGPUCORE void Gadgetron::reciprocal_inplace<double>( cuNDArray<double>* );
template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double> > Gadgetron::reciprocal_sqrt<double>( cuNDArray<double>* );
template EXPORTGPUCORE void Gadgetron::reciprocal_sqrt_inplace<double>( cuNDArray<double>* );
template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double> > Gadgetron::sgn<double>( cuNDArray<double>* );
template EXPORTGPUCORE void Gadgetron::sgn_inplace<double>( cuNDArray<double>* );
template EXPORTGPUCORE void Gadgetron::clear<double>( cuNDArray<double>* );
template EXPORTGPUCORE void Gadgetron::fill<double>( cuNDArray<double>*, double );
template EXPORTGPUCORE void Gadgetron::clamp<double>( cuNDArray<double>*, double, double );
template EXPORTGPUCORE void Gadgetron::clamp_min<double>( cuNDArray<double>*, double );
template EXPORTGPUCORE void Gadgetron::clamp_max<double>( cuNDArray<double>*, double );
template EXPORTGPUCORE void Gadgetron::normalize<double>( cuNDArray<double>*, double );
template EXPORTGPUCORE void Gadgetron::shrink1<double>( cuNDArray<double>*, double );
template EXPORTGPUCORE void Gadgetron::shrinkd<double> ( cuNDArray<double>*, cuNDArray<double>*, double );

/*template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float> > Gadgetron::abs< std::complex<float> >( cuNDArray< std::complex<float> >* );
template EXPORTGPUCORE boost::shared_ptr< cuNDArray< std::complex<float> > > Gadgetron::sqrt< std::complex<float> >( cuNDArray< std::complex<float> >* );
template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float> > Gadgetron::abs_square< std::complex<float> >( cuNDArray< std::complex<float> >* );
template EXPORTGPUCORE void Gadgetron::sqrt_inplace< std::complex<float> >( cuNDArray< std::complex<float> >* );
template EXPORTGPUCORE boost::shared_ptr< cuNDArray< std::complex<float> > > Gadgetron::square< std::complex<float> >( cuNDArray< std::complex<float> >* );
template EXPORTGPUCORE void Gadgetron::square_inplace< std::complex<float> >( cuNDArray< std::complex<float> >* );
template EXPORTGPUCORE boost::shared_ptr< cuNDArray< std::complex<float> > > Gadgetron::reciprocal< std::complex<float> >( cuNDArray< std::complex<float> >* );
template EXPORTGPUCORE void Gadgetron::reciprocal_inplace< std::complex<float> >( cuNDArray< std::complex<float> >* );
template EXPORTGPUCORE boost::shared_ptr< cuNDArray< std::complex<float> > > Gadgetron::reciprocal_sqrt< std::complex<float> >( cuNDArray< std::complex<float> >* );
template EXPORTGPUCORE void Gadgetron::reciprocal_sqrt_inplace< std::complex<float> >( cuNDArray< std::complex<float> >* );
template EXPORTGPUCORE void Gadgetron::clear< std::complex<float> >( cuNDArray< std::complex<float> >* );
template EXPORTGPUCORE void Gadgetron::fill< std::complex<float> >( cuNDArray< std::complex<float> >*, std::complex<float> );
template EXPORTGPUCORE void Gadgetron::normalize< std::complex<float> >( cuNDArray< std::complex<float> >*, float );
template EXPORTGPUCORE void Gadgetron::shrink1< std::complex<float> >( cuNDArray< std::complex<float> >*, float );
template EXPORTGPUCORE void Gadgetron::shrinkd< std::complex<float> > ( cuNDArray< std::complex<float> >*, cuNDArray<float>*, float );

template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double> > Gadgetron::abs< std::complex<double> >( cuNDArray< std::complex<double> >* );
template EXPORTGPUCORE boost::shared_ptr< cuNDArray< std::complex<double> > > Gadgetron::sqrt< std::complex<double> >( cuNDArray< std::complex<double> >* );
template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double> > Gadgetron::abs_square< std::complex<double> >( cuNDArray< std::complex<double> >* );
template EXPORTGPUCORE void Gadgetron::sqrt_inplace< std::complex<double> >( cuNDArray< std::complex<double> >* );
template EXPORTGPUCORE boost::shared_ptr< cuNDArray< std::complex<double> > > Gadgetron::square< std::complex<double> >( cuNDArray< std::complex<double> >* );
template EXPORTGPUCORE void Gadgetron::square_inplace< std::complex<double> >( cuNDArray< std::complex<double> >* );
template EXPORTGPUCORE boost::shared_ptr< cuNDArray< std::complex<double> > > Gadgetron::reciprocal< std::complex<double> >( cuNDArray< std::complex<double> >* );
template EXPORTGPUCORE void Gadgetron::reciprocal_inplace< std::complex<double> >( cuNDArray< std::complex<double> >* );
template EXPORTGPUCORE boost::shared_ptr< cuNDArray< std::complex<double> > > Gadgetron::reciprocal_sqrt< std::complex<double> >( cuNDArray< std::complex<double> >* );
template EXPORTGPUCORE void Gadgetron::reciprocal_sqrt_inplace< std::complex<double> >( cuNDArray< std::complex<double> >* );
template EXPORTGPUCORE void Gadgetron::clear< std::complex<double> >( cuNDArray< std::complex<double> >* );
template EXPORTGPUCORE void Gadgetron::fill< std::complex<double> >( cuNDArray< std::complex<double> >*, std::complex<double> );
template EXPORTGPUCORE void Gadgetron::normalize< std::complex<double> >( cuNDArray< std::complex<double> >*, double );
template EXPORTGPUCORE void Gadgetron::shrink1< std::complex<double> >( cuNDArray< std::complex<double> >*, double );
template EXPORTGPUCORE void Gadgetron::shrinkd< std::complex<double> > ( cuNDArray< std::complex<double> >*, cuNDArray<double>*, double );
*/
template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float> > Gadgetron::abs< complext<float> >( cuNDArray< complext<float> >* );
template EXPORTGPUCORE boost::shared_ptr< cuNDArray< complext<float> > > Gadgetron::sqrt< complext<float> >( cuNDArray< complext<float> >* );
template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float> > Gadgetron::abs_square< complext<float> >( cuNDArray< complext<float> >* );
template EXPORTGPUCORE void Gadgetron::sqrt_inplace< complext<float> >( cuNDArray< complext<float> >* );
template EXPORTGPUCORE boost::shared_ptr< cuNDArray< complext<float> > > Gadgetron::square< complext<float> >( cuNDArray< complext<float> >* );
template EXPORTGPUCORE void Gadgetron::square_inplace< complext<float> >( cuNDArray< complext<float> >* );
template EXPORTGPUCORE boost::shared_ptr< cuNDArray< complext<float> > > Gadgetron::reciprocal< complext<float> >( cuNDArray< complext<float> >* );
template EXPORTGPUCORE void Gadgetron::reciprocal_inplace< complext<float> >( cuNDArray< complext<float> >* );
template EXPORTGPUCORE boost::shared_ptr< cuNDArray< complext<float> > > Gadgetron::reciprocal_sqrt< complext<float> >( cuNDArray< complext<float> >* );
template EXPORTGPUCORE void Gadgetron::reciprocal_sqrt_inplace< complext<float> >( cuNDArray< complext<float> >* );
template EXPORTGPUCORE void Gadgetron::clear< complext<float> >( cuNDArray< complext<float> >* );
template EXPORTGPUCORE void Gadgetron::fill< complext<float> >( cuNDArray< complext<float> >*, complext<float> );
template EXPORTGPUCORE void Gadgetron::normalize< complext<float> >( cuNDArray< complext<float> >*, float );
template EXPORTGPUCORE void Gadgetron::shrink1< complext<float> >( cuNDArray< complext<float> >*, float );
template EXPORTGPUCORE void Gadgetron::shrinkd< complext<float> > ( cuNDArray< complext<float> >*, cuNDArray<float>*, float );

template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double> > Gadgetron::abs< complext<double> >( cuNDArray< complext<double> >* );
template EXPORTGPUCORE boost::shared_ptr< cuNDArray< complext<double> > > Gadgetron::sqrt< complext<double> >( cuNDArray< complext<double> >* );
template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double> > Gadgetron::abs_square< complext<double> >( cuNDArray< complext<double> >* );
template EXPORTGPUCORE void Gadgetron::sqrt_inplace< complext<double> >( cuNDArray< complext<double> >* );
template EXPORTGPUCORE boost::shared_ptr< cuNDArray< complext<double> > > Gadgetron::square< complext<double> >( cuNDArray< complext<double> >* );
template EXPORTGPUCORE void Gadgetron::square_inplace< complext<double> >( cuNDArray< complext<double> >* );
template EXPORTGPUCORE boost::shared_ptr< cuNDArray< complext<double> > > Gadgetron::reciprocal< complext<double> >( cuNDArray< complext<double> >* );
template EXPORTGPUCORE void Gadgetron::reciprocal_inplace< complext<double> >( cuNDArray< complext<double> >* );
template EXPORTGPUCORE boost::shared_ptr< cuNDArray< complext<double> > > Gadgetron::reciprocal_sqrt< complext<double> >( cuNDArray< complext<double> >* );
template EXPORTGPUCORE void Gadgetron::reciprocal_sqrt_inplace< complext<double> >( cuNDArray< complext<double> >* );
template EXPORTGPUCORE void Gadgetron::clear< complext<double> >( cuNDArray< complext<double> >* );
template EXPORTGPUCORE void Gadgetron::fill< complext<double> >( cuNDArray< complext<double> >*, complext<double> );
template EXPORTGPUCORE void Gadgetron::normalize< complext<double> >( cuNDArray< complext<double> >*, double );
template EXPORTGPUCORE void Gadgetron::shrink1< complext<double> >( cuNDArray< complext<double> >*, double );
template EXPORTGPUCORE void Gadgetron::shrinkd< complext<double> > ( cuNDArray< complext<double> >*, cuNDArray<double>*, double );

//template EXPORTGPUCORE boost::shared_ptr< cuNDArray< std::complex<float> > > Gadgetron::real_to_std_complex<float>( cuNDArray<float>* );
template EXPORTGPUCORE boost::shared_ptr< cuNDArray< complext<float> > > Gadgetron::real_to_complext<float>( cuNDArray<float>* );
//template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float> > Gadgetron::real<float>( cuNDArray< std::complex<float> >* );
template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float> > Gadgetron::real<float>( cuNDArray< complext<float> >* );
//template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float> > Gadgetron::imag<float>( cuNDArray< std::complex<float> >* );
template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float> > Gadgetron::imag<float>( cuNDArray< complext<float> >* );

//template EXPORTGPUCORE boost::shared_ptr< cuNDArray< std::complex<double> > > Gadgetron::real_to_std_complex<double>( cuNDArray<double>* );
template EXPORTGPUCORE boost::shared_ptr< cuNDArray< complext<double> > > Gadgetron::real_to_complext<double>( cuNDArray<double>* );
//template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double> > Gadgetron::real<double>( cuNDArray< std::complex<double> >* );
template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double> > Gadgetron::real<double>( cuNDArray< complext<double> >* );
//template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double> > Gadgetron::imag<double>( cuNDArray< std::complex<double> >* );
template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double> > Gadgetron::imag<double>( cuNDArray< complext<double> >* );
