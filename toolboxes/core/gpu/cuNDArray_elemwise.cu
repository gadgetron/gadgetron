#include "cuNDArray_elemwise.h"
#include "cuNDArray_operators.h"
#include "cuNDArray_blas.h"
#include "complext.h"

#include <complex>
#include <thrust/functional.h>

using namespace Gadgetron;
//using namespace std;

template<typename T> struct cuNDA_abs : public thrust::unary_function<T,typename realType<T>::Type>
{
  __device__ typename Gadgetron::realType<T>::Type operator()(const T &x) const {return abs(x);}
};

template<class T> boost::shared_ptr< cuNDArray<typename realType<T>::Type> > 
Gadgetron::abs( const cuNDArray<T> *x )
{ 
  if( x == 0x0 )
    throw std::runtime_error("Gadgetron::abs(): Invalid input array");
   
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
    throw std::runtime_error("Gadgetron::abs_inplace(): Invalid input array");
   
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
Gadgetron::abs_square( const cuNDArray<T> *x )
{ 
  if( x == 0x0 )
    throw std::runtime_error("Gadgetron::abs_square(): Invalid input array");
   
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
Gadgetron::sqrt( const cuNDArray<T> *x )
{ 
  if( x == 0x0 )
    throw std::runtime_error("Gadgetron::sqrt(): Invalid input array");
   
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
    throw std::runtime_error("Gadgetron::sqrt_inplace(): Invalid input array");
   
  thrust::device_ptr<T> xPtr = x->get_device_ptr();
  thrust::transform(xPtr,xPtr+x->get_number_of_elements(),xPtr,cuNDA_sqrt<T>());
}
 
template<typename T> struct cuNDA_square : public thrust::unary_function<T,T>
{
  __device__ T operator()(const T &x) const {return x*x;}
};

template<class T> boost::shared_ptr< cuNDArray<T> > Gadgetron::square( const cuNDArray<T> *x )
{ 
  if( x == 0x0 )
    throw std::runtime_error("Gadgetron::square(): Invalid input array");
   
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
    throw std::runtime_error("Gadgetron::square_inplace(): Invalid input array");
   
  thrust::device_ptr<T> xPtr = x->get_device_ptr();
  thrust::transform(xPtr,xPtr+x->get_number_of_elements(),xPtr,cuNDA_square<T>());
}  

template<typename T> struct cuNDA_reciprocal : public thrust::unary_function<T,T>
{
  __device__ T operator()(const T &x) const {return T(1)/x;}
};

template<class T> boost::shared_ptr< cuNDArray<T> > Gadgetron::reciprocal( const cuNDArray<T> *x )
{ 
  if( x == 0x0 )
    throw std::runtime_error("Gadgetron::reciprocal(): Invalid input array");
   
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
    throw std::runtime_error("Gadgetron::reciprocal_inplace(): Invalid input array");
   
  thrust::device_ptr<T> xPtr = x->get_device_ptr();
  thrust::transform(xPtr,xPtr+x->get_number_of_elements(),xPtr,cuNDA_reciprocal<T>());
}  
 
template<typename T> struct cuNDA_reciprocal_sqrt : public thrust::unary_function<T,T>
{
  __device__ T operator()(const T &x) const {return T(1)/sqrt(x);}
};

template<class T> boost::shared_ptr< cuNDArray<T> > Gadgetron::reciprocal_sqrt( const cuNDArray<T> *x )
{ 
  if( x == 0x0 )
    throw std::runtime_error("Gadgetron::reciprocal_sqrt(): Invalid input array");
   
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
    throw std::runtime_error("Gadgetron::reciprocal_sqrt_inplace(): Invalid input array");
   
  thrust::device_ptr<T> xPtr = x->get_device_ptr();
  thrust::transform(xPtr,xPtr+x->get_number_of_elements(),xPtr,cuNDA_reciprocal_sqrt<T>());
}  

template<typename T> struct cuNDA_sgn : public thrust::unary_function<T,T>
{
  __device__ T operator()(const T &x) const {return sgn(x);}
};

template<class T> boost::shared_ptr< cuNDArray<T> > Gadgetron::sgn( const cuNDArray<T> *x )
{ 
  if( x == 0x0 )
    throw std::runtime_error("Gadgetron::sgn(): Invalid input array");
   
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
    throw std::runtime_error("Gadgetron::sgn_inplace(): Invalid input array");
   
  thrust::device_ptr<T> xPtr = x->get_device_ptr();
  thrust::transform(xPtr,xPtr+x->get_number_of_elements(),xPtr,cuNDA_sgn<T>());
}  
 
template<typename T> struct cuNDA_real : public thrust::unary_function<T,typename realType<T>::Type>
{
  __device__ typename realType<T>::Type operator()(const T &x) const {return real(x);}
};

template<class T> boost::shared_ptr< cuNDArray<typename realType<T>::Type> > 
Gadgetron::real( const cuNDArray<T> *x )
{ 
  if( x == 0x0 )
    throw std::runtime_error("Gadgetron::real(): Invalid input array");
   
  boost::shared_ptr< cuNDArray<typename realType<T>::Type> > result(new cuNDArray<typename realType<T>::Type>());
  result->create(x->get_dimensions());
  thrust::device_ptr<typename realType<T>::Type> resPtr = result->get_device_ptr();
  thrust::device_ptr<T> xPtr = x->get_device_ptr();
  thrust::transform(xPtr,xPtr+x->get_number_of_elements(),resPtr,cuNDA_real<T>());
  return result;
}

template <typename T> struct cuNDA_imag : public thrust::unary_function<T,typename realType<T>::Type>
{
  __device__ typename realType<T>::Type operator()(const T &x) const {return imag(x);}
};

template<class T> boost::shared_ptr< cuNDArray<typename realType<T>::Type> > 
Gadgetron::imag( const cuNDArray<T> *x )
{ 
  if( x == 0x0 )
    throw std::runtime_error("Gadgetron::imag(): Invalid input array");
   
  boost::shared_ptr< cuNDArray<typename realType<T>::Type> > result(new cuNDArray<typename realType<T>::Type>());
  result->create(x->get_dimensions());
  thrust::device_ptr<typename realType<T>::Type> resPtr = result->get_device_ptr();
  thrust::device_ptr<T> xPtr = x->get_device_ptr();
  thrust::transform(xPtr,xPtr+x->get_number_of_elements(),resPtr,cuNDA_imag<T>());
  return result;
}

template <typename T> struct cuNDA_conj : public thrust::unary_function<T,T>
{
  __device__ T operator()(const T &x) const {return conj(x);}
};

template<class T> boost::shared_ptr< cuNDArray<T> > 
Gadgetron::conj( const cuNDArray<T> *x )
{ 
  if( x == 0x0 )
    throw std::runtime_error("Gadgetron::conj(): Invalid input array");
   
  boost::shared_ptr< cuNDArray<T> > result(new cuNDArray<T>());
  result->create(x->get_dimensions());
  thrust::device_ptr<T> resPtr = result->get_device_ptr();
  thrust::device_ptr<T> xPtr = x->get_device_ptr();
  thrust::transform(xPtr,xPtr+x->get_number_of_elements(),resPtr,cuNDA_conj<T>());
  return result;
}

template <typename T> struct cuNDA_real_to_complex : public thrust::unary_function<typename realType<T>::Type,T>
{
  __device__ T operator()(const typename realType<T>::Type &x) const {return T(x);}
};

template<class T> boost::shared_ptr< cuNDArray<T> > 
Gadgetron::real_to_complex( const cuNDArray<typename realType<T>::Type> *x )
{
  if( x == 0x0 )
    throw std::runtime_error("Gadgetron::real_to_complex(): Invalid input array");
   
  boost::shared_ptr< cuNDArray<T> > result(new cuNDArray<T>());
  result->create(x->get_dimensions());
  thrust::device_ptr<T> resPtr = result->get_device_ptr();
  thrust::device_ptr<typename realType<T>::Type> xPtr = x->get_device_ptr();
  thrust::transform(xPtr,xPtr+x->get_number_of_elements(),resPtr,cuNDA_real_to_complex<T>());
  return result;
}

template <typename T,typename T2> struct cuNDA_convert_to : public thrust::unary_function<T,T2>
{
  __device__ T2 operator()(T &x) const {return T2(x);}
};

template <typename T,typename T2> struct cuNDA_convert_to<complext<T>,complext<T2> > : public thrust::unary_function<complext<T>,complext<T2> >
{
  __device__ complext<T2> operator()(complext<T> &x) const {return complext<T2>(x._real,x._imag);}
};

template<class T, class T2> boost::shared_ptr< cuNDArray<T2> >
Gadgetron::convert_to( const cuNDArray<T> *x )
{
  if( x == 0x0 )
    throw std::runtime_error("Gadgetron::convert_to(): Invalid input array");

  boost::shared_ptr< cuNDArray<T2> > result(new cuNDArray<T2>());
  result->create(x->get_dimensions());
  thrust::device_ptr<T2> resPtr = result->get_device_ptr();
  thrust::device_ptr<T> xPtr = x->get_device_ptr();
  thrust::transform(xPtr,xPtr+x->get_number_of_elements(),resPtr,cuNDA_convert_to<T,T2>());
  return result;
}

template<class T, class T2> void
Gadgetron::convert_to( cuNDArray<T> *x ,cuNDArray<T2> * y)
{
  if( x == 0x0 || !x->dimensions_equal(y))
    throw std::runtime_error("Gadgetron::convert_to(): Invalid input array");
  thrust::device_ptr<T2> resPtr = y->get_device_ptr();
  thrust::device_ptr<T> xPtr = x->get_device_ptr();
  thrust::transform(xPtr,xPtr+x->get_number_of_elements(),resPtr,cuNDA_convert_to<T,T2>());
}

template<class T> void Gadgetron::clear( cuNDArray<T> *x )
{
  if( x == 0x0 )
    throw std::runtime_error("Gadgetron::clear(): Invalid input array");

  if ( x->get_number_of_elements() > 0 )
  {
    cudaMemset(x->get_data_ptr(),0,sizeof(T)*x->get_number_of_elements());
  }
}

template<class T> void 
Gadgetron::fill( cuNDArray<T> *x, T val )
{ 
  if( x == 0x0 )
    throw std::runtime_error("Gadgetron::fill(): Invalid input array");
  
  thrust::device_ptr<T> devPtr = x->get_device_ptr();
  thrust::fill(devPtr,devPtr+x->get_number_of_elements(),val);
}  

template<typename T> struct cuNDA_clamp : public thrust::unary_function<T,T>
{
  cuNDA_clamp( T _min, T _max, T _min_val, T _max_val ) : min(_min), max(_max),min_val(_min_val), max_val(_max_val) {}
  __device__ T operator()(const T &x) const 
  {
    if( x < min ) return min_val;
    else if ( x >= max) return max_val;
    else return x;
  }
  T min, max;
  T min_val, max_val;
};

template<typename T> struct cuNDA_clamp< complext<T> > : public thrust::unary_function< complext<T>, complext<T> >
{
	cuNDA_clamp( T _min, T _max, complext<T> _min_val, complext<T> _max_val ) : min(_min), max(_max),min_val(_min_val), max_val(_max_val) {}
  __device__ complext<T> operator()(const complext<T> &x) const 
  {
    if( real(x) < min ) return min_val;
    else if ( real(x) >= max) return max_val;
    else return complext<T>(real(x));
  }
  T min, max;
  complext<T> min_val, max_val;
};

template<class T> void 
Gadgetron::clamp( cuNDArray<T> *x, typename realType<T>::Type min, typename realType<T>::Type max, T min_val, T max_val)
{ 
  if( x == 0x0 )
    throw std::runtime_error("Gadgetron::clamp(): Invalid input array");
   
  thrust::device_ptr<T> xPtr = x->get_device_ptr();
  thrust::transform(xPtr,xPtr+x->get_number_of_elements(),xPtr,cuNDA_clamp<T>(min, max,min_val, max_val));
}  

template<class T> void 
Gadgetron::clamp( cuNDArray<T> *x, typename realType<T>::Type min, typename realType<T>::Type max)
{
    clamp(x,min,max,T(min),T(max));
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

template<typename T> struct cuNDA_clamp_min< complext<T> > : public thrust::unary_function< complext<T>, complext<T> >
{
  cuNDA_clamp_min( T _min ) : min(_min) {}
  __device__ complext<T> operator()(const complext<T> &x) const 
  {
    if( real(x) < min ) return complext<T>(min);
    else return complext<T>(real(x));
  }
  T min;
};

template<class T> void 
Gadgetron::clamp_min( cuNDArray<T> *x, typename realType<T>::Type min )
{ 
  if( x == 0x0 )
    throw std::runtime_error("Gadgetron::clamp_min(): Invalid input array");
   
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

template<typename T> struct cuNDA_clamp_max< complext<T> > : public thrust::unary_function< complext<T>, complext<T> >
{
  cuNDA_clamp_max( T _max ) : max(_max) {}
  __device__ complext<T> operator()(const complext<T> &x) const 
  {
    if( real(x) > max ) return complext<T>(max);
    else return complext<T>(real(x));
  }
  T max;
};

template<class T> void 
Gadgetron::clamp_max( cuNDArray<T> *x, typename realType<T>::Type max )
{ 
  if( x == 0x0 )
    throw std::runtime_error("Gadgetron::clamp_max(): Invalid input array");
   
  thrust::device_ptr<T> xPtr = x->get_device_ptr();
  thrust::transform(xPtr,xPtr+x->get_number_of_elements(),xPtr,cuNDA_clamp_max<T>(max));
}  

template<class T> void 
Gadgetron::normalize( cuNDArray<T> *x, typename realType<T>::Type val )
{
  if( x == 0x0 )
    throw std::runtime_error("Gadgetron::normalize(): Invalid input array");
  
  size_t max_idx = amax(x);
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
Gadgetron::shrink1( cuNDArray<T> *x, typename realType<T>::Type gamma, cuNDArray<T> *out )
{ 
  if( x == 0x0 )
    throw std::runtime_error("Gadgetron::shrink1(): Invalid input array");
  
  thrust::device_ptr<T> xPtr = x->get_device_ptr();
  thrust::device_ptr<T> outPtr = (out == 0x0) ? x->get_device_ptr() : out->get_device_ptr();
  thrust::transform(xPtr,xPtr+x->get_number_of_elements(),outPtr,cuNDA_shrink1<T>(gamma));
}

template<typename T> struct cuNDA_pshrink : public thrust::unary_function<T,T>
{
  cuNDA_pshrink( typename realType<T>::Type _gamma, typename realType<T>::Type _p ) : gamma(_gamma),p(_p) {}
  __device__ T operator()(const T &x) const {
    typename realType<T>::Type absX = abs(x);
    T sgnX = (absX <= typename realType<T>::Type(0)) ? T(0) : x/absX;
    return sgnX*max(absX-gamma*pow(absX,p-1), typename realType<T>::Type(0));
  }
  typename realType<T>::Type gamma;
  typename realType<T>::Type p;
};

template<class T> void
Gadgetron::pshrink( cuNDArray<T> *x, typename realType<T>::Type gamma,typename realType<T>::Type p, cuNDArray<T> *out )
{
  if( x == 0x0 )
    throw std::runtime_error("Gadgetron::shrink1(): Invalid input array");

  thrust::device_ptr<T> xPtr = x->get_device_ptr();
  thrust::device_ptr<T> outPtr = (out == 0x0) ? x->get_device_ptr() : out->get_device_ptr();
  thrust::transform(xPtr,xPtr+x->get_number_of_elements(),outPtr,cuNDA_pshrink<T>(gamma,p));
}  

template<typename T> struct cuNDA_shrinkd : public thrust::binary_function<T,typename realType<T>::Type,T>
{
  cuNDA_shrinkd( typename realType<T>::Type _gamma ) : gamma(_gamma) {}
  __device__ T operator()(const T &x, const typename realType<T>::Type &s) const {
  	T xs = (s <= typename realType<T>::Type(0)) ? T(0) : x/s;
    return xs*max(s-gamma,typename realType<T>::Type(0));
  }
  typename realType<T>::Type gamma;
};

template<class T> void 
Gadgetron::shrinkd( cuNDArray<T> *x, cuNDArray<typename realType<T>::Type> *s, typename realType<T>::Type gamma, cuNDArray<T> *out )
{ 
  if( x == 0x0 )
    throw std::runtime_error("Gadgetron::shrinkd(): Invalid input array");
  
  thrust::device_ptr<T> xPtr = x->get_device_ptr();
  thrust::device_ptr<T> outPtr = (out == 0x0) ? x->get_device_ptr() : out->get_device_ptr();
  thrust::device_ptr<typename realType<T>::Type> sPtr = s->get_device_ptr();
  thrust::transform(xPtr,xPtr+x->get_number_of_elements(),sPtr,outPtr,cuNDA_shrinkd<T>(gamma));
}  


template<typename T> struct cuNDA_pshrinkd : public thrust::binary_function<T,typename realType<T>::Type,T>
{
  cuNDA_pshrinkd( typename realType<T>::Type _gamma,typename realType<T>::Type _p ) : gamma(_gamma), p(_p) {}
  __device__ T operator()(const T &x, const typename realType<T>::Type &s) const {
    return x/s*max(s-gamma*pow(s,p-1),typename realType<T>::Type(0));
  }
  typename realType<T>::Type gamma;
  typename realType<T>::Type p;
};

template<class T> void
Gadgetron::pshrinkd( cuNDArray<T> *x, cuNDArray<typename realType<T>::Type> *s, typename realType<T>::Type gamma,typename realType<T>::Type p, cuNDArray<T> *out )
{
  if( x == 0x0 )
    throw std::runtime_error("Gadgetron::shrinkd(): Invalid input array");

  thrust::device_ptr<T> xPtr = x->get_device_ptr();
  thrust::device_ptr<T> outPtr = (out == 0x0) ? x->get_device_ptr() : out->get_device_ptr();
  thrust::device_ptr<typename realType<T>::Type> sPtr = s->get_device_ptr();
  thrust::transform(xPtr,xPtr+x->get_number_of_elements(),sPtr,outPtr,cuNDA_pshrinkd<T>(gamma,p));
}

//
// Instantiation
//

template boost::shared_ptr< cuNDArray<float> > Gadgetron::abs<float>( const cuNDArray<float>* );
template void Gadgetron::abs_inplace<float>( cuNDArray<float>* );
template boost::shared_ptr< cuNDArray<float> > Gadgetron::abs_square<float>( const cuNDArray<float>* );
template boost::shared_ptr< cuNDArray<float> > Gadgetron::sqrt<float>( const cuNDArray<float>* );
template void Gadgetron::sqrt_inplace<float>( cuNDArray<float>* );
template boost::shared_ptr< cuNDArray<float> > Gadgetron::square<float>( const cuNDArray<float>* );
template void Gadgetron::square_inplace<float>( cuNDArray<float>* );
template boost::shared_ptr< cuNDArray<float> > Gadgetron::reciprocal<float>( const cuNDArray<float>* );
template void Gadgetron::reciprocal_inplace<float>( cuNDArray<float>* );
template boost::shared_ptr< cuNDArray<float> > Gadgetron::reciprocal_sqrt<float>( const cuNDArray<float>* );
template void Gadgetron::reciprocal_sqrt_inplace<float>( cuNDArray<float>* );
template boost::shared_ptr< cuNDArray<float> > Gadgetron::sgn<float>( const cuNDArray<float>* );
template void Gadgetron::sgn_inplace<float>( cuNDArray<float>* );
template void Gadgetron::clear<float>( cuNDArray<float>* );
template void Gadgetron::fill<float>( cuNDArray<float>*, float );
template void Gadgetron::clamp<float>( cuNDArray<float>*, float, float );
template void Gadgetron::clamp<float>( cuNDArray<float>*, float, float, float,float );
template void Gadgetron::clamp_min<float>( cuNDArray<float>*, float );
template void Gadgetron::clamp_max<float>( cuNDArray<float>*, float );
template void Gadgetron::normalize<float>( cuNDArray<float>*, float );
template void Gadgetron::shrink1<float>( cuNDArray<float>*, float, cuNDArray<float>* );
template void Gadgetron::pshrink<float>( cuNDArray<float>*, float,float, cuNDArray<float>* );
template void Gadgetron::shrinkd<float> ( cuNDArray<float>*, cuNDArray<float>*, float, cuNDArray<float>* );
template void Gadgetron::pshrinkd<float> ( cuNDArray<float>*, cuNDArray<float>*, float,float, cuNDArray<float>* );

template boost::shared_ptr< cuNDArray<double> > Gadgetron::abs<double>( const cuNDArray<double>* );
template void Gadgetron::abs_inplace<double>( cuNDArray<double>* );
template boost::shared_ptr< cuNDArray<double> > Gadgetron::abs_square<double>( const cuNDArray<double>* );
template boost::shared_ptr< cuNDArray<double> > Gadgetron::sqrt<double>( const cuNDArray<double>* );
template void Gadgetron::sqrt_inplace<double>( cuNDArray<double>* );
template boost::shared_ptr< cuNDArray<double> > Gadgetron::square<double>( const cuNDArray<double>* );
template void Gadgetron::square_inplace<double>( cuNDArray<double>* );
template boost::shared_ptr< cuNDArray<double> > Gadgetron::reciprocal<double>( const cuNDArray<double>* );
template void Gadgetron::reciprocal_inplace<double>( cuNDArray<double>* );
template boost::shared_ptr< cuNDArray<double> > Gadgetron::reciprocal_sqrt<double>( const cuNDArray<double>* );
template void Gadgetron::reciprocal_sqrt_inplace<double>( cuNDArray<double>* );
template boost::shared_ptr< cuNDArray<double> > Gadgetron::sgn<double>( const cuNDArray<double>* );
template void Gadgetron::sgn_inplace<double>( cuNDArray<double>* );
template void Gadgetron::clear<double>( cuNDArray<double>* );
template void Gadgetron::fill<double>( cuNDArray<double>*, double );
template void Gadgetron::clamp<double>( cuNDArray<double>*, double, double );
template void Gadgetron::clamp<double>( cuNDArray<double>*, double, double, double, double );
template void Gadgetron::clamp_min<double>( cuNDArray<double>*, double );
template void Gadgetron::clamp_max<double>( cuNDArray<double>*, double );
template void Gadgetron::normalize<double>( cuNDArray<double>*, double );
template void Gadgetron::shrink1<double>( cuNDArray<double>*, double, cuNDArray<double>* );
template void Gadgetron::pshrink<double>( cuNDArray<double>*, double,double, cuNDArray<double>* );
template void Gadgetron::shrinkd<double> ( cuNDArray<double>*, cuNDArray<double>*, double, cuNDArray<double>* );
template void Gadgetron::pshrinkd<double> ( cuNDArray<double>*, cuNDArray<double>*, double,double, cuNDArray<double>* );


template void Gadgetron::fill<bool>( cuNDArray<bool>*, bool );
/*template boost::shared_ptr< cuNDArray<float> > Gadgetron::abs< std::complex<float> >( const cuNDArray< std::complex<float> >* );
template boost::shared_ptr< cuNDArray< std::complex<float> > > Gadgetron::sqrt< std::complex<float> >( const cuNDArray< std::complex<float> >* );
template boost::shared_ptr< cuNDArray<float> > Gadgetron::abs_square< std::complex<float> >( const cuNDArray< std::complex<float> >* );
template void Gadgetron::sqrt_inplace< std::complex<float> >( cuNDArray< std::complex<float> >* );
template boost::shared_ptr< cuNDArray< std::complex<float> > > Gadgetron::square< std::complex<float> >( const cuNDArray< std::complex<float> >* );
template void Gadgetron::square_inplace< std::complex<float> >( cuNDArray< std::complex<float> >* );
template boost::shared_ptr< cuNDArray< std::complex<float> > > Gadgetron::reciprocal< std::complex<float> >( const cuNDArray< std::complex<float> >* );
template void Gadgetron::reciprocal_inplace< std::complex<float> >( cuNDArray< std::complex<float> >* );
template boost::shared_ptr< cuNDArray< std::complex<float> > > Gadgetron::reciprocal_sqrt< std::complex<float> >( const cuNDArray< std::complex<float> >* );
template void Gadgetron::reciprocal_sqrt_inplace< std::complex<float> >( cuNDArray< std::complex<float> >* );
template void Gadgetron::clear< std::complex<float> >( cuNDArray< std::complex<float> >* );
template void Gadgetron::fill< std::complex<float> >( cuNDArray< std::complex<float> >*, std::complex<float> );
template void Gadgetron::normalize< std::complex<float> >( cuNDArray< std::complex<float> >*, float );
template void Gadgetron::shrink1< std::complex<float> >( cuNDArray< std::complex<float> >*, float );
template void Gadgetron::shrinkd< std::complex<float> > ( cuNDArray< std::complex<float> >*, cuNDArray<float>*, float );

template boost::shared_ptr< cuNDArray<double> > Gadgetron::abs< std::complex<double> >( const cuNDArray< std::complex<double> >* );
template boost::shared_ptr< cuNDArray< std::complex<double> > > Gadgetron::sqrt< std::complex<double> >( const cuNDArray< std::complex<double> >* );
template boost::shared_ptr< cuNDArray<double> > Gadgetron::abs_square< std::complex<double> >( const cuNDArray< std::complex<double> >* );
template void Gadgetron::sqrt_inplace< std::complex<double> >( cuNDArray< std::complex<double> >* );
template boost::shared_ptr< cuNDArray< std::complex<double> > > Gadgetron::square< std::complex<double> >( const cuNDArray< std::complex<double> >* );
template void Gadgetron::square_inplace< std::complex<double> >( cuNDArray< std::complex<double> >* );
template boost::shared_ptr< cuNDArray< std::complex<double> > > Gadgetron::reciprocal< std::complex<double> >( const cuNDArray< std::complex<double> >* );
template void Gadgetron::reciprocal_inplace< std::complex<double> >( cuNDArray< std::complex<double> >* );
template boost::shared_ptr< cuNDArray< std::complex<double> > > Gadgetron::reciprocal_sqrt< std::complex<double> >( const cuNDArray< std::complex<double> >* );
template void Gadgetron::reciprocal_sqrt_inplace< std::complex<double> >( cuNDArray< std::complex<double> >* );
template void Gadgetron::clear< std::complex<double> >( cuNDArray< std::complex<double> >* );
template void Gadgetron::fill< std::complex<double> >( cuNDArray< std::complex<double> >*, std::complex<double> );
template void Gadgetron::normalize< std::complex<double> >( cuNDArray< std::complex<double> >*, double );
template void Gadgetron::shrink1< std::complex<double> >( cuNDArray< std::complex<double> >*, double );
template void Gadgetron::shrinkd< std::complex<double> > ( cuNDArray< std::complex<double> >*, cuNDArray<double>*, double );
*/
template boost::shared_ptr< cuNDArray<float> > Gadgetron::abs< complext<float> >( const cuNDArray< complext<float> >* );
template void Gadgetron::abs_inplace<complext<float> >(cuNDArray<complext<float> >*);
template boost::shared_ptr< cuNDArray< complext<float> > > Gadgetron::sqrt< complext<float> >( const cuNDArray< complext<float> >* );
template boost::shared_ptr< cuNDArray<float> > Gadgetron::abs_square< complext<float> >( const cuNDArray< complext<float> >* );
template void Gadgetron::sqrt_inplace< complext<float> >( cuNDArray< complext<float> >* );
template boost::shared_ptr< cuNDArray< complext<float> > > Gadgetron::square< complext<float> >( const cuNDArray< complext<float> >* );
template void Gadgetron::square_inplace< complext<float> >( cuNDArray< complext<float> >* );
template boost::shared_ptr< cuNDArray< complext<float> > > Gadgetron::reciprocal< complext<float> >( const cuNDArray< complext<float> >* );
template void Gadgetron::reciprocal_inplace< complext<float> >( cuNDArray< complext<float> >* );
template boost::shared_ptr< cuNDArray< complext<float> > > Gadgetron::reciprocal_sqrt< complext<float> >( const cuNDArray< complext<float> >* );
template void Gadgetron::reciprocal_sqrt_inplace< complext<float> >( cuNDArray< complext<float> >* );
template boost::shared_ptr< cuNDArray<complext<float> > > Gadgetron::sgn<complext<float> >( const cuNDArray<complext<float> >* );
template void Gadgetron::sgn_inplace<complext<float> >( cuNDArray<complext<float> >* );
template void Gadgetron::clear< complext<float> >( cuNDArray< complext<float> >* );
template void Gadgetron::fill< complext<float> >( cuNDArray< complext<float> >*, complext<float> );
template void Gadgetron::clamp< complext<float> >( cuNDArray< complext<float> >*, float, float );
template void Gadgetron::clamp_min< complext<float> >( cuNDArray< complext<float> >*, float );
template void Gadgetron::clamp_max< complext< float> >( cuNDArray<complext<float> >*, float );
template void Gadgetron::normalize< complext<float> >( cuNDArray< complext<float> >*, float );
template void Gadgetron::shrink1< complext<float> >( cuNDArray< complext<float> >*, float, cuNDArray< complext<float> >* );
template void Gadgetron::pshrink< complext<float> >( cuNDArray< complext<float> >*, float,float, cuNDArray< complext<float> >* );
template void Gadgetron::shrinkd< complext<float> > ( cuNDArray< complext<float> >*, cuNDArray<float>*, float, cuNDArray< complext<float> >* );
template void Gadgetron::pshrinkd< complext<float> > ( cuNDArray< complext<float> >*, cuNDArray<float>*, float,float, cuNDArray< complext<float> >* );

template boost::shared_ptr< cuNDArray<double> > Gadgetron::abs< complext<double> >( const cuNDArray< complext<double> >* );
template boost::shared_ptr< cuNDArray< complext<double> > > Gadgetron::sqrt< complext<double> >( const cuNDArray< complext<double> >* );
template boost::shared_ptr< cuNDArray<double> > Gadgetron::abs_square< complext<double> >( const cuNDArray< complext<double> >* );
template void Gadgetron::sqrt_inplace< complext<double> >( cuNDArray< complext<double> >* );
template boost::shared_ptr< cuNDArray< complext<double> > > Gadgetron::square< complext<double> >( const cuNDArray< complext<double> >* );
template void Gadgetron::square_inplace< complext<double> >( cuNDArray< complext<double> >* );
template boost::shared_ptr< cuNDArray< complext<double> > > Gadgetron::reciprocal< complext<double> >( const cuNDArray< complext<double> >* );
template void Gadgetron::reciprocal_inplace< complext<double> >( cuNDArray< complext<double> >* );
template boost::shared_ptr< cuNDArray< complext<double> > > Gadgetron::reciprocal_sqrt< complext<double> >( const cuNDArray< complext<double> >* );
template void Gadgetron::reciprocal_sqrt_inplace< complext<double> >( cuNDArray< complext<double> >* );
template boost::shared_ptr< cuNDArray<complext<double> > > Gadgetron::sgn<complext<double> >( const cuNDArray<complext<double> >* );
template void Gadgetron::sgn_inplace<complext<double> >( cuNDArray<complext<double> >* );
template void Gadgetron::clear< complext<double> >( cuNDArray< complext<double> >* );
template void Gadgetron::fill< complext<double> >( cuNDArray< complext<double> >*, complext<double> );
template void Gadgetron::clamp< complext<double> >( cuNDArray< complext<double> >*, double, double );
template void Gadgetron::clamp_min< complext<double> >( cuNDArray< complext<double> >*, double );
template void Gadgetron::clamp_max< complext<double> >( cuNDArray<complext<double> >*, double );
template void Gadgetron::normalize< complext<double> >( cuNDArray< complext<double> >*, double );
template void Gadgetron::shrink1< complext<double> >( cuNDArray< complext<double> >*, double, cuNDArray< complext<double> >* );
template void Gadgetron::pshrink< complext<double> >( cuNDArray< complext<double> >*, double, double, cuNDArray< complext<double> >* );
template void Gadgetron::shrinkd< complext<double> > ( cuNDArray< complext<double> >*, cuNDArray<double>*, double, cuNDArray< complext<double> >* );
template void Gadgetron::pshrinkd< complext<double> > ( cuNDArray< complext<double> >*, cuNDArray<double>*, double,double, cuNDArray< complext<double> >* );

template boost::shared_ptr< cuNDArray<float> > Gadgetron::real<float>( const cuNDArray<float>* );
template boost::shared_ptr< cuNDArray<float> > Gadgetron::imag<float>( const cuNDArray<float>* );
template boost::shared_ptr< cuNDArray<float> > Gadgetron::conj<float>( const cuNDArray<float>* );
template boost::shared_ptr< cuNDArray<float> > Gadgetron::real<float_complext>( const cuNDArray<float_complext>* );
template boost::shared_ptr< cuNDArray<float> > Gadgetron::imag<float_complext>( const cuNDArray<float_complext>* );
template boost::shared_ptr< cuNDArray<float_complext> > Gadgetron::conj<float_complext>( const cuNDArray<float_complext>* );
template boost::shared_ptr< cuNDArray<float_complext> > Gadgetron::real_to_complex<float_complext>( const cuNDArray<float>* );

template boost::shared_ptr< cuNDArray<double> > Gadgetron::real<double>( const cuNDArray<double>* );
template boost::shared_ptr< cuNDArray<double> > Gadgetron::imag<double>( const cuNDArray<double>* );
template boost::shared_ptr< cuNDArray<double> > Gadgetron::conj<double>( const cuNDArray<double>* );
template boost::shared_ptr< cuNDArray<double> > Gadgetron::real<double_complext>( const cuNDArray<double_complext>* );
template boost::shared_ptr< cuNDArray<double> > Gadgetron::imag<double_complext>( const cuNDArray<double_complext>* );
template boost::shared_ptr< cuNDArray<double_complext> > Gadgetron::conj<double_complext>( const cuNDArray<double_complext>* );
template boost::shared_ptr< cuNDArray<double_complext> > Gadgetron::real_to_complex<double_complext>( const cuNDArray<double>* );

template boost::shared_ptr< cuNDArray<double> > Gadgetron::convert_to<float,double>( const cuNDArray<float>* );
template boost::shared_ptr< cuNDArray<float> > Gadgetron::convert_to<double,float>( const cuNDArray<double>* );
template boost::shared_ptr< cuNDArray<double_complext> > Gadgetron::convert_to<float_complext,double_complext>( const cuNDArray<float_complext>* );
template boost::shared_ptr< cuNDArray<float_complext> > Gadgetron::convert_to<double_complext,float_complext>( const cuNDArray<double_complext>* );

template void Gadgetron::convert_to<float,double>( cuNDArray<float>*,cuNDArray<double>* );
template void Gadgetron::convert_to<double,float>( cuNDArray<double>*, cuNDArray<float>* );
template void Gadgetron::convert_to<float_complext,double_complext>( cuNDArray<float_complext>*,cuNDArray<double_complext>*  );
template void Gadgetron::convert_to<double_complext,float_complext>( cuNDArray<double_complext>*, cuNDArray<float_complext>*);
