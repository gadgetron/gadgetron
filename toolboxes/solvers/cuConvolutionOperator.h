#pragma once

#include "convolutionOperator.h"
#include "cuMatrixOperator_macros.h"
#include "vector_td.h"
#include "cuNDArray.h"
#include "cuNDFFT.h"
#include "ndarray_vector_td_utilities.h"

template <class REAL> class cuConvolutionOperator : public convolutionOperator<REAL, cuNDArray<typename complext<REAL>::Type> >
{
  
public:
  
  cuConvolutionOperator( int device = -1 ) : convolutionOperator<REAL, cuNDArray<typename complext<REAL>::Type> >() { set_device(device); }
  virtual ~cuConvolutionOperator() {}
  
  virtual int mult_M( cuNDArray<typename complext<REAL>::Type> *in, cuNDArray<typename complext<REAL>::Type> *out, bool accumulate = false )
  {
    set_device();
    int res = convolutionOperator<REAL, cuNDArray<typename complext<REAL>::Type> >::mult_M( in, out, accumulate );
    restore_device();
    return res;
  }

  virtual int mult_MH( cuNDArray<typename complext<REAL>::Type> *in, cuNDArray<typename complext<REAL>::Type> *out, bool accumulate = false )
  {
    set_device();
    int res = convolutionOperator<REAL, cuNDArray<typename complext<REAL>::Type> >::mult_MH( in, out, accumulate );
    restore_device();
    return res;
  }

  virtual int mult_MH_M( cuNDArray<typename complext<REAL>::Type> *in, cuNDArray<typename complext<REAL>::Type> *out, bool accumulate = false )
  {
    set_device();
    int res = convolutionOperator<REAL, cuNDArray<typename complext<REAL>::Type> >::mult_MH_M( in, out, accumulate );
    restore_device();
    return res;
  }

  virtual bool operator_fft( bool forwards_transform, cuNDArray<typename complext<REAL>::Type> *image )
  {
    int res;

    if( forwards_transform )
      res = cuNDFFT<typename complext<REAL>::Type>().fft(image);
    else
      res = cuNDFFT<typename complext<REAL>::Type>().ifft(image);
    
    if( res<0 )
      return false;
    else
      return true;    
  }
  
  virtual bool operator_scale( bool conjugate_kernel, cuNDArray<typename complext<REAL>::Type> *kernel, cuNDArray<typename complext<REAL>::Type> *image )
  {
    if( conjugate_kernel )
      return cuNDA_scale_conj( kernel, image );
    else
      return cuNDA_scale( kernel, image );
  }  
  
  virtual bool operator_xpy( cuNDArray<typename complext<REAL>::Type> *x, cuNDArray<typename complext<REAL>::Type> *y )
  {
    return cuNDA_axpy<typename complext<REAL>::Type>( get_one<typename complext<REAL>::Type>(), x, y );
  }

  DECLARE_MATRIX_OPERATOR_DEVICE_SUPPORT(cuConvolutionOperator)
};
