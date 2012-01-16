#pragma once

#include "convolutionOperator.h"
#include "cuMatrixOperator_macros.h"
#include "cuNDArray.h"
#include "cuNDFFT.h"
#include "vector_td_utilities.h"
#include "ndarray_vector_td_utilities.h"

template <class REAL, unsigned int D> class cuConvolutionOperator : public convolutionOperator<REAL, cuNDArray<typename complext<REAL>::Type>, D >
{
  
public:
  
  cuConvolutionOperator( int device = -1 ) : convolutionOperator<REAL, cuNDArray<typename complext<REAL>::Type>, D>() { set_device(device); }
  virtual ~cuConvolutionOperator() {}
  
  virtual int mult_M( cuNDArray<typename complext<REAL>::Type> *in, cuNDArray<typename complext<REAL>::Type> *out, bool accumulate = false )
  {
    set_device();
    int res = convolutionOperator<REAL, cuNDArray<typename complext<REAL>::Type>, D >::mult_M( in, out, accumulate );
    restore_device();
    return res;
  }

  virtual int mult_MH( cuNDArray<typename complext<REAL>::Type> *in, cuNDArray<typename complext<REAL>::Type> *out, bool accumulate = false )
  {
    set_device();
    int res = convolutionOperator<REAL, cuNDArray<typename complext<REAL>::Type>,D>::mult_MH( in, out, accumulate );
    restore_device();
    return res;
  }

  virtual int mult_MH_M( cuNDArray<typename complext<REAL>::Type> *in, cuNDArray<typename complext<REAL>::Type> *out, bool accumulate = false )
  {
    set_device();
    int res = convolutionOperator<REAL, cuNDArray<typename complext<REAL>::Type>, D>::mult_MH_M( in, out, accumulate );
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
  
  virtual bool operator_scale( cuNDArray<typename complext<REAL>::Type> *kernel, cuNDArray<typename complext<REAL>::Type> *image, bool conjugate_kernel )
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

  virtual bool operator_mirror( cuNDArray<typename complext<REAL>::Type> *in, cuNDArray<typename complext<REAL>::Type> *out )
  {
    return cuNDA_origin_mirror<typename complext<REAL>::Type,D>( in, out, false );
  }

  virtual bool operator_expand( cuNDArray<typename complext<REAL>::Type> *in, cuNDArray<typename complext<REAL>::Type> *out )
  {
    return cuNDA_expand_with_zero_fill<typename complext<REAL>::Type,D>( in, out );
  }
  
  virtual bool operator_crop( cuNDArray<typename complext<REAL>::Type> *in, cuNDArray<typename complext<REAL>::Type> *out ) 
  {
    typename uintd<D>::Type offset = vector_to_uintd<D>(*(in->get_dimensions().get()))>>2;
    return cuNDA_crop<typename complext<REAL>::Type,D>( offset, in, out );
  }
  
  DECLARE_MATRIX_OPERATOR_DEVICE_SUPPORT(cuConvolutionOperator)
};
