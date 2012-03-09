#pragma once

#include "convolutionOperator.h"
#include "cuMatrixOperator_macros.h"
#include "cuNDArray.h"
#include "cuNDFFT.h"
#include "vector_td_utilities.h"
#include "ndarray_vector_td_utilities.h"


template <class REAL, unsigned int D> class cuConvolutionOperator 
	: public convolutionOperator<REAL, cuNDArray<complext<REAL> >, D >
{
  
public:
  
  cuConvolutionOperator( int device = -1 ) : convolutionOperator<REAL, cuNDArray<complext<REAL> >, D>() { set_device(device); }
  virtual ~cuConvolutionOperator() {}
  
  virtual int mult_M( cuNDArray<complext<REAL> > *in, cuNDArray<complext<REAL> > *out, bool accumulate = false )
  {
    _set_device();
    int res = convolutionOperator<REAL, cuNDArray<complext<REAL> >, D >::mult_M( in, out, accumulate );
    _restore_device();

    return res;
  }

  virtual int mult_MH( cuNDArray<complext<REAL> > *in, cuNDArray<complext<REAL> > *out, bool accumulate = false )
  {
    _set_device();
    int res = convolutionOperator<REAL, cuNDArray<complext<REAL> >,D>::mult_MH( in, out, accumulate );
    _restore_device();
    return res;
  }

  virtual int mult_MH_M( cuNDArray<complext<REAL> > *in, cuNDArray<complext<REAL> > *out, bool accumulate = false )
  {
    _set_device();
    int res = convolutionOperator<REAL, cuNDArray<complext<REAL> >, D>::mult_MH_M( in, out, accumulate );
    _restore_device();
    return res;
  }

  virtual bool operator_fft( bool forwards_transform, cuNDArray<complext<REAL> > *image )
  {
    int res;

    if( forwards_transform )
      res = cuNDFFT<complext<REAL> >().fft(image);
    else
      res = cuNDFFT<complext<REAL> >().ifft(image);
    
    if( res<0 )
      return false;
    else
      return true;    
  }
  
  virtual bool operator_scale( cuNDArray<complext<REAL> > *kernel, cuNDArray<complext<REAL> > *image, bool conjugate_kernel )
  {
    if( conjugate_kernel )
      return cuNDA_scale_conj( kernel, image );
    else
      return cuNDA_scale( kernel, image );
  }  
  
  virtual bool operator_xpy( cuNDArray<complext<REAL> > *x, cuNDArray<complext<REAL> > *y )
  {
    return cuNDA_axpy<complext<REAL> >( complext<REAL>(1), x, y );
  }

  virtual bool operator_mirror( cuNDArray<complext<REAL> > *in, cuNDArray<complext<REAL> > *out )
  {
    return cuNDA_origin_mirror<complext<REAL>,D>( in, out, false );
  }

  virtual bool operator_expand( cuNDArray<complext<REAL> > *in, cuNDArray<complext<REAL> > *out )
  {
    return cuNDA_expand_with_zero_fill<complext<REAL>,D>( in, out );
  }
  
  virtual bool operator_crop( cuNDArray<complext<REAL> > *in, cuNDArray<complext<REAL> > *out )
  {
    typename uintd<D>::Type offset = vector_to_uintd<D>(*(in->get_dimensions().get()))>>2;
    return cuNDA_crop<complext<REAL>,D>( offset, in, out );
  }
  
  DECLARE_MATRIX_OPERATOR_DEVICE_SUPPORT(cuConvolutionOperator)
};
