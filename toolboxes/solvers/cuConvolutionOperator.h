#pragma once
#include "ndarray_vector_td_utilities.h"
#include "convolutionOperator.h"
#include "cuLinearOperator_macros.h"
#include "cuNDArray.h"
#include "cuNDFFT.h"
#include "vector_td_utilities.h"



template <class REAL, unsigned int D> class cuConvolutionOperator 
	: public convolutionOperator<cuNDArray<complext<REAL> >, D >
{
  
public:
  
  cuConvolutionOperator() : convolutionOperator<cuNDArray<complext<REAL> >, D>() {  }
  virtual ~cuConvolutionOperator() {}


  virtual void operator_fft( bool forwards_transform, cuNDArray<complext<REAL> > *image )
  {


    if( forwards_transform )
      cuNDFFT<complext<REAL> >().fft(image);
    else
      cuNDFFT<complext<REAL> >().ifft(image);

  }
  

    
  virtual boost::shared_ptr< linearOperator<cuNDArray< complext<REAL> > > > clone()
  {
    return linearOperator< cuNDArray< complext<REAL> > >::clone(this);
  }


};
