#pragma once

#include "partialDerivativeOperator.h"
#include "cuPartialDerivativeOperator.h"
#include "hoCuNDArray.h"
#include "cudaDeviceManager.h"
#include "cuNDArray_operators.h"
#include "cuNDArray_elemwise.h"

#include "hoPartialDerivativeOperator.h"

namespace Gadgetron{

  template <class T, unsigned int D> class hoCuPartialDerivativeOperator :
    public linearOperator<hoCuNDArray<T> >
  {
  public: 
  
    hoCuPartialDerivativeOperator() : 
      linearOperator<hoCuNDArray<T> >(),dev(),hoDev(),_dimension(0) {}
  
    hoCuPartialDerivativeOperator( unsigned int dimension ) : 
      linearOperator<hoCuNDArray<T> >(),dev(dimension),hoDev(dimension), _dimension(dimension){ }

    virtual ~hoCuPartialDerivativeOperator() {}

    //TODO: Generalize to work if we can fit just the 1 single dimension on the gpu
    virtual void mult_M(hoCuNDArray<T>* in, hoCuNDArray<T>* out, bool accumulate)
    {
      size_t free = cudaDeviceManager::Instance()->getFreeMemory();

      if( free/sizeof(T) < in->get_number_of_elements()*2)
	throw std::runtime_error("hoCuPartialDerivativeOperator: not enough device memory");
      cuNDArray<T> cuIn(in);
      cuNDArray<T> cuOut(out->get_dimensions());

      if (accumulate) cuOut =cuNDArray<T>(out);

      dev.mult_M(&cuIn,&cuOut,accumulate);

      cudaMemcpy(out->get_data_ptr(),cuOut.get_data_ptr(),out->get_number_of_elements()*sizeof(T),cudaMemcpyDeviceToHost);
    	//hoDev.mult_M(in,out,accumulate);
    }

    //TODO: Generalize to work if we can fit just the 1 single dimension on the gpu
    virtual void mult_MH(hoCuNDArray<T>* in, hoCuNDArray<T>* out, bool accumulate)
    {

      size_t free = cudaDeviceManager::Instance()->getFreeMemory();

      if( free/sizeof(T) < in->get_number_of_elements()*2)
	throw std::runtime_error("hoCuPartialDerivativeOperator: not enough device memory");
      cuNDArray<T> cuIn(in);
      cuNDArray<T> cuOut(out->get_dimensions());

      if (accumulate) cuOut =cuNDArray<T>(out);

      dev.mult_MH(&cuIn,&cuOut,accumulate);

      cudaMemcpy(out->get_data_ptr(),cuOut.get_data_ptr(),out->get_number_of_elements()*sizeof(T),cudaMemcpyDeviceToHost);

    	//hoDev.mult_MH(in,out,accumulate);
    }

    //TODO: Generalize to work if we can fit just the 1 single dimension on the gpu
    virtual void mult_MH_M(hoCuNDArray<T>* in, hoCuNDArray<T>* out, bool accumulate)
    {

      size_t free = cudaDeviceManager::Instance()->getFreeMemory();

      if( free/sizeof(T) < in->get_number_of_elements()*2)
	throw std::runtime_error("hoCuPartialDerivativeOperator: not enough device memory");
      cuNDArray<T> cuIn(in);
      cuNDArray<T> cuOut(out->get_dimensions());

      if (accumulate) cuOut =cuNDArray<T>(out);

      dev.mult_MH_M(&cuIn,&cuOut,accumulate);

      cudaMemcpy(out->get_data_ptr(),cuOut.get_data_ptr(),out->get_number_of_elements()*sizeof(T),cudaMemcpyDeviceToHost);

    	//hoDev.mult_MH_M(in,out,accumulate);
    }

  protected:
    cuPartialDerivativeOperator<T,D> dev;
    hoPartialDerivativeOperator<T,D> hoDev;
    unsigned int _dimension;
  };
}
