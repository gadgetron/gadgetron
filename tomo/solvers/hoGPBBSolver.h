#pragma once

#include "gpBBSolver.h"
#include "hoGTBLAS.h"
namespace Gadgetron{

template<class T> class hoGPBBSolver: public gpBBSolver<hoNDArray<T> >{
public:
	typedef typename realType<T>::type REAL;
	hoGPBBSolver():gpBBSolver<hoNDArray<T> >(){

	}

protected:
  virtual void solver_non_negativity_filter(hoNDArray<T> *x,hoNDArray<T> *g)
  {
    T* x_ptr = x->get_data_ptr();
    T* g_ptr = g->get_data_ptr();
    for (int i =0; i < x->get_number_of_elements(); i++){
      if ( real(x_ptr[i]) < REAL(0) && real(g_ptr[i]) > 0) g_ptr[i]=T(0);
    }

  }

  virtual void solver_reciprocal_clamp( hoNDArray<T>* x_arr, REAL threshold){
  	T* x = x_arr->get_data_ptr();
  	for (int i = 0;  i  < x_arr->get_number_of_elements(); i++){
  		if (real(x[i]) < threshold) x[i]= T(0);
  				  else x[i]= T(1)/x[i];
  	}
  }
};
}
