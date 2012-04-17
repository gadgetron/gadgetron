#pragma once

#include "lwSolver.h"
#include "cuNDArray.h"
#include "ndarray_vector_td_utilities.h"

#include <iostream>

template <class REAL, class T> class cuLwSolver 
  : public lwSolver< REAL, T, cuNDArray<T> >
{
public:
  
  cuLwSolver() : lwSolver< REAL, T, cuNDArray<T> >() { set_device(-1); }
  virtual ~cuLwSolver() {}
  
  virtual bool solver_clear( cuNDArray<T> *x )
  {
    return cuNDA_clear<T>(x);
  }

  virtual bool solver_axpy( T a, cuNDArray<T> *x, cuNDArray<T> *y )
  {
    return cuNDA_axpy<T>(a,x,y);
  }

  virtual bool set_device( int device )
  { 
    device_ = device;
    
    if( device<0 ){
      
      int old_device;  
      
      if( cudaGetDevice( &old_device ) != cudaSuccess ){
	std::cerr << "cuLwSolver::set_device: unable to get current device." << std::endl ;
	return false;
      }
      
      device_ = old_device;
    }
    
    return true;
  }
    
protected:
  int device_;
  int old_device_;
};
