#pragma once

#include "lwSolver.h"
#include "cuNDArray.h"
#include "cuNDArray_blas.h"

#include <iostream>

namespace Gadgetron{

  template <class T> class cuLwSolver
    : public lwSolver<cuNDArray<T> >
  {
  public:
  
    cuLwSolver() : lwSolver< cuNDArray<T> >() { set_device(-1); }
    virtual ~cuLwSolver() {}
  
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
}
