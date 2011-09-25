#pragma once

#include "sbSolver.h"
#include "cuNDArray.h"
#include "real_utilities.h"
#include "vector_td_utilities.h"
#include "ndarray_vector_td_utilities.h"

template <class REAL, class T> class cuSBSolver : public sbSolver< REAL, T, cuNDArray<T> >
{
public:
  
  cuSBSolver( int device=-1 ) : sbSolver< REAL, T, cuNDArray<T> >() { set_device(device); }
  virtual ~cuSBSolver() {}

  virtual void solver_error( std::string err )
  {
    sbSolver< REAL, T, cuNDArray<T> >::solver_error(err);
    cudaSetDevice(old_device_);
  }

  virtual bool solver_clear( cuNDArray<T> *x )
  {
    return cuNDA_clear<T>(x);
  }

  virtual bool solver_scal( T a, cuNDArray<T> *x )
  {
    return cuNDA_scal<T>(a,x);
  }

  virtual bool solver_axpy( T a, cuNDArray<T> *x, cuNDArray<T> *y )
  {
    return cuNDA_axpy<T>(a,x,y);
  }

  virtual bool solver_shrink( REAL reciprocal_scale, cuNDArray<T> *in, cuNDArray<T> *out )
  {
    return cuNDA_shrink<REAL,T>( reciprocal_scale, in, out );    
  }

  virtual bool set_device( int device )
  { 
    device_ = device;
    
    if( device<0 ){
      
      int old_device;  
      
      if( cudaGetDevice( &old_device ) != cudaSuccess ){
	std::cerr << "cuCGSolver::set_device: unable to get current device." << std::endl ;
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
