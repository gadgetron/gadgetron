#pragma once

#include "cgSolver.h"
#include "cuNDArray.h"
#include "cuCGPreconditioner.h"
#include "real_utilities.h"
#include "vector_td_utilities.h"
#include "ndarray_vector_td_utilities.h"

#include <iostream>


template <class REAL, class T> class cuCGSolver : public cgSolver< REAL, T, cuNDArray<T> >
{
public:

  //  cuCGSolver() : cgSolver< REAL, T, cuNDArray<T> >() { set_device(); }  
  cuCGSolver( int device=-1 ) : cgSolver< REAL, T, cuNDArray<T> >() { set_device(device); }

  virtual ~cuCGSolver() {}

  virtual bool pre_solve(cuNDArray<T> **rhs)
  {
    // Query and set device no
    //
    if( cudaGetDevice( &old_device_ ) != cudaSuccess ){
      solver_error( "cuCGSolver::pre_solve: unable to get current device." );
      return false;
    }
    //
    if( device_ != old_device_ && cudaSetDevice(device_) != cudaSuccess) {
      solver_error( "cuCGSolver:solve: unable to set specified device" );
      return false;
    }
  
    // Transfer arrays to compute device if necessary
    if( device_ != (*rhs)->get_device() ){
      new_rhs = new cuNDArray<T>(*(*rhs));    
      if( !new_rhs ){
	solver_error( "cuCGSolver::pre_solve: failed to copy rhs to the specified compute device." );
	return false;
      }    
      *rhs = new_rhs;
    }
    else
      new_rhs = 0x0;
  
    return true;
  }

  virtual bool post_solve(cuNDArray<T>**)
  {
    // Clear temporary rhs if allocated in 'pre_solve'
    if( new_rhs )
      delete new_rhs;
  
    if( device_ != old_device_ && cudaSetDevice(old_device_) != cudaSuccess) {
      solver_error( "cuCGSolver::solve: unable to restore device no" );
      return false;
    }
    return true;
  }

  virtual void solver_error( std::string err )
  {
    cgSolver< REAL, T, cuNDArray<T> >::solver_error(err);
    cudaSetDevice(old_device_);
  }

  virtual T solver_dot( cuNDArray<T> *x, cuNDArray<T> *y )
  {
    return cuNDA_dot<T>(x,y);
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
  cuNDArray<T> *new_rhs;
  
};
