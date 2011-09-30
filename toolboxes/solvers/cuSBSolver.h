#pragma once

#include "sbSolver.h"
#include "cuNDArray.h"
#include "real_utilities.h"
#include "vector_td_utilities.h"
#include "ndarray_vector_td_utilities.h"

template <class REAL, class T> class cuSBSolver : public sbSolver< REAL, T, cuNDArray<REAL>, cuNDArray<T> >
{
public:
  
  cuSBSolver( int device=-1 ) : sbSolver< REAL, T, cuNDArray<REAL>, cuNDArray<T> >() { set_device(device); }
  virtual ~cuSBSolver() {}

  virtual void solver_error( std::string err )
  {
    sbSolver< REAL, T, cuNDArray<REAL>, cuNDArray<T> >::solver_error(err);
    cudaSetDevice(old_device_);
  }

  virtual bool solver_clear( cuNDArray<T> *x )
  {
    return cuNDA_clear<T>(x);
  }

  virtual bool solver_clear( cuNDArray<REAL> *x )
  {
    return cuNDA_clear<REAL>(x);
  }

  virtual bool solver_scal( T a, cuNDArray<T> *x )
  {
    return cuNDA_scal<T>(a,x);
  }

  virtual bool solver_sqrt( cuNDArray<REAL> *x )
  {
    return cuNDA_sqrt<REAL>(x);
  }

  virtual bool solver_axpy( T a, cuNDArray<T> *x, cuNDArray<T> *y )
  {
    return cuNDA_axpy<T>(a,x,y);
  }

  virtual bool solver_axpy( REAL a, cuNDArray<REAL> *x, cuNDArray<REAL> *y )
  {
    return cuNDA_axpy<REAL>(a,x,y);
  }

  virtual REAL solver_asum( cuNDArray<T> *x )
  {
    return cuNDA_asum<REAL,T>(x);
  }

  virtual boost::shared_ptr< cuNDArray<REAL> > solver_norm( cuNDArray<T> *x )
  {
    return cuNDA_norm<REAL,T>(x);
  }

  virtual boost::shared_ptr< cuNDArray<REAL> > solver_norm_squared( cuNDArray<T> *x )
  {
    return cuNDA_norm_squared<REAL,T>(x);
  }
  
  virtual bool solver_shrink1( REAL reciprocal_scale, cuNDArray<T> *in, cuNDArray<T> *out )
  {
    return cuNDA_shrink1<REAL,T>( reciprocal_scale, in, out );    
  }

  virtual bool solver_shrinkd( REAL reciprocal_scale, cuNDArray<REAL> *s_k, cuNDArray<T> *in, cuNDArray<T> *out )
  {
    return cuNDA_shrinkd<REAL,T>( reciprocal_scale, s_k, in, out );    
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
