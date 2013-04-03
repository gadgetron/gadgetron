#pragma once

#include "cuOpticalFlowSolver.h"

template<class REAL, unsigned int D> class EXPORTGPUREG cuCKOpticalFlowSolver 
  : public cuOpticalFlowSolver<REAL, D>
{
  
 public:

  // Constructors / destructors
  //
  
  cuCKOpticalFlowSolver() : cuOpticalFlowSolver<REAL,D>(){ 
    alpha_ = REAL(0.1); 
    beta_ = REAL(0.01); 
  } 
  
  virtual ~cuCKOpticalFlowSolver() {}
  
  // Set the regularization weight
  //
  
  inline void set_alpha( REAL alpha ) { alpha_ = alpha; }
  inline void set_beta( REAL beta ) { beta_ = beta; }
  
 protected:  
  virtual boost::shared_ptr< cuNDArray<REAL> > 
    core_solver( cuNDArray<REAL> *gradient_image, cuNDArray<REAL> *stencil_image );  
  
 protected:
  REAL alpha_;
  REAL beta_;
};
