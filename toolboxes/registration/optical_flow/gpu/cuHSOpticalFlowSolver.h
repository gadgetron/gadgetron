#pragma once

#include "cuOpticalFlowSolver.h"

template<class REAL, unsigned int D> class EXPORTGPUREG cuHSOpticalFlowSolver 
  : public cuOpticalFlowSolver<REAL, D>
{
  
 public:

  // Constructors / destructors
  //
  
  cuHSOpticalFlowSolver() : cuOpticalFlowSolver<REAL,D>(){ 
    alpha_ = REAL(0.1); 
  } 
  
  virtual ~cuHSOpticalFlowSolver() {}
  
  // Set the regularization weight
  //
  
  inline void set_alpha( REAL alpha ) { alpha_ = alpha; }
  
 protected:  
  virtual boost::shared_ptr< cuNDArray<REAL> > 
    core_solver( cuNDArray<REAL> *gradient_image, cuNDArray<REAL> *stencil_image );  
  
 protected:
  REAL alpha_;
};
