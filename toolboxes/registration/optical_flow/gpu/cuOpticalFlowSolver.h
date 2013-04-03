#pragma once

#include "multiresRegistrationSolver.h"
#include "cuLinearResampleOperator.h"
#include "cuNDArray.h"

template<class REAL, unsigned int D> class EXPORTGPUREG cuOpticalFlowSolver 
  : public multiresRegistrationSolver< REAL, cuNDArray<REAL> >
{
  
 public:
  
  cuOpticalFlowSolver() : multiresRegistrationSolver< REAL, cuNDArray<REAL> >(){ 
    limit_ = (REAL) 0.01;
    set_interpolator( boost::shared_ptr< cuLinearResampleOperator<REAL,REAL,D> >
		      (new cuLinearResampleOperator<REAL,REAL,D>) );
  } 
  
  virtual ~cuOpticalFlowSolver() {}
  
  // Set termination threshold
  inline void set_limit( REAL limit ) { limit_ = limit; }
  
 protected:
  
  // The actual solver
  //
  
  virtual boost::shared_ptr< cuNDArray<REAL> > 
    core_solver( cuNDArray<REAL> *gradient_image, cuNDArray<REAL> *stencil_image ) = 0;  
  
  // Shared functionality
  //
  
  virtual bool normalize( cuNDArray<REAL> *image );
  virtual bool compute( cuNDArray<REAL> *fixed_image, cuNDArray<REAL> *moving_image, cuNDArray<REAL> *stencil_image, 
			boost::shared_ptr< cuNDArray<REAL> > &result_in_out );  
  virtual boost::shared_ptr< cuNDArray<REAL> > downsample( cuNDArray<REAL> *image );
  virtual boost::shared_ptr< cuNDArray<REAL> > upsample( cuNDArray<REAL> *displacements );
  virtual boost::shared_ptr< cuNDArray<REAL> > grad( cuNDArray<REAL> *fixed_image, cuNDArray<REAL> *moving_image );  
  virtual bool setup_grid( dim3 *blockDim, dim3* gridDim, unsigned int number_of_elements, 
			   unsigned int num_batches = 1, bool used_2d_blocks = false, unsigned int num_unknowns = D );
  
 protected:
  REAL limit_;
};
