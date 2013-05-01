/** \file cuOpticalFlowSolver.h
    \brief Abstract class for a GPU-based optical flow registration solver.

    cuOpticalFlowSolver is derived from class opticalFlowSolver 
    and implements the computation of the spatial and temporal gradients.
    A pure virtual function is expected to implement the specific algorithm (Horn-Schunck, Cornelius-Kanade).
*/

#pragma once

#include "cuNDArray.h"
#include "cuNDArray_operators.h"
#include "cuNDArray_elemwise.h"
#include "cuRegistration_utils.h"
#include "opticalFlowSolver.h"
#include "gpureg_export.h"

namespace Gadgetron{

  template<class REAL, unsigned int D> class EXPORTGPUREG cuOpticalFlowSolver 
    : public opticalFlowSolver< cuNDArray<REAL>,D >
  {  
  public:
  
    cuOpticalFlowSolver() : opticalFlowSolver< cuNDArray<REAL>,D >() {}   
    virtual ~cuOpticalFlowSolver() {}
    
  protected:

    // Inherited and still pure virtual...
    virtual boost::shared_ptr< cuNDArray<REAL> > core_solver( cuNDArray<REAL> *gradient_image, cuNDArray<REAL> *stencil_image ) = 0;      

    // General tool to set up the block/grid dimensions
    void setup_grid( dim3 *blockDim, dim3* gridDim, unsigned int number_of_elements, 
		     unsigned int num_batches = 1, bool use_2d_blocks = false, unsigned int num_unknowns = D);  
 
    // GPU-based computation of the spatial and temporal image gradient
    //
    
    virtual void core_grad_spatial( REAL *fixed_image, REAL *moving_image, REAL *gradient_image, 
				    typename uintd<D>::Type matrix_size_moving, 
				    unsigned int number_of_batches_fixed, 
				    unsigned int number_of_batches_moving );
    
    virtual void core_grad_temporal( REAL *fixed_image, REAL *moving_image, REAL *gradient_image, 
				     typename uintd<D>::Type matrix_size_moving, 
				     unsigned int number_of_batches_fixed, 
				     unsigned int number_of_batches_moving );
  };  
}
