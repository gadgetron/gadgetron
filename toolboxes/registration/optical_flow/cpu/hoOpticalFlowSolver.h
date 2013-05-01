/** \file hoOpticalFlowSolver.h
    \brief Abstract class for a CPU-based optical flow registration solver.

    hoOpticalFlowSolver is derived from class opticalFlowSolver 
    and implements the computation of the spatial and temporal gradients.
    A pure virtual function is expected to implement the specific algorithm (Horn-Schunck, Cornelius-Kanade).
*/

#pragma once

#include "hoNDArray.h"
#include "hoNDArray_operators.h"
#include "hoNDArray_elemwise.h"
#include "hoRegistration_utils.h"
#include "opticalFlowSolver.h"
#include "cpureg_export.h"

namespace Gadgetron{

  template<class REAL, unsigned int D> class EXPORTCPUREG hoOpticalFlowSolver 
    : public opticalFlowSolver< hoNDArray<REAL>,D >
  {  
  public:
  
    hoOpticalFlowSolver() : opticalFlowSolver< hoNDArray<REAL>,D >() {}   
    virtual ~hoOpticalFlowSolver() {}
    
  protected:

    // Inherited and still pure virtual...
    virtual boost::shared_ptr< hoNDArray<REAL> > core_solver( hoNDArray<REAL> *gradient_image, hoNDArray<REAL> *stencil_image ) = 0;      

    // CPU-based computation of the spatial and temporal image gradient
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
