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
#include "cuNDArray_utils.h"
#include "cuNDArray_blas.h"
#include "opticalFlowSolver.h"

namespace Gadgetron{

  template<class T, unsigned int D> class cuOpticalFlowSolver
    : public opticalFlowSolver< cuNDArray<T>,D >
  {
  public:

    cuOpticalFlowSolver() : opticalFlowSolver< cuNDArray<T>,D >() {}
    virtual ~cuOpticalFlowSolver() {}

  protected:

    // General tool to set up the block/grid dimensions
    //

    void setup_grid( dim3 *blockDim, dim3* gridDim, unsigned int number_of_elements,
                     unsigned int num_batches = 1, bool use_2d_blocks = false, unsigned int num_unknowns = D);

    // GPU-based computation of the spatial and temporal image gradient
    //

    virtual void core_grad_spatial( T *fixed_image, T *moving_image, T *gradient_image,
                                    typename uint64d<D>::Type matrix_size_moving,
                                    size_t number_of_batches_fixed,
                                    size_t number_of_batches_moving );

    virtual void core_grad_temporal( T *fixed_image, T *moving_image, T *gradient_image,
                                     typename uint64d<D>::Type matrix_size_moving,
                                     size_t number_of_batches_fixed,
                                     size_t number_of_batches_moving );
  };
}
