/** \file opticalFlowSolver.h
    \brief Abstract class defining an optical flow registration solver.

    Pure virtual functions are expected to do the actual work 
    - on the CPU and GPU respectively.
*/

#pragma once

#include "multiresRegistrationSolver.h"
#include "resampleOperator.h"
#include "vector_td_utilities.h"

#include <algorithm>

namespace Gadgetron{

  template<class ARRAY_TYPE_REAL, unsigned int D> class opticalFlowSolver 
    : public multiresRegistrationSolver<ARRAY_TYPE_REAL,D>
  {  
  protected:
    typedef typename ARRAY_TYPE_REAL::element_type REAL;

  public:

    opticalFlowSolver() : multiresRegistrationSolver<ARRAY_TYPE_REAL,D>(){ 
      limit_ = REAL(0.01);
    } 

    virtual ~opticalFlowSolver() {}

    // Set termination threshold
    inline void set_limit( REAL limit ) { limit_ = limit; }

  protected:

    // Inherited from the multiresolution solver
    //

    virtual void compute( ARRAY_TYPE_REAL *fixed_image, ARRAY_TYPE_REAL *moving_image, ARRAY_TYPE_REAL *stencil_image, 
                          boost::shared_ptr<ARRAY_TYPE_REAL> &result_in_out )
    {
      // Test the validity of the input images
      //

      if( !fixed_image || !moving_image ){
        throw std::runtime_error("opticalFlowSolver::compute(): illegal input array received.");
      }

      if( prod(from_std_vector<size_t,D>(*fixed_image->get_dimensions().get())) != 
          prod(from_std_vector<size_t,D>(*moving_image->get_dimensions().get())) ){
        throw std::runtime_error("opticalFlowSolver::compute(): core image dimensions (excluding batches) mismatch.");
      }

      if( stencil_image && 
          prod(from_std_vector<size_t,D>(*fixed_image->get_dimensions().get())) != 
          prod(from_std_vector<size_t,D>(*stencil_image->get_dimensions().get())) ){
        throw std::runtime_error("opticalFlowSolver::compute(): stencil image dimensions mismatch fixed/moving image dimensions.");
      }

      if( result_in_out.get() && 
          !( result_in_out->get_number_of_dimensions() > D ||
             result_in_out->get_size(result_in_out->get_number_of_dimensions()-1) == D )){
        throw std::runtime_error("opticalFlowSolver::compute(): input displacements dimensionality mismatch");
      }


      boost::shared_ptr<ARRAY_TYPE_REAL> grad_image = grad( fixed_image, moving_image );


      // Invoke core solver (e.g. Horn-Schunk, Cornelius-Kanade, ...)
      //

      boost::shared_ptr<ARRAY_TYPE_REAL> displacements = core_solver( grad_image.get(), stencil_image );

      // If an input vector field was provided then our result should be added element-wise
      // 

      if( result_in_out.get() ){
        *result_in_out += *displacements;
      }
      else{    
        result_in_out = displacements;
      }
    }

    // Compute the gradient
    //

    virtual boost::shared_ptr<ARRAY_TYPE_REAL> grad( ARRAY_TYPE_REAL *fixed_image, ARRAY_TYPE_REAL *moving_image )
    {
      // Sanity checks
      //

      if( !fixed_image || !moving_image ){
        throw std::runtime_error("opticalFlowSolver::grad(): illegal input received.");
      }

      if( !((moving_image->get_number_of_elements() % fixed_image->get_number_of_elements()) == 0 ||
            (fixed_image->get_number_of_elements() % moving_image->get_number_of_elements()) == 0 )){
        throw std::runtime_error("opticalFlowSolver::grad(): fixed/moving image dimensions mismatch.");
      }

      // Determine dimension size of the gradient field:
      // D spatial dimensions plus one temporal dimension
      //

      std::vector<size_t> grad_dims;

      (fixed_image->get_number_of_elements()<moving_image->get_number_of_elements() )
        ? grad_dims = *moving_image->get_dimensions() : grad_dims = *fixed_image->get_dimensions();

      grad_dims.push_back(D+1); 

      boost::shared_ptr<ARRAY_TYPE_REAL> grad_image(new ARRAY_TYPE_REAL(grad_dims));

      // Setup for the spatial partial derivatives
      //

      typename uint64d<D>::Type matrix_size_fixed = from_std_vector<size_t,D>( *fixed_image->get_dimensions() );
      typename uint64d<D>::Type matrix_size_moving = from_std_vector<size_t,D>( *moving_image->get_dimensions() );

      if( matrix_size_fixed != matrix_size_moving ){
        throw std::runtime_error("opticalFlowSolver::grad(): fixed/moving image dimensions mismatch (2).");
      }

      // Ignoring the batch dimensions the fixed and moving images have the same number of elements
      //

      size_t number_of_elements = prod(matrix_size_moving);
      size_t number_of_batches_fixed = 1;
      size_t number_of_batches_moving = 1;

      for( size_t d=D; d<fixed_image->get_number_of_dimensions(); d++ ){
        number_of_batches_fixed *= fixed_image->get_size(d);
      }

      for( size_t d=D; d<moving_image->get_number_of_dimensions(); d++ ){
        number_of_batches_moving *= moving_image->get_size(d);
      }

      // Compute spatial partial derivatives
      //

      core_grad_spatial( fixed_image->get_data_ptr(), moving_image->get_data_ptr(), grad_image->get_data_ptr(), 
                         matrix_size_moving, number_of_batches_fixed, number_of_batches_moving );

      // Compute temporal partial derivatives
      //

      core_grad_temporal( fixed_image->get_data_ptr(), moving_image->get_data_ptr(), 
                          grad_image->get_data_ptr()+number_of_elements*std::max(number_of_batches_moving, number_of_batches_fixed)*D, 
                          matrix_size_moving, number_of_batches_fixed, number_of_batches_moving );

      return grad_image;
    }

    // The actual work is being done in these functions, to be implemented on both host and device
    //

    virtual boost::shared_ptr<ARRAY_TYPE_REAL> core_solver( ARRAY_TYPE_REAL *gradient_image, ARRAY_TYPE_REAL *stencil_image ) = 0;      

    virtual void core_grad_spatial( REAL *fixed_image, REAL *moving_image, REAL *gradient_image, 
                                    typename uint64d<D>::Type matrix_size_moving, 
                                    size_t number_of_batches_fixed, 
                                    size_t number_of_batches_moving ) = 0;

    virtual void core_grad_temporal( REAL *fixed_image, REAL *moving_image, REAL *gradient_image, 
                                     typename uint64d<D>::Type matrix_size_moving, 
                                     size_t number_of_batches_fixed, 
                                     size_t number_of_batches_moving ) = 0;

  protected:
    REAL limit_;
  };
}
