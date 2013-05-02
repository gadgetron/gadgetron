/** \file multiresRegistrationSolver.h
    Abstract class defining a multiresolution registration solver.
    Pure virtual functions are expected to do the actual work.
*/

#pragma once

#include "registrationSolver.h"

namespace Gadgetron{

  template<class ARRAY_TYPE_REAL, unsigned int D> class multiresRegistrationSolver 
    : public registrationSolver<ARRAY_TYPE_REAL>
  {
  protected:    
      typedef typename ARRAY_TYPE_REAL::element_type REAL;

  public:

    multiresRegistrationSolver() : registrationSolver<ARRAY_TYPE_REAL>(){
      num_multires_levels_ = 0;
      max_num_iterations_per_level_ = 500;
    }
    
    virtual ~multiresRegistrationSolver() {}
  
    // Utilities to specify the registration settings
    //
  
    virtual void set_num_multires_levels( unsigned int levels ) { 
      num_multires_levels_ = levels; }

    virtual void set_max_num_iterations_per_level( unsigned int iterations ) { 
      max_num_iterations_per_level_ = iterations; }
  
    // 
    // The main solver interface
    //

    virtual boost::shared_ptr<ARRAY_TYPE_REAL> 
    solve( registrationData<ARRAY_TYPE_REAL> *rd )
    {
      return registrationSolver<ARRAY_TYPE_REAL>::solve(rd);
    }

    virtual boost::shared_ptr<ARRAY_TYPE_REAL> solve( ARRAY_TYPE_REAL *fixed_image, 
						      ARRAY_TYPE_REAL *moving_image, 
						      bool input_normalization_allowed = false  )
    {
      // Some initial validity tests
      //
    
      if( !fixed_image || !moving_image ){
	BOOST_THROW_EXCEPTION( runtime_error("multiresRegistrationSolver::solve : invalid input pointer."));
      }

      if( !this->interpolator_.get() ){
	BOOST_THROW_EXCEPTION( runtime_error("multiresRegistrationSolver::solve : interpolator not set."));
      }
      
      // Normalize the input
      //
    
      ARRAY_TYPE_REAL *normalized_fixed;
      ARRAY_TYPE_REAL *normalized_moving;
    
      boost::shared_ptr<ARRAY_TYPE_REAL> 
	garbage_collector_fixed, garbage_collector_moving; 
    
      if( input_normalization_allowed ){
	normalized_fixed = fixed_image;
	normalized_moving = moving_image;
      }
      else{
	normalized_fixed = new ARRAY_TYPE_REAL(*fixed_image);
	normalized_moving = new ARRAY_TYPE_REAL(*moving_image);
	garbage_collector_fixed = boost::shared_ptr<ARRAY_TYPE_REAL>(normalized_fixed);
	garbage_collector_moving = boost::shared_ptr<ARRAY_TYPE_REAL>(normalized_moving);
      }
      
      normalize(normalized_fixed, REAL(1));
      normalize(normalized_moving, REAL(1));
    
      // Invoke multi-resolution solver
      // 
    
      if( this->output_mode_ >= registrationSolver<ARRAY_TYPE_REAL>::OUTPUT_VERBOSE ) {
	std::cout << std::endl << "Starting multiresolution registration " <<  std::endl;
      }
    
      boost::shared_ptr<ARRAY_TYPE_REAL> result = 
	solveMultiRes( num_multires_levels_, normalized_fixed, normalized_moving, this->stencil_.get() );
    
      return result;
    }
  
  protected:

    // Pure virtual fuctions to be implemented in a subclass
    //

    virtual void compute( ARRAY_TYPE_REAL *fixed_image, ARRAY_TYPE_REAL *moving_image, 
			  ARRAY_TYPE_REAL *stencil_image, boost::shared_ptr<ARRAY_TYPE_REAL> &result ) = 0;
  
    // The recursive multi-resolution solver
    //

    virtual boost::shared_ptr<ARRAY_TYPE_REAL> 
      solveMultiRes( unsigned int res_level, 
		     ARRAY_TYPE_REAL *fixed_image, 
		     ARRAY_TYPE_REAL *moving_image, 
		     ARRAY_TYPE_REAL *stencil_image )
    {
    
      boost::shared_ptr<ARRAY_TYPE_REAL> result;
      
      if (res_level>0){

	// 
	// We are not yet at the end of the multi-resolution chain
	//
      
	// Downsample input images (and stencil if provided)
	//

	boost::shared_ptr<ARRAY_TYPE_REAL> fixed_image_lowres = downsample<REAL,D>(fixed_image);
	boost::shared_ptr<ARRAY_TYPE_REAL> moving_image_lowres = downsample<REAL,D>(moving_image);
	boost::shared_ptr<ARRAY_TYPE_REAL> stencil_image_lowres = 
	  ((stencil_image) ? downsample<REAL,D>(stencil_image) : boost::shared_ptr<ARRAY_TYPE_REAL>());
      
	// Compute displacement field at the downsampled resolution
	//

	boost::shared_ptr<ARRAY_TYPE_REAL> result_lowres = 
	  solveMultiRes( res_level-1, fixed_image_lowres.get(), moving_image_lowres.get(), stencil_image_lowres.get() );
      
	// Clean up low resolution image data
	//
      
	fixed_image_lowres.reset();
	moving_image_lowres.reset();
	stencil_image_lowres.reset();

	// Upsample lowres results to current resolution
	//
      
	result = upsample<REAL,D>(result_lowres.get());
	*result *= REAL(2); // To adjust the flow vectors to the fact that the resolution is now twice as high

	// Clean up low resolution result
	//

	result_lowres.reset();

	// Some output to track our progress at runtime
	//

	if( this->output_mode_ >= registrationSolver<ARRAY_TYPE_REAL>::OUTPUT_VERBOSE ) {
	  std::cout << std::endl << "Multiresolution level " << res_level;
	}
      
	// Use estimated (lowres) motion to compute displacements at the current resolution
	//
	
	compute( fixed_image, moving_image, stencil_image, result );      
      }
      else{
	
	// 
	// We are now at the end of the multi-resolution chain
	//
      
	// Some output to track our progress at runtime
	//
      
	if( this->output_mode_ >= registrationSolver<ARRAY_TYPE_REAL>::OUTPUT_VERBOSE ) {
	  std::cout << std::endl << "Multiresolution level " << res_level << " (lowest)";
	}
      
	// Compute displacements at the current resolution (no estimate can be provided)
	//
      
	compute( fixed_image, moving_image, stencil_image, result );
      }
      
      return result;
    }
    
  protected:
    unsigned int num_multires_levels_;
    unsigned int max_num_iterations_per_level_;
  };
}
