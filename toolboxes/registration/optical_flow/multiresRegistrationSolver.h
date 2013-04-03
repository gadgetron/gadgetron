#pragma once

#include "registrationSolver.h"

template<class REAL, class ARRAY_TYPE_REAL> class EXPORTGPUREG multiresRegistrationSolver 
  : public registrationSolver<REAL,ARRAY_TYPE_REAL>
{
public:
  
  multiresRegistrationSolver() : registrationSolver<REAL,ARRAY_TYPE_REAL>(){
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

  using registrationSolver<REAL,ARRAY_TYPE_REAL>::solve;
  
  virtual boost::shared_ptr<ARRAY_TYPE_REAL> solve
    ( ARRAY_TYPE_REAL *fixed_image, ARRAY_TYPE_REAL *moving_image, 
      bool input_normalization_allowed = false  )
  {
    // Some initial validity tests
    //
    
    if( !fixed_image || !moving_image ){
      this->solver_error( "multiresRegistrationSolver::solve : invalid input pointer" );
      return boost::shared_ptr<ARRAY_TYPE_REAL>();
    }

    if( !this->interpolator_.get() ){
      this->solver_error( "multiresRegistrationSolver::solve : interpolator not set" );
      return boost::shared_ptr<ARRAY_TYPE_REAL>();
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
    
    if( !normalize(normalized_fixed) || !normalize(normalized_moving) ){
      this->solver_error( "multiresRegistrationSolver::solve : failed to normalize input images" );
      return boost::shared_ptr<ARRAY_TYPE_REAL>();
    }   
    
    // Invoke multi-resolution solver
    // 
    
    if( this->output_mode_ >= registrationSolver<REAL,ARRAY_TYPE_REAL>::OUTPUT_VERBOSE ) {
      std::cout << std::endl << "Starting multiresolution registration " <<  std::endl;
    }
    
    boost::shared_ptr<ARRAY_TYPE_REAL> result = 
      solveMultiRes( num_multires_levels_, normalized_fixed, normalized_moving, this->stencil_.get() );
    
    return result;
  }
  
protected:

  // Pure virtual fuctions to be implemented in a subclass
  //

  virtual bool normalize( ARRAY_TYPE_REAL *image ) = 0;
  virtual bool compute( ARRAY_TYPE_REAL *fixed_image, ARRAY_TYPE_REAL *moving_image, ARRAY_TYPE_REAL *stencil_image,
			boost::shared_ptr<ARRAY_TYPE_REAL> &result ) = 0;
  virtual boost::shared_ptr<ARRAY_TYPE_REAL> downsample( ARRAY_TYPE_REAL *image ) = 0;
  virtual boost::shared_ptr<ARRAY_TYPE_REAL> upsample( ARRAY_TYPE_REAL *displacements ) = 0;  
  
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

      boost::shared_ptr<ARRAY_TYPE_REAL> fixed_image_lowres = downsample(fixed_image);
      boost::shared_ptr<ARRAY_TYPE_REAL> moving_image_lowres = downsample(moving_image);
      boost::shared_ptr<ARRAY_TYPE_REAL> stencil_image_lowres = 
	((stencil_image) ? downsample(stencil_image) : boost::shared_ptr<ARRAY_TYPE_REAL>());
      
      // Compute displacement field at the downsampled resolution
      //

      boost::shared_ptr<ARRAY_TYPE_REAL> result_lowres = 
	solveMultiRes( res_level-1, fixed_image_lowres.get(), moving_image_lowres.get(), stencil_image_lowres.get() );
      
      if( result_lowres.get() == 0x0 ){
	this->solver_error( "multiresRegistrationSolver::solveMultiRes : Recursive solver failed" );
	return boost::shared_ptr<ARRAY_TYPE_REAL>();
      }
      
      // Clean up low resolution image data
      //
      
      fixed_image_lowres.reset();
      moving_image_lowres.reset();
      stencil_image_lowres.reset();

      // Upsample lowres results to current resolution
      //
      
      result = upsample(result_lowres.get());
     
      // Clean up low resolution result
      //

      result_lowres.reset();

      // Some output to track our progress at runtime
      //

      if( this->output_mode_ >= registrationSolver<REAL,ARRAY_TYPE_REAL>::OUTPUT_VERBOSE ) {
	std::cout << std::endl << "Multiresolution level " << res_level;
      }
      
      // Use estimated (lowres) motion to compute displacements at the current resolution
      //

      if( !compute( fixed_image, moving_image, stencil_image, result )){
	this->solver_error( "multiresRegistrationSolver::solveMultiRes : compute failed (1)" );
	return boost::shared_ptr<ARRAY_TYPE_REAL>();
      }      
    }  
    else{
      
      // 
      // We are now at the end of the multi-resolution chain
      //
      
      // Some output to track our progress at runtime
      //
      
      if( this->output_mode_ >= registrationSolver<REAL,ARRAY_TYPE_REAL>::OUTPUT_VERBOSE ) {
	std::cout << std::endl << "Multiresolution level " << res_level << " (lowest)";
      }
      
      // Compute displacements at the current resolution (no estimate can be provided)
      //
      
      if( !compute( fixed_image, moving_image, stencil_image, result )){
	this->solver_error( "multiresRegistrationSolver::solveMultiRes : compute failed (2)" );
	return boost::shared_ptr<ARRAY_TYPE_REAL>();
      }      	
    }
    
    return result;
  }
  
protected:
  unsigned int num_multires_levels_;
  unsigned int max_num_iterations_per_level_;
};
