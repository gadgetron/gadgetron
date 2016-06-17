/*
  An implementation of the constrained solver of the paper
  "The Split Bregman Method for L1-Regularized Problems" by Tom Goldstein and Stanley Osher. 
  Siam J. Imaging Sciences. Vol. 2, No. 2, pp. 323-343.
*/

#pragma once

#include "sbSolver.h"

namespace Gadgetron{

  template<class ARRAY_TYPE_REAL,
	   class ARRAY_TYPE_ELEMENT, 
	   class INNER_SOLVER>
  class sbcSolver : public sbSolver<ARRAY_TYPE_REAL, ARRAY_TYPE_ELEMENT, INNER_SOLVER>
  {
  protected:
    typedef typename ARRAY_TYPE_REAL::element_type REAL;
    
  public:
  
    sbcSolver() : sbSolver<ARRAY_TYPE_REAL, ARRAY_TYPE_ELEMENT, INNER_SOLVER>() {}
    virtual ~sbcSolver() {}
    
    virtual boost::shared_ptr<ARRAY_TYPE_ELEMENT> solve( ARRAY_TYPE_ELEMENT *_f )
    {
      // Check if everything is set up right
      //
      this->validate_solver();

      // Define u_k
      //
      boost::shared_ptr<ARRAY_TYPE_ELEMENT> u_k( new ARRAY_TYPE_ELEMENT(this->encoding_operator_->get_domain_dimensions()));

      // Use x0 (if provided) as starting estimate
      if(this->get_x0().get())
	*u_k = *(this->get_x0());
      else 
	clear(u_k.get());


      // Normalize and _then_ initialize (the order matters)
      //
      
      boost::shared_ptr<ARRAY_TYPE_ELEMENT> f(new ARRAY_TYPE_ELEMENT(*_f));
      REAL normalization_factor = this->normalize_data( f.get() );
      boost::shared_ptr<ARRAY_TYPE_ELEMENT> f_k(new ARRAY_TYPE_ELEMENT(*f));
      this->initialize( normalization_factor );
        
      // Outer loop
      //

      for( unsigned int outer_iteration=0; outer_iteration<this->outer_iterations_; outer_iteration++ ) {
      
	if( this->output_mode_ >= solver<ARRAY_TYPE_ELEMENT, ARRAY_TYPE_ELEMENT>::OUTPUT_VERBOSE )
	  GDEBUG_STREAM(std::endl << "SBC outer loop iteration " << outer_iteration << std::endl << std::endl);
	
	// Invoke the core solver
	//
	
	this->core( this->tolerance_, this->inner_iterations_, 1, f_k, u_k );

	// Update f_k
	//

	ARRAY_TYPE_ELEMENT encoded_image(f->get_dimensions());
	this->encoding_operator_->mult_M( u_k.get(), &encoded_image );
	encoded_image -= *f;

	if( this->tolerance_ > REAL(0) || this->output_mode_ >= solver<ARRAY_TYPE_ELEMENT, ARRAY_TYPE_ELEMENT>::OUTPUT_VERBOSE ){
	
	  REAL delta = nrm2(&encoded_image);
	
	  if( this->output_mode_ >= solver<ARRAY_TYPE_ELEMENT, ARRAY_TYPE_ELEMENT>::OUTPUT_VERBOSE )
	    GDEBUG_STREAM("Squared residual norm (outer loop): " << delta*delta << std::endl << std::endl);
	  
	  if( delta < this->tolerance_ )
	    break;
	}

	*f_k -= encoded_image;
      
      } // end of outer loop
        
      // Clean up memory occupied by the operator container and inner solver
      this->deinitialize();

      // Undo normalization
      *u_k /= normalization_factor;
      
      // ... and return the result
      return u_k;
    }  
  };
}
