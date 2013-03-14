/*
  An implementation of the "Constrained CS Optimization Algorithm" of the paper
  "The Split Bregman Method for L1-Regularized Problems" by Tom Goldstein and Stanley Osher. 
  Siam J. Imaging Sciences. Vol. 2, No. 2, pp. 323-343.
*/

#pragma once

#include "sbSolver.h"

namespace Gadgetron{
template<
	  class ARRAY_TYPE_REAL, 
	  class ARRAY_TYPE_ELEMENT, 
	  class INNER_SOLVER>
class sbcSolver 
  : public sbSolver<ARRAY_TYPE_REAL, ARRAY_TYPE_ELEMENT, INNER_SOLVER>
{
  	typedef typename ARRAY_TYPE_ELEMENT::element_type ELEMENT_TYPE;
		 typedef typename realType<ELEMENT_TYPE>::type REAL;
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
    boost::shared_ptr<ARRAY_TYPE_ELEMENT> u_k( new ARRAY_TYPE_ELEMENT(this->encoding_operator_->get_domain_dimensions()) );


    // Use x0 (if provided) as starting estimate
    if( this->get_x0().get() )
      *u_k = *(this->get_x0());
    else
    	u_k->clear();

    //this->get_inner_solver()->set_x0( u_k );

    // Normalize (a copy of) the input data
    //

    boost::shared_ptr<ARRAY_TYPE_ELEMENT> f(new ARRAY_TYPE_ELEMENT(*_f));



    REAL image_scale;
    this->normalize( f, image_scale );
    
    // Initialize f_k
    //

    boost::shared_ptr<ARRAY_TYPE_ELEMENT> f_k( new ARRAY_TYPE_ELEMENT(*f) );


    // Initialze d, b, and p_M
    //

    boost::shared_array< boost::shared_ptr<ARRAY_TYPE_ELEMENT> > d_k, b_k, p_M;
    
    this->initialize( image_scale, d_k, b_k, p_M );
        
    // Outer loop
    //

    for( unsigned int outer_iteration=0; outer_iteration<this->outer_iterations_; outer_iteration++ ) {
      
      // Invoke the core solver
      //

      this->core( this->tolerance_, this->inner_iterations_, 1, f_k, u_k, d_k, b_k, p_M );

      // Update f_k
      //

      ARRAY_TYPE_ELEMENT encoded_image(f->get_dimensions());

      this->encoding_operator_->mult_M( u_k.get(), &encoded_image );
      encoded_image -= *f;


      if( this->tolerance_ > REAL(0) || this->output_mode_ >= solver<ARRAY_TYPE_ELEMENT, ARRAY_TYPE_ELEMENT>::OUTPUT_VERBOSE ){
	
      	REAL delta = asum(&encoded_image);
	
      	if( this->output_mode_ >= solver<ARRAY_TYPE_ELEMENT, ARRAY_TYPE_ELEMENT>::OUTPUT_VERBOSE )
      		std::cout << "Residual (outer loop): " << delta << std::endl << std::endl;

      	if( delta < this->tolerance_ )
      		break;
      }
      *f_k -=  encoded_image;
      
    } // end of outer loop
    
    // Undo the intermediate scaling of u_k ...
    //

    this->undo_normalization( u_k, 1/(image_scale) );
    
    // Clean up memory occupied by the operator container and inner solver
    this->deinitialize();
    
    // ... and return the result
    //    

    return u_k;
  }  
};
}
