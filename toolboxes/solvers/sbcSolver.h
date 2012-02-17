/*
  An implementation of the "Constrained CS Optimization Algorithm" of the paper
  "The Split Bregman Method for L1-Regularized Problems" by Tom Goldstein and Stanley Osher. 
  Siam J. Imaging Sciences. Vol. 2, No. 2, pp. 323�343.
*/

#pragma once

#include "sbSolver.h"
//#include "ndarray_vector_td_utilities.h"

template <class REAL, class ELEMENT_TYPE, class ARRAY_TYPE_REAL, class ARRAY_TYPE_ELEMENT> 
class sbcSolver : public sbSolver<REAL, ELEMENT_TYPE, ARRAY_TYPE_REAL, ARRAY_TYPE_ELEMENT>
{
public:

  sbcSolver( int output_mode = solver<ARRAY_TYPE_ELEMENT, ARRAY_TYPE_ELEMENT>::OUTPUT_SILENT ) : sbSolver<REAL, ELEMENT_TYPE, ARRAY_TYPE_REAL, ARRAY_TYPE_ELEMENT>( output_mode ) { 
    this->tolerance_ = REAL(0);
    this->outer_iterations_ = 10;
    this->inner_iterations_ = 5;
  }

  virtual ~sbcSolver() {}

  virtual bool solver_clear_real( ARRAY_TYPE_REAL* ) = 0;
  virtual bool solver_clear_element( ARRAY_TYPE_ELEMENT* ) = 0;
  virtual bool solver_sqrt( ARRAY_TYPE_REAL* ) = 0;
  virtual bool solver_scal( ELEMENT_TYPE, ARRAY_TYPE_ELEMENT* ) = 0;
  virtual bool solver_axpy_real( REAL, ARRAY_TYPE_REAL*, ARRAY_TYPE_REAL* ) = 0;
  virtual bool solver_axpy_element( ELEMENT_TYPE, ARRAY_TYPE_ELEMENT*, ARRAY_TYPE_ELEMENT* ) = 0;
  virtual REAL solver_asum( ARRAY_TYPE_ELEMENT* ) = 0;
  virtual boost::shared_ptr<ARRAY_TYPE_REAL> solver_abs( ARRAY_TYPE_ELEMENT* ) = 0;
  virtual boost::shared_ptr<ARRAY_TYPE_REAL> solver_norm( ARRAY_TYPE_ELEMENT* ) = 0;
  virtual bool solver_shrink1( REAL, ARRAY_TYPE_ELEMENT*, ARRAY_TYPE_ELEMENT* ) = 0;
  virtual bool solver_shrinkd( REAL, ARRAY_TYPE_REAL*, ARRAY_TYPE_ELEMENT*, ARRAY_TYPE_ELEMENT* ) = 0;

  virtual boost::shared_ptr<ARRAY_TYPE_ELEMENT> solve( ARRAY_TYPE_ELEMENT *_f )
  {
    // Check if everything is set up right
    //
    if( !this->validate() ){
      this->solver_error( "sbcSolver::solve : setup failed validation");
      return boost::shared_ptr<ARRAY_TYPE_ELEMENT>();
    }
    
    // Initialze d and b
    //
    boost::shared_array< boost::shared_ptr<ARRAY_TYPE_ELEMENT> > d_k;
    boost::shared_array< boost::shared_ptr<ARRAY_TYPE_ELEMENT> > b_k;
    if( !this->initialize( d_k, b_k ) ){
      this->solver_error( "sbcSolver::solve : failed to initialize");
      return boost::shared_ptr<ARRAY_TYPE_ELEMENT>();
    }
    
    // Make a copy of _f before normalization
    //
    boost::shared_ptr<ARRAY_TYPE_ELEMENT> f( new ARRAY_TYPE_ELEMENT(*_f) );
    if( !f->get_data_ptr() ){
      this->solver_error( "sbcSolver::solve : memory allocation of f failed" );
      return boost::shared_ptr<ARRAY_TYPE_ELEMENT>();
    }

    // Define u_k
    //
    boost::shared_ptr<ARRAY_TYPE_ELEMENT> u_k;

    // Normalize the data
    //
    ELEMENT_TYPE image_scale = ELEMENT_TYPE(0);
    if( !this->normalize( f, u_k, image_scale ) ){
      this->solver_error( "sbcSolver::solve : normalization failed" );
      return boost::shared_ptr<ARRAY_TYPE_ELEMENT>();
    }
    
    // Initialize f_k
    //
    boost::shared_ptr<ARRAY_TYPE_ELEMENT> f_k( new ARRAY_TYPE_ELEMENT(*f) );
    if( !f_k->get_data_ptr() ){
      this->solver_error( "sbcSolver::solve : memory allocation of f_k failed" );
      return boost::shared_ptr<ARRAY_TYPE_ELEMENT>();
    }

    boost::shared_ptr<ARRAY_TYPE_ELEMENT> muEHf(new ARRAY_TYPE_ELEMENT());
    if( muEHf.get() ) muEHf->create( this->image_dims_.get() );
    
    if( !muEHf.get() || !muEHf->get_data_ptr() ){
      this->solver_error( "sbcSolver::solve : memory allocation of muEHf failed" );
      return boost::shared_ptr<ARRAY_TYPE_ELEMENT>();
    }
    
    // Outer loop
    //
    for( unsigned int outer_iteration=0; outer_iteration<this->outer_iterations_; outer_iteration++ ) {

      // Keep a copy of E^H f for subsequent rhs computations in the inner solver
      //             
      if( this->encoding_operator_->mult_MH( f_k.get(), muEHf.get() ) < 0 ){
	this->solver_error( "sbcSolver::solve : adjoint encoding operation failed on f_k" );
	return boost::shared_ptr<ARRAY_TYPE_ELEMENT>();
      }

      // Scale EHf with mu
      //
      if( !solver_scal( REAL(this->encoding_operator_->get_weight()), muEHf.get() ) ){
	this->solver_error( "sbcSolver::solve : error scaling EHf with mu" );
	return boost::shared_ptr<ARRAY_TYPE_ELEMENT>();
      }
      
      // Invoke the core solver
      //
      if( !core( this->tolerance_, this->inner_iterations_, 1, f_k, muEHf, u_k, d_k, b_k ) ){
      	this->solver_error( "sbcSolver::solve : core solver failed" );
      	return boost::shared_ptr<ARRAY_TYPE_ELEMENT>();
      } 

      // Update f_k
      //
      ARRAY_TYPE_ELEMENT encoded_image;
      if( encoded_image.create( f->get_dimensions().get() ) < 0 ){
        this->solver_error( "sbcSolver::solve : memory allocation for encoded image failed" );
        return boost::shared_ptr<ARRAY_TYPE_ELEMENT>();      
      }

      if( this->encoding_operator_->mult_M( u_k.get(), &encoded_image ) < 0 ){
        this->solver_error( "sbcSolver::solve : computation of encoded image failed" );
        return boost::shared_ptr<ARRAY_TYPE_ELEMENT>();
      }

      if( !this->solver_axpy_element( ELEMENT_TYPE(0)-ELEMENT_TYPE(1), f.get(), &encoded_image )){ // notice the deliberate sign manipulation
        this->solver_error( "sbcSolver::solve : computation of update argument to f_k failed" );
        return boost::shared_ptr<ARRAY_TYPE_ELEMENT>();
      }

      if( this->tolerance_ > REAL(0) || this->output_mode_ >= solver<ARRAY_TYPE_ELEMENT, ARRAY_TYPE_ELEMENT>::OUTPUT_VERBOSE ){
	
	REAL delta = cuNDA_asum<REAL>(&encoded_image);
	
	if( this->output_mode_ >= solver<ARRAY_TYPE_ELEMENT, ARRAY_TYPE_ELEMENT>::OUTPUT_VERBOSE )
	  std::cout << std::endl << "u_k delta (outer loop): " << delta << std::endl << std::endl;

	//static REAL min_delta = (REAL) 1e9;
	//if( delta < min_delta ){
	//  min_delta = delta;
	// boost::shared_ptr<ARRAY_TYPE_REAL > outm = solver_abs( u_k.get() );
	// boost::shared_ptr<hoNDArray<REAL> > out = outm->to_host();
	//  write_nd_array( out.get(), "cur_min.real");
	//}

	if( delta < this->tolerance_ )
	  break;
      }
      
      if( !solver_axpy_element( ELEMENT_TYPE(0)-ELEMENT_TYPE(1), &encoded_image, f_k.get() )){ // notice the deliberate sign manipulation
	this->solver_error( "sbcSolver::solve : computation of update argument to f_k failed" );
	return boost::shared_ptr<ARRAY_TYPE_ELEMENT>();
      }
      
    } // end of outer loop
    
    // Undo the intermediate scaling of u_k ...
    //
    if( !undo_normalization( u_k, 1/(image_scale) )){
      this->solver_error( "sbcSolver::solve : unable to undo normalization" );
      return boost::shared_ptr<ARRAY_TYPE_ELEMENT>();
    } 
    
    // ... and return the result
    //    
    return u_k;
  }  
};
