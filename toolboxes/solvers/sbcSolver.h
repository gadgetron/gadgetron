/*
  An implementation of the "Constrained CS Optimization Algorithm" of the paper
  "The Split Bregman Method for L1-Regularized Problems" by Tom Goldstein and Stanley Osher. 
  Siam J. Imaging Sciences. Vol. 2, No. 2, pp. 323–343.
*/

#pragma once

#include "sbSolver.h"

template <class REAL, class ELEMENT_TYPE, class ARRAY_TYPE_REAL, class ARRAY_TYPE_ELEMENT> 
class sbcSolver : public sbSolver<REAL, ELEMENT_TYPE, ARRAY_TYPE_REAL, ARRAY_TYPE_ELEMENT>
{
public:

  sbcSolver( int output_mode = solver<ARRAY_TYPE_ELEMENT>::OUTPUT_SILENT ) : sbSolver<REAL, ELEMENT_TYPE, ARRAY_TYPE_REAL, ARRAY_TYPE_ELEMENT>( output_mode ) { 
    this->tolerance_ = get_zero<REAL>();
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
  virtual boost::shared_ptr<ARRAY_TYPE_REAL> solver_norm( ARRAY_TYPE_ELEMENT* ) = 0;
  virtual boost::shared_ptr<ARRAY_TYPE_REAL> solver_norm_squared( ARRAY_TYPE_ELEMENT* ) = 0;
  virtual bool solver_shrink1( REAL, ARRAY_TYPE_ELEMENT*, ARRAY_TYPE_ELEMENT* ) = 0;
  virtual bool solver_shrinkd( REAL, ARRAY_TYPE_REAL*, ARRAY_TYPE_ELEMENT*, ARRAY_TYPE_ELEMENT* ) = 0;

  virtual boost::shared_ptr<ARRAY_TYPE_ELEMENT> solve( ARRAY_TYPE_ELEMENT *_f )
  {
    // Some tests to see if we are ready to go...
    //
    if( !this->inner_solver_.get() ){
      this->solver_error( "sbcSolver::solve : inner solver has not been set" );
      return boost::shared_ptr<ARRAY_TYPE_ELEMENT>();
    }

    if( !this->encoding_operator_.get() ){
      this->solver_error( "sbcSolver::solve : encoding operator has not been set" );
      return boost::shared_ptr<ARRAY_TYPE_ELEMENT>();
    }

    if( this->regularization_operators_.size() == 0 && this->regularization_group_operators_.size() == 0 ){
      this->solver_error( "sbcSolver::solve : at least one matrix regularizer must be added" );
      return boost::shared_ptr<ARRAY_TYPE_ELEMENT>();
    }

    if( !this->image_dims_.get() ){
      this->solver_error( "sbcSolver::solve : image dimensions have not been set" );
      return boost::shared_ptr<ARRAY_TYPE_ELEMENT>();
    }

    // Make a copy of _f - we are going to "normalize" the input shortly
    //
    ARRAY_TYPE_ELEMENT f(*_f);
    if( !f.get_data_ptr() ){
      this->solver_error( "sbcSolver::solve : memory allocation of f failed" );
      return boost::shared_ptr<ARRAY_TYPE_ELEMENT>();
    }

    // Initialize u_k to E^H f 
    //
    boost::shared_ptr<ARRAY_TYPE_ELEMENT> u_k = boost::shared_ptr<ARRAY_TYPE_ELEMENT>(new ARRAY_TYPE_ELEMENT());
    if( u_k.get() ) u_k->create( this->image_dims_.get() );

    if( !u_k.get() || !u_k->get_data_ptr() ){
      this->solver_error( "sbcSolver::solve : memory allocation of u_k failed" );
      return boost::shared_ptr<ARRAY_TYPE_ELEMENT>();
    }    

    if( this->encoding_operator_->mult_MH( &f, u_k.get() ) < 0 ){
      this->solver_error( "sbcSolver::solve : adjoint encoding operation failed on f" );
      return boost::shared_ptr<ARRAY_TYPE_ELEMENT>();
    }

    // Normalize u_k and f
    //
    ELEMENT_TYPE image_scale=get_one<ELEMENT_TYPE>();
    /*    {
      // Normalize to an average energy of "one intensity unit per image element"
      REAL sum = solver_asum( u_k.get() );
      image_scale = mul<REAL>(( (REAL) (u_k->get_number_of_elements())/sum), get_one<ELEMENT_TYPE>() );
      solver_scal( image_scale, u_k.get() );
      solver_scal( image_scale, &f );
      }*/

    // Initialize f_k
    //
    ARRAY_TYPE_ELEMENT f_k(f);
    if( !f_k.get_data_ptr() ){
      this->solver_error( "sbcSolver::solve : memory allocation of f_k failed" );
      return boost::shared_ptr<ARRAY_TYPE_ELEMENT>();
    }

    // Initialize d_k and b_k arrays to zero images
    //
    //

    unsigned int num_reg_operators = this->regularization_operators_.size();
    for( unsigned int i=0; i<this->regularization_group_operators_.size(); i++ )
      num_reg_operators += this->regularization_group_operators_[i].size();

    boost::shared_array< boost::shared_ptr<ARRAY_TYPE_ELEMENT> > d_k( new boost::shared_ptr<ARRAY_TYPE_ELEMENT>[num_reg_operators] );
    boost::shared_array< boost::shared_ptr<ARRAY_TYPE_ELEMENT> > b_k( new boost::shared_ptr<ARRAY_TYPE_ELEMENT>[num_reg_operators] );

    if( !d_k.get() || !b_k.get() ){
      this->solver_error( "sbcSolver::solve : memory allocation of d_k or b_k failed" );
      return boost::shared_ptr<ARRAY_TYPE_ELEMENT>();
    }

    for( unsigned int i=0; i<num_reg_operators; i++ ){

      d_k[i] = boost::shared_ptr<ARRAY_TYPE_ELEMENT>(new ARRAY_TYPE_ELEMENT());

      if( d_k[i] ) d_k[i]->create( this->image_dims_.get() );

      if( !d_k[i]->get_data_ptr() ){
	this->solver_error( "sbSolver::solve : memory allocation of d_k failed" );
	return boost::shared_ptr<ARRAY_TYPE_ELEMENT>();
      }

      if( !solver_clear_element( d_k[i].get() )){
	this->solver_error( "sbSolver::solve : failed to clear internal memory buffer d_k" );
	return boost::shared_ptr<ARRAY_TYPE_ELEMENT>();
      }

      b_k[i] = boost::shared_ptr<ARRAY_TYPE_ELEMENT>(new ARRAY_TYPE_ELEMENT());

      if( b_k[i] ) b_k[i]->create( this->image_dims_.get() );

      if( !b_k[i]->get_data_ptr() ){
	this->solver_error( "sbSolver::solve : memory allocation of b_k failed" );
	return boost::shared_ptr<ARRAY_TYPE_ELEMENT>();
      }
      
      if( !solver_clear_element( b_k[i].get() )){
	this->solver_error( "sbSolver::solve : failed to clear internal memory buffer b_k" );
	return boost::shared_ptr<ARRAY_TYPE_ELEMENT>();
      } 
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
      if( this->encoding_operator_->mult_MH( &f_k, muEHf.get() ) < 0 ){
	this->solver_error( "sbcSolver::solve : adjoint encoding operation failed on f" );
	return boost::shared_ptr<ARRAY_TYPE_ELEMENT>();
      }

      // Scale EHf with mu
      //
      if( !solver_scal( mul<REAL>(this->encoding_operator_->get_weight(), get_one<ELEMENT_TYPE>()), muEHf.get() ) ){
	this->solver_error( "sbcSolver::solve : error scaling EHf with mu" );
	return boost::shared_ptr<ARRAY_TYPE_ELEMENT>();
      }
      
      // Invoke the core solver
      //
      if( !solve_core( this->tolerance_, this->inner_iterations_, 1, f_k, *muEHf.get(), u_k, d_k, b_k ) ){
	this->solver_error( "sbcSolver::solve : core solver failed" );
	return boost::shared_ptr<ARRAY_TYPE_ELEMENT>();
      } 

      // Update f_k
      //
      ARRAY_TYPE_ELEMENT encoded_image;
      if( encoded_image.create( f.get_dimensions().get() ) < 0 ){
        this->solver_error( "sbcSolver::solve : memory allocation for encoded image failed" );
        return boost::shared_ptr<ARRAY_TYPE_ELEMENT>();      
      }

      if( this->encoding_operator_->mult_M( u_k.get(), &encoded_image ) < 0 ){
        this->solver_error( "sbSolver::solve : computation of encoded image failed" );
        return boost::shared_ptr<ARRAY_TYPE_ELEMENT>();
      }

      if( !this->solver_axpy_element( get_zero<ELEMENT_TYPE>()-get_one<ELEMENT_TYPE>(), &f, &encoded_image )){ // notice the deliberate sign manipulation
        this->solver_error( "sbSolver::solve : computation of update argument to f_k failed" );
        return boost::shared_ptr<ARRAY_TYPE_ELEMENT>();
      }

      if( this->tolerance_ > get_zero<REAL>() || this->output_mode_ >= solver<ARRAY_TYPE_ELEMENT>::OUTPUT_VERBOSE ){
	
	REAL delta = cuNDA_asum<REAL>(&encoded_image);
	
	if( this->output_mode_ >= solver<ARRAY_TYPE_ELEMENT>::OUTPUT_VERBOSE )
	  std::cout << std::endl << "u_k delta (outer loop): " << delta << std::endl << std::endl;

	if( delta < this->tolerance_ )
	  break;
      }
      
      if( !solver_axpy_element( get_zero<ELEMENT_TYPE>()-get_one<ELEMENT_TYPE>(), &encoded_image, &f_k )){ // notice the deliberate sign manipulation
	this->solver_error( "sbSolver::solve : computation of update argument to f_k failed" );
	return boost::shared_ptr<ARRAY_TYPE_ELEMENT>();
      }
      
    } // end of outer loop
    
    // Undo intermediate scaling of u_k ...
    solver_scal( reciprocal<ELEMENT_TYPE>(image_scale), u_k.get() );
    
    // ... and return the result
    //
    
    return u_k;
  }

};
