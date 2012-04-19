/*
  An implementation of the "Generalized Landweber Solver" based on the paper
  Acceleration and Filtering in the Generalized Landweber Iteration Using a Variable Shaping Matrix
  by Pan, Yagle, Clinthorne, and Rogers. IEEE TMI 1993; 12(2):278-286.
*/

#pragma once

#include "linearSolver.h"

#include <vector>
#include <iostream>

template <class REAL, class ELEMENT_TYPE, class ARRAY_TYPE> class lwSolver 
  : public linearSolver<REAL, ELEMENT_TYPE, ARRAY_TYPE>
{
  
public:

  // Constructor
  lwSolver() : linearSolver<REAL, ELEMENT_TYPE, ARRAY_TYPE>() { 
    iterations_ = 10;
  }
  
  // Destructor
  virtual ~lwSolver() {}
  
  // Set/get maximally allowed number of iterations
  virtual void set_max_iterations( unsigned int iterations ) { iterations_ = iterations; }
  virtual unsigned int get_max_iterations() { return iterations_; }  

  // Pure virtual functions defining core solver functionality
  // Implemented on the host/device respectively in a derived class

  virtual bool solver_clear( ARRAY_TYPE* ) = 0;
  virtual bool solver_axpy( ELEMENT_TYPE, ARRAY_TYPE*, ARRAY_TYPE* ) = 0;

  // Inherited solver interface
  virtual boost::shared_ptr<ARRAY_TYPE> solve( ARRAY_TYPE *b )
  {
   
    // Initial validity checks
    //

    std::vector<unsigned int> image_dims = this->encoding_operator_->get_domain_dimensions();

    if( image_dims.size() == 0 ){
      this->solver_error( "Error: lwSolver::solve : domain dimensions not set on encoding operator" );
      return boost::shared_ptr<ARRAY_TYPE>();
    }
        
    // Allocate solution array.
    // Clear or set to x0 if provided
    //

    boost::shared_ptr<ARRAY_TYPE> x( new ARRAY_TYPE() );
    if( this->get_x0().get() )
      *x = *(this->get_x0());
    else{
      if( !x->create( &image_dims )){
	this->solver_error( "Error: lwSolver::solve : failed to allocate solution" );
	return boost::shared_ptr<ARRAY_TYPE>();
      }
      if( !solver_clear( x.get() )){
	this->solver_error( "Error: lwSolver::solve : failed to clear solution" );
	return boost::shared_ptr<ARRAY_TYPE>();	
      }      
    }
    
    // Allocate some temporary storage
    ARRAY_TYPE tmp_M, tmp_MH, tmp_MH_M_acc;
    if( !tmp_M.create( b.get_dimensions().get() ) || 
	!tmp_MH.create( &image_dims ) || 
	!tmp_acc.create( &image_dims ) ){
      this->solver_error( "Error: lwSolver::solve : failed to allocate temporary storage" );
      return boost::shared_ptr<ARRAY_TYPE>();
    }

    /*    
    // Compute adjoint of encoding operator
    if( this->encoding_operator_->mult_MH( d, &tmp_MH ) < 0 ){
      this->solver_error( "Error: lwSolver::solve : failed to compute adjoint of encoding operator" );
      return boost::shared_ptr<ARRAY_TYPE>();
    }
    */

    // Main solver iteration loop
    //
    
    for( unsigned int iteration=0; iteration<iterations_; iteration++ ){

      // Clear accumulation buffer
      if( !solver_clear( &tmp_acc )){
	this->solver_error( "Error: lwSolver::solve : failed to clear accumulation buffer" );
	return boost::shared_ptr<ARRAY_TYPE>();	
      }      
      
      // Apply encoding operator to current solution
      if( this->encoding_operator_->mult_M( x.get(), &tmp_M ) < 0 ) {
	this->solver_error( "Error: lwSolver::solve : failed to apply encoding operator" );
	return boost::shared_ptr<ARRAY_TYPE>();
      }
      
      // Accumulate
      if( !solver_axpy(-this->encoding_operator_->get_weight(), &tmp_MH_M, &tmp_acc )) {
	this->solver_error( "Error: lwSolver::solve : failed to accumulate (1)" );
	return boost::shared_ptr<ARRAY_TYPE>();
      }

      // Loop over operators
      for( unsigned int i=0; i<this->regularization_operators_.size(); i++){

	// Compute encoding operator adjoint
	if( this->regularization_operators_[i]->mult_MH_M( x.get(), &tmp_MH_M ) < 0 ) {
	  this->solver_error( "Error: lwSolver::solve : failed to apply operator" );
	  return boost::shared_ptr<ARRAY_TYPE>();
	}
	
	// Accumulate
	if( !solver_axpy(-this->regularization_operators_[i]->get_weight(), &tmp_MH_M, &tmp_acc )) {
	  this->solver_error( "Error: lwSolver::solve : failed to accumulate (2)" );
	  return boost::shared_ptr<ARRAY_TYPE>();
	}
      }
      
      // Update x
      if( !solver_axpy(REAL(1), &tmp_acc, x.get() )) {
	this->solver_error( "Error: lwSolver::solve : failed to update x" );
	return boost::shared_ptr<ARRAY_TYPE>();
      }      
    }
    
    return x;
  }
  
protected:

  // Maximum number of iterations
  unsigned int iterations_;
};
