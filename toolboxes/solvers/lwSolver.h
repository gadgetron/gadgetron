/*
  An implementation of the "Generalized Landweber Solver" based on the paper
  "Theory and methods related to the singular-function expansion and Landweber's iteration..."
  by O.N. Strand, Siam J Numer, Anal. 1974;11(4):798-825.
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
    iterations_ = 3;
    alpha_ = REAL(1);
  }
  
  // Destructor
  virtual ~lwSolver() {}
  
  // Set/get maximally allowed number of iterations
  virtual void set_max_iterations( unsigned int iterations ) { iterations_ = iterations; }
  virtual unsigned int get_max_iterations() { return iterations_; }  

  // Set/get alpha.
  // Optimally set alpha to 1/(sigma^2), sigma being the largest singular value of the "sum of operators"
  virtual void set_alpha( REAL alpha ) { alpha_ = alpha; }
  virtual REAL get_alpha() { return alpha_; }  

  // Pure virtual functions defining core solver functionality
  // Implemented on the host/device respectively in a derived class

  virtual bool solver_clear( ARRAY_TYPE* ) = 0;
  virtual bool solver_scale( REAL, ARRAY_TYPE* ) = 0;
  virtual bool solver_axpy( ELEMENT_TYPE, ARRAY_TYPE*, ARRAY_TYPE* ) = 0;
  virtual REAL solver_asum( ARRAY_TYPE* ) = 0;

  // Inherited solver interface
  virtual boost::shared_ptr<ARRAY_TYPE> solve( ARRAY_TYPE *b )
  {   
    // Initial validity checks
    //

    std::vector<unsigned int> image_dims = *this->encoding_operator_->get_domain_dimensions();

    if( image_dims.size() == 0 ){
      this->solver_error( "Error: lwSolver::solve : domain dimensions not set on encoding operator" );
      return boost::shared_ptr<ARRAY_TYPE>();
    }
        
    // Allocate solution array.
    // Clear or set to x0 if provided
    //

    boost::shared_ptr<ARRAY_TYPE> x( new ARRAY_TYPE() );
    if( this->get_x0().get() ){
      *x = *(this->get_x0());
    }
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

    ARRAY_TYPE x_prev;

    // Main solver iteration loop
    //
    
    for( unsigned int iteration=0; iteration<iterations_; iteration++ ){
      
      // Keep previous x for convergence reporting
      // 

      if( this->output_mode_ >= solver<ARRAY_TYPE, ARRAY_TYPE>::OUTPUT_VERBOSE ){
	x_prev = *x;
	if( !x_prev.get_data_ptr() ){
	  this->solver_error( "Error: sbSolver::core : assignment to x_prev failed" );
	  return boost::shared_ptr<ARRAY_TYPE>();	
	}
      }
      
      // Compute residual image, i.e. A^T(b-Ax_k)
      //
      
      boost::shared_ptr<ARRAY_TYPE> r = compute_residual_image( x.get(), b );

      if( !r.get() ){
	this->solver_error( "Error: lwSolver::solve : failed to compute residual image" );
	return boost::shared_ptr<ARRAY_TYPE>();	
      }      

      // Multiply residual with shaping matrix
      //

      boost::shared_ptr<ARRAY_TYPE> rr = apply_shaping_matrix( r.get() );
      
      if( !rr.get() ){
	this->solver_error( "Error: lwSolver::solve : failed to apply shaping matrix" );
	return boost::shared_ptr<ARRAY_TYPE>();	
      }      

      // Update x
      //

      if( !solver_axpy( get_alpha(), rr.get(), x.get() )) {
	this->solver_error( "Error: lwSolver::solve : failed to update x" );
	return boost::shared_ptr<ARRAY_TYPE>();
      }
      
      if( this->output_mode_ >= solver<ARRAY_TYPE, ARRAY_TYPE>::OUTPUT_VERBOSE ){
	if( !solver_axpy( ELEMENT_TYPE(-1), x.get(), &x_prev )){
	  this->solver_error( "Error: sbSolver::core : error computing delta x" );
	  return boost::shared_ptr<ARRAY_TYPE>();
	}
	std::cout << " iteration: " << iteration << ", delta x: " << solver_asum(&x_prev) << std::endl;
      }      
    }
    
    return x;
  }

protected:
  virtual boost::shared_ptr<ARRAY_TYPE> compute_residual_image( ARRAY_TYPE *x, ARRAY_TYPE *b )
  {    
    // Allocate some temporary storage and the esult array
    //
    
    ARRAY_TYPE tmp_M, tmp_acc;
    boost::shared_ptr<ARRAY_TYPE> res = boost::shared_ptr<ARRAY_TYPE>(new ARRAY_TYPE);

    if( !res->create( x->get_dimensions().get() )) {
      this->solver_error( "Error: lwSolver::compute_residual_image : failed to allocate result" );
      return boost::shared_ptr<ARRAY_TYPE>();
    }
    
    if( !tmp_M.create( b->get_dimensions().get() ) || 
	!tmp_acc.create( b->get_dimensions().get() ) ){
      this->solver_error( "Error: lwSolver::compute_residual_image : failed to allocate temporary storage" );
      return boost::shared_ptr<ARRAY_TYPE>();
    }
        
    // Clear accumulation buffer to b
    tmp_acc = *b;
    
    // Apply encoding operator to current solution
    if( this->encoding_operator_->mult_M( x, &tmp_M ) < 0 ) {
      this->solver_error( "Error: lwSolvercompute_residual_image : failed to apply encoding operator" );
      return boost::shared_ptr<ARRAY_TYPE>();
    }
    
    // Find residual
    if( !solver_axpy(REAL(-1), &tmp_M, &tmp_acc )) {
      this->solver_error( "Error: lwSolver::compute_residual_image : failed to accumulate (1)" );
      return boost::shared_ptr<ARRAY_TYPE>();
    }    
    
    // Adjoint residual    
    if( this->encoding_operator_->mult_MH( &tmp_acc, res.get() ) < 0 ) {
      this->solver_error( "Error: lwSolvercompute_residual_image : failed to apply adjoint encoding operator" );
      return boost::shared_ptr<ARRAY_TYPE>();
    }
 
    return res;
  }
  
  boost::shared_ptr<ARRAY_TYPE> apply_shaping_matrix( ARRAY_TYPE *r )
  {
    //
    // Apply 6th order polynomial F(lambda) -- see paper referenced at top
    //
    
    // The input residual r is modified (it is an internal implementation variable anyway)
    //

    // Memory allocation
    std::vector<unsigned int> image_dims = *this->encoding_operator_->get_domain_dimensions();
    boost::shared_ptr<ARRAY_TYPE> res = boost::shared_ptr<ARRAY_TYPE>(new ARRAY_TYPE());
    if( !res->create( &image_dims ) ){
      this->solver_error( "Error: lwSolver::apply_shaping_matrix : failed to allocate storage" );
      return boost::shared_ptr<ARRAY_TYPE>();
    }
    
    // Handle 0th order   
    *res = *r;
    if( !solver_scale( REAL(31.5), res.get() ) ){
      this->solver_error( "Error: lwSolver::apply_shaping_matrix : failed to initialize result" );
      return boost::shared_ptr<ARRAY_TYPE>();
    }

    // Handle 1th order
    if( !apply_shape_matrix_mult_MH_M( r, res.get(), REAL(-315) )){
      this->solver_error( "Error: lwSolver::apply_shaping_matrix : 1th order computation failed" );
      return boost::shared_ptr<ARRAY_TYPE>();
    }
    
    // Handle 2th order
    if( !apply_shape_matrix_mult_MH_M( r, res.get(), REAL(1443.75) )){
      this->solver_error( "Error: lwSolver::apply_shaping_matrix : 2th order computation failed" );
      return boost::shared_ptr<ARRAY_TYPE>();
    }    

    // Handle 3th order
    if( !apply_shape_matrix_mult_MH_M( r, res.get(), REAL(-3465) )){
      this->solver_error( "Error: lwSolver::apply_shaping_matrix : 3th order computation failed" );
      return boost::shared_ptr<ARRAY_TYPE>();
    }
    
    // Handle 4th order
    if( !apply_shape_matrix_mult_MH_M( r, res.get(), REAL(4504.5) )){
      this->solver_error( "Error: lwSolver::apply_shaping_matrix : 4th order computation failed" );
      return boost::shared_ptr<ARRAY_TYPE>();
    }
    
    // Handle 5th order
    if( !apply_shape_matrix_mult_MH_M( r, res.get(), REAL(-3003) )){
      this->solver_error( "Error: lwSolver::apply_shaping_matrix : 5th order computation failed" );
      return boost::shared_ptr<ARRAY_TYPE>();
    }
    
    // Handle 6th order
    if( !apply_shape_matrix_mult_MH_M( r, res.get(), REAL(804.375) )){
      this->solver_error( "Error: lwSolver::apply_shaping_matrix : 6th order computation failed" );
      return boost::shared_ptr<ARRAY_TYPE>();
    }
    
    // Return result
    return res;
  }

  bool apply_shape_matrix_mult_MH_M( ARRAY_TYPE *r, ARRAY_TYPE *acc, REAL w )
  {
    // Temporary storage
    std::vector<unsigned int> image_dims = *this->encoding_operator_->get_domain_dimensions();
    ARRAY_TYPE tmp_MH_M, tmp_acc;
    if( !tmp_MH_M.create( &image_dims ) ||
	!tmp_acc.create( &image_dims ) ){
      this->solver_error( "Error: lwSolver::apply_shaping_matrix_mult_MH_M : failed to allocate storage" );
      return false;
    }
    
    // Apply encoding operator
    if( this->encoding_operator_->mult_MH_M( r, &tmp_MH_M ) < 0 ) {
      this->solver_error( "Error: lwSolver::apply_shaping_matrix_mult_MH_M : failed to apply encoding operator (1)" );
      return false;
    }
    
    // Accumulate for overall result
    if( !solver_axpy(get_alpha()*w*this->encoding_operator_->get_weight(), &tmp_MH_M, acc )) {
      this->solver_error( "Error: lwSolver::apply_shaping_matrix_mult_MH_M : failed to accumulate (1)" );
      return false;
    }

    // Accumulate for intermediate (MH_M)^i
    tmp_acc = tmp_MH_M;
    if( !solver_scale(get_alpha()*this->encoding_operator_->get_weight(), &tmp_acc )) {
      this->solver_error( "Error: lwSolver::apply_shaping_matrix_mult_MH_M : failed to accumulate (1)" );
      return false;
    }
    
    // Loop over operators
    for( unsigned int i=0; i<this->regularization_operators_.size(); i++){
      
      // Compute operator mult_MH_M
      if( this->regularization_operators_[i]->mult_MH_M( r, &tmp_MH_M ) < 0 ) {
	this->solver_error( "Error: lwSolver::apply_shaping_matrix_mult_MH_M : failed to apply regularization operator" );
	return false;
      }
      
      // Accumulate
      if( !solver_axpy(get_alpha()*w*this->regularization_operators_[i]->get_weight(), &tmp_MH_M, acc )) {
	this->solver_error( "Error: lwSolver::apply_shaping_matrix_mult_MH_M : failed to accumulate (2)" );
	return false;
      }

      // Accumulate for intermediate (MH_M)^i
      if( !solver_axpy(get_alpha()*this->encoding_operator_->get_weight(), &tmp_MH_M, &tmp_acc )) {
	this->solver_error( "Error: lwSolver::apply_shaping_matrix_mult_MH_M : failed to accumulate (1)" );
	return false;
      }
    }
    
    // Update r
    *r = tmp_acc;

    return true;
  }
  
protected:
  
  // Maximum number of iterations
  unsigned int iterations_;
  REAL alpha_;
};
