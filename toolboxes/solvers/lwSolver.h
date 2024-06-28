/*
  An implementation of the "Generalized Landweber Solver" based on the paper
  "Theory and methods related to the singular-function expansion and Landweber's iteration..."
  by O.N. Strand, Siam J Numer, Anal. 1974;11(4):798-825.
*/

#pragma once

#include "linearOperatorSolver.h"

#include <vector>
#include <iostream>

namespace Gadgetron{

  template <class ARRAY_TYPE> class lwSolver
    : public linearOperatorSolver< ARRAY_TYPE>
  {

  protected:
    typedef typename ARRAY_TYPE::element_type ELEMENT_TYPE;
    typedef typename realType<ELEMENT_TYPE>::type REAL;

  public:

    // Constructor
    lwSolver() : linearOperatorSolver<ARRAY_TYPE>() {
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

    // Inherited solver interface
    virtual boost::shared_ptr<ARRAY_TYPE> solve( ARRAY_TYPE *b )
    {   
      // Initial validity checks
      //

      std::vector<unsigned int> image_dims = this->encoding_operator_->get_domain_dimensions();

      if( image_dims.size() == 0 ){
	throw std::runtime_error("Error: lwSolver::solve : domain dimensions not set on encoding operator" );
      }
        
      // Allocate solution array.
      // Clear or set to x0 if provided
      //

      boost::shared_ptr<ARRAY_TYPE> x( new ARRAY_TYPE() );
      if( this->get_x0().get() ){
	*x = *(this->get_x0());
      }
      else{
    	x->create( &image_dims );
    	x->clear();
      }    

      ARRAY_TYPE x_prev;

      // Main solver iteration loop
      //
    
      for( unsigned int iteration=0; iteration<iterations_; iteration++ ){
      
	// Keep previous x for convergence reporting
	// 

	if( this->output_mode_ >= solver<ARRAY_TYPE, ARRAY_TYPE>::OUTPUT_VERBOSE ){
	  x_prev = *x;

	}
      
	// Compute residual image, i.e. A^T(b-Ax_k)
	//
      
	boost::shared_ptr<ARRAY_TYPE> r = compute_residual_image( x.get(), b );

	// Multiply residual with shaping matrix
	//

	boost::shared_ptr<ARRAY_TYPE> rr = apply_shaping_matrix( r.get() );

	// Update x
	//
	axpy( get_alpha(), rr.get(), x.get() );
      
	if( this->output_mode_ >= solver<ARRAY_TYPE, ARRAY_TYPE>::OUTPUT_VERBOSE ){
	  axpy( ELEMENT_TYPE(-1), x.get(), &x_prev );
	  GDEBUG_STREAM(" iteration: " << iteration << ", delta x: " << solver_asum(&x_prev) << std::endl);
	}      
      }
    
      return x;
    }

  protected:
    virtual boost::shared_ptr<ARRAY_TYPE> compute_residual_image( ARRAY_TYPE *x, ARRAY_TYPE *b )
    {    
      // Allocate some temporary storage and the esult array
      //
    

      boost::shared_ptr<ARRAY_TYPE> res = boost::shared_ptr<ARRAY_TYPE>(new ARRAY_TYPE(x->get_dimensions()));

      ARRAY_TYPE tmp_M(b->get_dimensions());
      ARRAY_TYPE tmp_acc(b->get_dimensions());
        
      // Clear accumulation buffer to b
      tmp_acc = *b;
    
      // Apply encoding operator to current solution
      this->encoding_operator_->mult_M( x, &tmp_M );
    
      // Find residual
      axpy(REAL(-1), &tmp_M, &tmp_acc );
    
      // Adjoint residual    
      this->encoding_operator_->mult_MH( &tmp_acc, res.get());
      // Apply encoding operator weight
      *res *= this->encoding_operator_->get_weight();
    
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
      std::vector<unsigned int> image_dims = this->encoding_operator_->get_domain_dimensions();
      boost::shared_ptr<ARRAY_TYPE> res = boost::shared_ptr<ARRAY_TYPE>(new ARRAY_TYPE(&image_dims ));
    
      // Handle 0th order   
      *res = *r;
      *res *=  REAL(31.5);

      // Handle 1th order
      apply_shape_matrix_mult_MH_M( r, res.get(), REAL(-315) );
    
      // Handle 2th order
      apply_shape_matrix_mult_MH_M( r, res.get(), REAL(1443.75) );

      // Handle 3th order
      apply_shape_matrix_mult_MH_M( r, res.get(), REAL(-3465) );
    
      // Handle 4th order
      apply_shape_matrix_mult_MH_M( r, res.get(), REAL(4504.5) );
    
      // Handle 5th order
      apply_shape_matrix_mult_MH_M( r, res.get(), REAL(-3003) );
    
      // Handle 6th order
      apply_shape_matrix_mult_MH_M( r, res.get(), REAL(804.375) );
    
      // Return result
      return res;
    }

    void apply_shape_matrix_mult_MH_M( ARRAY_TYPE *r, ARRAY_TYPE *acc, REAL w )
    {
      // Temporary storage
      std::vector<unsigned int> image_dims = this->encoding_operator_->get_domain_dimensions();
      ARRAY_TYPE tmp_MH_M(&image_dims), tmp_acc(&image_dims);
    
      // Apply encoding operator
      this->encoding_operator_->mult_MH_M( r, &tmp_MH_M );
    
      // Accumulate for overall result
      axpy(get_alpha()*w*this->encoding_operator_->get_weight(), &tmp_MH_M, acc );

      // Accumulate for intermediate (MH_M)^i
      tmp_acc = tmp_MH_M;
      tmp_acc *= get_alpha()*this->encoding_operator_->get_weight();
    
      // Loop over operators
      for( unsigned int i=0; i<this->regularization_operators_.size(); i++){
      
	// Compute operator mult_MH_M
	this->regularization_operators_[i]->mult_MH_M( r, &tmp_MH_M );
      
	// Accumulate
	axpy(get_alpha()*w*this->regularization_operators_[i]->get_weight(), &tmp_MH_M, acc );

	// Accumulate for intermediate (MH_M)^i
	axpy(get_alpha()*this->encoding_operator_->get_weight(), &tmp_MH_M, &tmp_acc );
      }
    
      // Update r
      *r = tmp_acc;
    }
  
  protected:
  
    // Maximum number of iterations
    unsigned int iterations_;
    REAL alpha_;
  };
}
