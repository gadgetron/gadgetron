#pragma once

#include "solver.h"
#include "matrixOperator.h"
#include "vector_td_utilities.h"

#include <vector>
#include <iostream>

template <class REAL, class ELEMENT_TYPE, class ARRAY_TYPE> class sbSolver : public solver<ARRAY_TYPE>
{
public:

  sbSolver( int output_mode = solver<ARRAY_TYPE>::OUTPUT_SILENT ) : solver<ARRAY_TYPE>( output_mode ) { 
    outer_iterations_ = 10;
    inner_iterations_ = 10;
    regularization_operators_ = boost::shared_ptr< std::vector< boost::shared_ptr< matrixOperator<REAL, ARRAY_TYPE> > > > 
      ( new std::vector< boost::shared_ptr< matrixOperator<REAL, ARRAY_TYPE> > >() );
  }
  
  virtual ~sbSolver() {}

  virtual int set_solver( boost::shared_ptr< solver<ARRAY_TYPE> > solver ) {
    inner_solver_ = solver;
    return 0;
  }

  virtual int set_encoding_operator( boost::shared_ptr< matrixOperator<REAL, ARRAY_TYPE> > op ) {
    encoding_operator_ = op;
    return 0;
  }
  
  virtual int add_regularization_operator( boost::shared_ptr< matrixOperator<REAL, ARRAY_TYPE> > op ) {
    regularization_operators_->push_back(op);
    return 0;
  }
  
  virtual void set_outer_iterations( unsigned int iterations ) {
    outer_iterations_ = iterations;
  }

  virtual void set_inner_iterations( unsigned int iterations ) {
    inner_iterations_ = iterations;
  }

  virtual void set_image_dimensions( boost::shared_ptr< std::vector<unsigned int> > dims )
  {
    image_dims_ = dims;
  }

  virtual bool solver_clear( ARRAY_TYPE* ) = 0;
  virtual bool solver_scal( ELEMENT_TYPE, ARRAY_TYPE* ) = 0;
  virtual bool solver_axpy( ELEMENT_TYPE, ARRAY_TYPE*, ARRAY_TYPE* ) = 0;
  virtual bool solver_shrink( REAL, ARRAY_TYPE*, ARRAY_TYPE* ) = 0;

  virtual boost::shared_ptr<ARRAY_TYPE> solve( ARRAY_TYPE *f )
  {

    // Some tests to see if we are ready to go...
    //

    if( !inner_solver_.get() ){
      this->solver_error( "sbSolver::solve : inner solver has not been set" );
      return boost::shared_ptr<ARRAY_TYPE>();
    }
    
    if( !encoding_operator_.get() ){
      this->solver_error( "sbSolver::solve : encoding operator has not been set" );
      return boost::shared_ptr<ARRAY_TYPE>();
    }

    if( regularization_operators_->size() == 0 ){
      this->solver_error( "sbSolver::solve : at least one matrix regularizer must be added" );
      return boost::shared_ptr<ARRAY_TYPE>();
    }
    
    if( !image_dims_.get() ){
      this->solver_error( "sbSolver::solve : image dimensions have not been set" );
      return boost::shared_ptr<ARRAY_TYPE>();
    }

    // Initialize f_k to f
    //
    ARRAY_TYPE f_k(*f);
    if( !f_k.get_data_ptr() ){
      this->solver_error( "sbSolver::solve : memory allocation of f_k failed" );
      return boost::shared_ptr<ARRAY_TYPE>();
    }
    
    // Initialize u_k to E^H f 
    //
    ARRAY_TYPE *u_k = new ARRAY_TYPE();
    if( u_k ) u_k->create( image_dims_.get() );

    if( !u_k || !u_k->get_data_ptr() ){
      this->solver_error( "sbSolver::solve : memory allocation of u_k failed" );
      return boost::shared_ptr<ARRAY_TYPE>();
    }    

    if( encoding_operator_->mult_MH( f, u_k ) < 0 ){
      this->solver_error( "sbSolver::solve : adjoint encoding operation failed on f" );
      return boost::shared_ptr<ARRAY_TYPE>();
    }

    // Initialize d_k and b_k arrays to zero images
    //
    ARRAY_TYPE **d_k = new ARRAY_TYPE*[regularization_operators_->size()];
    ARRAY_TYPE **b_k = new ARRAY_TYPE*[regularization_operators_->size()];
    
    if( !d_k || !b_k ){
      this->solver_error( "sbSolver::solve : memory allocation of d_k or b_k failed (1)" );
      return boost::shared_ptr<ARRAY_TYPE>();
    }

    for( unsigned int i=0; i<regularization_operators_->size(); i++ ){
      
      d_k[i] = new ARRAY_TYPE();

      if( d_k[i] ) d_k[i]->create( image_dims_.get() );

      if( !d_k[i]->get_data_ptr() ){
	this->solver_error( "sbSolver::solve : memory allocation of d_k failed" );
	return boost::shared_ptr<ARRAY_TYPE>();
      }

      if( !solver_clear( d_k[i] )){
	this->solver_error( "sbSolver::solve : failed to clear internal memory buffer (1)" );
	return boost::shared_ptr<ARRAY_TYPE>();
      }
      
      b_k[i] = new ARRAY_TYPE();
      
      if( b_k[i] ) b_k[i]->create( image_dims_.get() );
      
      if( !b_k[i]->get_data_ptr() ){
	this->solver_error( "sbSolver::solve : memory allocation of b_k failed" );
	return boost::shared_ptr<ARRAY_TYPE>();
      }
      
      if( !solver_clear( d_k[i] )){
	this->solver_error( "sbSolver::solve : failed to clear internal memory buffer (2)" );
	return boost::shared_ptr<ARRAY_TYPE>();
      } 
    }
    
    // Outer loop
    //
    for( unsigned int outer_iteration=0; outer_iteration<outer_iterations_; outer_iteration++ ) {
      
      std::cout << std::endl << "SB outer loop iteration " << outer_iteration << std::endl;
      std::cout << "------------------------ " << std::endl << std::endl;

      // Inner loop
      //
      
      for( unsigned int inner_iteration=0; inner_iteration<inner_iterations_; inner_iteration++ ) {
	
	std::cout << std::endl << "SB inner loop iteration " << inner_iteration << std::endl;
	std::cout << "------------------------ " << std::endl << std::endl;

	// Form rhs for inner loop solver
	// 
	ARRAY_TYPE rhs;
	if( rhs.create( image_dims_.get() ) < 0 ){
	  this->solver_error( "sbSolver::solve : memory allocation of rhs failed" );
	  return boost::shared_ptr<ARRAY_TYPE>();
	}    
	
	if( encoding_operator_->mult_MH( &f_k, &rhs ) < 0 ){
	  this->solver_error( "sbSolver::solve : adjoint encoding operation failed on rhs" );
	  return boost::shared_ptr<ARRAY_TYPE>();
	}
    
	if( !solver_scal( mul<REAL>(encoding_operator_->get_weight(), get_one<ELEMENT_TYPE>()), &rhs ) ){
	  this->solver_error( "sbSolver::solve : could not apply weight to rhs" );
	  return boost::shared_ptr<ARRAY_TYPE>();
	}
	
	for( unsigned int i=0; i<regularization_operators_->size(); i++ ){
	  
	  ARRAY_TYPE tmp_diff, reg_out;
	  if( tmp_diff.create( image_dims_.get() ) < 0 || reg_out.create( image_dims_.get() ) < 0 ){
	    this->solver_error( "sbSolver::solve : memory allocation for regularization operator failed in rhs computation" );
	    return boost::shared_ptr<ARRAY_TYPE>();
	  }    
	  
	  tmp_diff = *d_k[i];

	  if( !solver_axpy( get_zero<ELEMENT_TYPE>()-get_one<ELEMENT_TYPE>(), b_k[i], &tmp_diff )){
	    this->solver_error( "sbSolver::solve : computation of regularization argument failed in rhs computation" );
	    return boost::shared_ptr<ARRAY_TYPE>();
	  }    
	  
	  if( (regularization_operators_->at(i))->mult_MH( &tmp_diff, &reg_out ) < 0 ){
	    this->solver_error( "sbSolver::solve : application of regularization operator failed in rhs computation" );
	    return boost::shared_ptr<ARRAY_TYPE>();
	  }    
	  
	  if( !solver_axpy( mul<REAL>(regularization_operators_->at(i)->get_weight(), get_one<ELEMENT_TYPE>()), &reg_out, &rhs )){
	    this->solver_error( "sbSolver::solve : accumulation in rhs computation failed" );
	    return boost::shared_ptr<ARRAY_TYPE>();
	  }    	  	  	  
	}
	
	// Solve for u_k
	//
	{
	  boost::shared_ptr<ARRAY_TYPE> tmp = inner_solver_->solve(&rhs);
	  
	  // Output change in u_k	
	  if( !solver_axpy( get_zero<ELEMENT_TYPE>()-get_one<ELEMENT_TYPE>(), tmp.get(), u_k )){
	    this->solver_error( "sbSolver::solve : error computing u_k delta" );
	    return boost::shared_ptr<ARRAY_TYPE>();
	  }	  
	  std::cout << std::endl << "u_k delta: " << cuNDA_asum<REAL>(u_k) << std::endl;      	
	  
	  // Update u_k
	  *u_k = *(tmp.get());
	}
	
	// Update d_k and b_k (using "1D shrinkage" for now)
	//
	for( unsigned int i=0; i<regularization_operators_->size(); i++ ){
	  
	  ARRAY_TYPE tmp_sum, reg_out;
	  if( tmp_sum.create( image_dims_.get() ) < 0 || reg_out.create( image_dims_.get() ) < 0 ){
	    this->solver_error( "sbSolver::solve : memory allocation for regularization operator failed in {d_k,b_k} update" );
	    return boost::shared_ptr<ARRAY_TYPE>();
	  }
	  
	  tmp_sum = *b_k[i];
	  
	  if( regularization_operators_->at(i)->mult_M( u_k, &reg_out ) < 0 ){
	    this->solver_error( "sbSolver::solve : application of regularization operator failed in {d_k,b_k} update" );
	    return boost::shared_ptr<ARRAY_TYPE>();
	  }
	  
	  if( !solver_axpy( get_one<ELEMENT_TYPE>(), &reg_out, &tmp_sum )){
	    this->solver_error( "sbSolver::solve : computation of shrinkage argument for d_k failed" );
	    return boost::shared_ptr<ARRAY_TYPE>();
	  }
	  
	  if( !solver_shrink( reciprocal<REAL>(regularization_operators_->at(i)->get_weight()), &tmp_sum, d_k[i] )){
	    this->solver_error( "sbSolver::solve : shrinkage of d_k failed" );
	    return boost::shared_ptr<ARRAY_TYPE>();
	  }
	  
	  // Update of d_k
	  if( !solver_axpy( get_zero<ELEMENT_TYPE>()-get_one<ELEMENT_TYPE>(), d_k[i], &reg_out )){
	    this->solver_error( "sbSolver::solve : computation of update argument to b_k failed" );
	    return boost::shared_ptr<ARRAY_TYPE>();
	  }

	  // Update of b_k
	  if( !solver_axpy( get_one<ELEMENT_TYPE>(), &reg_out, b_k[i] )){
	    this->solver_error( "sbSolver::solve : update of b_k failed" );
	    return boost::shared_ptr<ARRAY_TYPE>();
	  }
	}
      } // end of inner loop

      // Update f_k
      //
      ARRAY_TYPE encoded_image;
      if( encoded_image.create( f->get_dimensions().get() ) < 0 ){
	this->solver_error( "sbSolver::solve : memory allocation for encoded image failed" );
	return boost::shared_ptr<ARRAY_TYPE>();
      }
      
      if( encoding_operator_->mult_M( u_k, &encoded_image ) < 0 ){
	this->solver_error( "sbSolver::solve : computation of encoded image failed" );
	return boost::shared_ptr<ARRAY_TYPE>();
      }

      if( !solver_axpy( get_zero<ELEMENT_TYPE>()-get_one<ELEMENT_TYPE>(), f, &encoded_image )){ // notice the deliberate sign manipulation
	this->solver_error( "sbSolver::solve : computation of update argument to f_k failed" );
	return boost::shared_ptr<ARRAY_TYPE>();
      }

      // Output residual
      std::cout << std::endl << "Residual: " << cuNDA_asum<REAL>(&encoded_image) << std::endl;
      std::cout << "---------"  << std::endl;
      
      // and update f_k with residual
      if( !solver_axpy( get_zero<ELEMENT_TYPE>()-get_one<ELEMENT_TYPE>(), &encoded_image, &f_k )){ // notice the deliberate sign manipulation
	this->solver_error( "sbSolver::solve : computation of update argument to f_k failed" );
	return boost::shared_ptr<ARRAY_TYPE>();
      }      
    } // end of outer loop
    
    // Clean up
    //

    for( unsigned int i=0; i<regularization_operators_->size(); i++ ){      
      delete d_k[i];
      delete b_k[i];
    }

    delete[] d_k;
    delete[] b_k;

    // and return the result
    //

    return boost::shared_ptr<ARRAY_TYPE>(u_k);
  }

protected:
  unsigned int outer_iterations_, inner_iterations_;
  boost::shared_ptr< std::vector<unsigned int> > image_dims_;
  boost::shared_ptr< solver<ARRAY_TYPE> > inner_solver_;
  boost::shared_ptr< matrixOperator<REAL, ARRAY_TYPE> > encoding_operator_;
  boost::shared_ptr< std::vector< boost::shared_ptr< matrixOperator<REAL, ARRAY_TYPE> > > > regularization_operators_;
};
