#pragma once

#include "solver.h"
#include "matrixOperator.h"
#include "vector_td_utilities.h"

#include <vector>
#include <iostream>

template <class REAL, class ELEMENT_TYPE, class ARRAY_TYPE_REAL, class ARRAY_TYPE_ELEMENT> 
class sbSolver : public solver<ARRAY_TYPE_ELEMENT>
{
public:

  sbSolver( int output_mode = solver<ARRAY_TYPE_ELEMENT>::OUTPUT_SILENT ) : solver<ARRAY_TYPE_ELEMENT>( output_mode ) { 
    outer_iterations_ = 10;
    inner_iterations_ = 10;
  }
  
  virtual ~sbSolver() {}

  virtual int set_solver( boost::shared_ptr< solver<ARRAY_TYPE_ELEMENT> > solver ) {
    inner_solver_ = solver;
    return 0;
  }

  virtual int set_encoding_operator( boost::shared_ptr< matrixOperator<REAL, ARRAY_TYPE_ELEMENT> > op ) {
    encoding_operator_ = op;
    return 0;
  }
  
  virtual int add_regularization_operator( boost::shared_ptr< matrixOperator<REAL, ARRAY_TYPE_ELEMENT> > op ) {
    regularization_operators_.push_back(op);
    return 0;
  }

  virtual int add_regularization_group_operator( boost::shared_ptr< matrixOperator<REAL, ARRAY_TYPE_ELEMENT> > op ) {
    regularization_group_operators_.push_back(op);
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

  virtual bool solver_clear( ARRAY_TYPE_ELEMENT* ) = 0;
  virtual bool solver_clear( ARRAY_TYPE_REAL* ) = 0;
  virtual bool solver_sqrt( ARRAY_TYPE_REAL* ) = 0;
  virtual bool solver_scal( ELEMENT_TYPE, ARRAY_TYPE_ELEMENT* ) = 0;
  virtual bool solver_axpy( ELEMENT_TYPE, ARRAY_TYPE_ELEMENT*, ARRAY_TYPE_ELEMENT* ) = 0;
  virtual bool solver_axpy( REAL, ARRAY_TYPE_REAL*, ARRAY_TYPE_REAL* ) = 0;
  virtual REAL solver_asum( ARRAY_TYPE_ELEMENT* ) = 0;
  virtual boost::shared_ptr<ARRAY_TYPE_REAL> solver_norm( ARRAY_TYPE_ELEMENT* ) = 0;
  virtual boost::shared_ptr<ARRAY_TYPE_REAL> solver_norm_squared( ARRAY_TYPE_ELEMENT* ) = 0;
  virtual bool solver_shrink1( REAL, ARRAY_TYPE_ELEMENT*, ARRAY_TYPE_ELEMENT* ) = 0;
  virtual bool solver_shrinkd( REAL, ARRAY_TYPE_REAL*, ARRAY_TYPE_ELEMENT*, ARRAY_TYPE_ELEMENT* ) = 0;

  virtual boost::shared_ptr<ARRAY_TYPE_ELEMENT> solve( ARRAY_TYPE_ELEMENT *_f )
  {

    // Some tests to see if we are ready to go...
    //

    if( !inner_solver_.get() ){
      this->solver_error( "sbSolver::solve : inner solver has not been set" );
      return boost::shared_ptr<ARRAY_TYPE_ELEMENT>();
    }
    
    if( !encoding_operator_.get() ){
      this->solver_error( "sbSolver::solve : encoding operator has not been set" );
      return boost::shared_ptr<ARRAY_TYPE_ELEMENT>();
    }

    if( regularization_operators_.size() == 0 && regularization_group_operators_.size() == 0 ){
      this->solver_error( "sbSolver::solve : at least one matrix regularizer must be added" );
      return boost::shared_ptr<ARRAY_TYPE_ELEMENT>();
    }
    
    if( !image_dims_.get() ){
      this->solver_error( "sbSolver::solve : image dimensions have not been set" );
      return boost::shared_ptr<ARRAY_TYPE_ELEMENT>();
    }

    // Make a copy of f - we are going to "normalize" the input
    //
    ARRAY_TYPE_ELEMENT f(*_f);
    if( !f.get_data_ptr() ){
      this->solver_error( "sbSolver::solve : memory allocation of f failed" );
      return boost::shared_ptr<ARRAY_TYPE_ELEMENT>();
    }

    // Initialize u_k to E^H f 
    //
    boost::shared_ptr<ARRAY_TYPE_ELEMENT> u_k = boost::shared_ptr<ARRAY_TYPE_ELEMENT>(new ARRAY_TYPE_ELEMENT());
    if( u_k.get() ) u_k->create( image_dims_.get() );

    if( !u_k.get() || !u_k->get_data_ptr() ){
      this->solver_error( "sbSolver::solve : memory allocation of u_k failed" );
      return boost::shared_ptr<ARRAY_TYPE_ELEMENT>();
    }    

    if( encoding_operator_->mult_MH( &f, u_k.get() ) < 0 ){
      this->solver_error( "sbSolver::solve : adjoint encoding operation failed on f" );
      return boost::shared_ptr<ARRAY_TYPE_ELEMENT>();
    }

    // Normalize u_k and f
    //
    ELEMENT_TYPE image_scale;
    {
      // Normalize to an average energy of "one intensity unit per image element"
      REAL sum = solver_asum( u_k.get() );
      image_scale = mul<REAL>(( (REAL) (u_k->get_number_of_elements())/sum), get_one<ELEMENT_TYPE>() );
      solver_scal( image_scale, u_k.get() );
      solver_scal( image_scale, &f );
      std::cout << std::endl << "SCALE: " << norm<REAL>(image_scale);
    }

    // Initialize f_k to f
    //
    ARRAY_TYPE_ELEMENT f_k(f);
    if( !f_k.get_data_ptr() ){
      this->solver_error( "sbSolver::solve : memory allocation of f_k failed" );
      return boost::shared_ptr<ARRAY_TYPE_ELEMENT>();
    }

    // Initialize d_k and b_k arrays to zero images
    //

    unsigned int num_reg_operators = regularization_operators_.size() + regularization_group_operators_.size();

    boost::shared_array< boost::shared_ptr<ARRAY_TYPE_ELEMENT> > d_k = 
      boost::shared_array< boost::shared_ptr<ARRAY_TYPE_ELEMENT> >( new boost::shared_ptr<ARRAY_TYPE_ELEMENT>[num_reg_operators] );

    boost::shared_array< boost::shared_ptr<ARRAY_TYPE_ELEMENT> > b_k = 
      boost::shared_array< boost::shared_ptr<ARRAY_TYPE_ELEMENT> >( new boost::shared_ptr<ARRAY_TYPE_ELEMENT>[num_reg_operators] );
    
    if( !d_k.get() || !b_k.get() ){
      this->solver_error( "sbSolver::solve : memory allocation of d_k or b_k failed (1)" );
      return boost::shared_ptr<ARRAY_TYPE_ELEMENT>();
    }

    for( unsigned int i=0; i<num_reg_operators; i++ ){
      
      d_k[i] = boost::shared_ptr<ARRAY_TYPE_ELEMENT>(new ARRAY_TYPE_ELEMENT());

      if( d_k[i] ) d_k[i]->create( image_dims_.get() );

      if( !d_k[i]->get_data_ptr() ){
	this->solver_error( "sbSolver::solve : memory allocation of d_k failed" );
	return boost::shared_ptr<ARRAY_TYPE_ELEMENT>();
      }

      if( !solver_clear( d_k[i].get() )){
	this->solver_error( "sbSolver::solve : failed to clear internal memory buffer (1)" );
	return boost::shared_ptr<ARRAY_TYPE_ELEMENT>();
      }
      
      b_k[i] = boost::shared_ptr<ARRAY_TYPE_ELEMENT>(new ARRAY_TYPE_ELEMENT());
      
      if( b_k[i] ) b_k[i]->create( image_dims_.get() );
      
      if( !b_k[i]->get_data_ptr() ){
	this->solver_error( "sbSolver::solve : memory allocation of b_k failed" );
	return boost::shared_ptr<ARRAY_TYPE_ELEMENT>();
      }
      
      if( !solver_clear( b_k[i].get() )){
	this->solver_error( "sbSolver::solve : failed to clear internal memory buffer (2)" );
	return boost::shared_ptr<ARRAY_TYPE_ELEMENT>();
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
	ARRAY_TYPE_ELEMENT rhs;
	if( rhs.create( image_dims_.get() ) < 0 ){
	  this->solver_error( "sbSolver::solve : memory allocation of rhs failed" );
	  return boost::shared_ptr<ARRAY_TYPE_ELEMENT>();
	}    
	
	// Add encoding operator
	if( encoding_operator_->mult_MH( &f_k, &rhs ) < 0 ){
	  this->solver_error( "sbSolver::solve : adjoint encoding operation failed on rhs" );
	  return boost::shared_ptr<ARRAY_TYPE_ELEMENT>();
	}
    	
	unsigned int operator_idx = 0;
	
	// Add regularization operators
	for( unsigned int i=0; i<regularization_operators_.size(); i++ ){
	  
	  ARRAY_TYPE_ELEMENT tmp_diff, reg_out;
	  if( tmp_diff.create( image_dims_.get() ) < 0 || reg_out.create( image_dims_.get() ) < 0 ){
	    this->solver_error( "sbSolver::solve : memory allocation for regularization operator failed in rhs computation" );
	    return boost::shared_ptr<ARRAY_TYPE_ELEMENT>();
	  }    
	  
	  tmp_diff = *d_k[operator_idx];

	  if( !solver_axpy( get_zero<ELEMENT_TYPE>()-get_one<ELEMENT_TYPE>(), b_k[operator_idx].get(), &tmp_diff )){
	    this->solver_error( "sbSolver::solve : computation of regularization argument failed in rhs computation" );
	    return boost::shared_ptr<ARRAY_TYPE_ELEMENT>();
	  }    
	  
	  if( regularization_operators_.at(i)->mult_MH( &tmp_diff, &reg_out ) < 0 ){
	    this->solver_error( "sbSolver::solve : application of regularization operator failed in rhs computation" );
	    return boost::shared_ptr<ARRAY_TYPE_ELEMENT>();
	  }    
	  
	  if( !solver_axpy( mul<REAL>(regularization_operators_.at(i)->get_weight(), get_one<ELEMENT_TYPE>()), &reg_out, &rhs )){
	    this->solver_error( "sbSolver::solve : accumulation in rhs computation failed" );
	    return boost::shared_ptr<ARRAY_TYPE_ELEMENT>();
	  }
	  operator_idx++;
	}	

	// Add regularization groups
	for( unsigned int i=0; i<regularization_group_operators_.size(); i++ ){
	  
	  ARRAY_TYPE_ELEMENT tmp_diff, reg_out;
	  if( tmp_diff.create( image_dims_.get() ) < 0 || reg_out.create( image_dims_.get() ) < 0 ){
	    this->solver_error( "sbSolver::solve : memory allocation for regularization operator failed in rhs computation" );
	    return boost::shared_ptr<ARRAY_TYPE_ELEMENT>();
	  }    
	  
	  tmp_diff = *d_k[operator_idx];
	  
	  if( !solver_axpy( get_zero<ELEMENT_TYPE>()-get_one<ELEMENT_TYPE>(), b_k[operator_idx].get(), &tmp_diff )){
	    this->solver_error( "sbSolver::solve : computation of regularization argument failed in rhs computation" );
	    return boost::shared_ptr<ARRAY_TYPE_ELEMENT>();
	  }    
	  
	  if( regularization_group_operators_.at(i)->mult_MH( &tmp_diff, &reg_out ) < 0 ){
	    this->solver_error( "sbSolver::solve : application of regularization operator failed in rhs computation" );
	    return boost::shared_ptr<ARRAY_TYPE_ELEMENT>();
	  }    
	  
	  if( !solver_axpy( mul<REAL>(regularization_group_operators_.at(i)->get_weight(), get_one<ELEMENT_TYPE>()), &reg_out, &rhs )){
	    this->solver_error( "sbSolver::solve : accumulation in rhs computation failed" );
	    return boost::shared_ptr<ARRAY_TYPE_ELEMENT>();
	  }
	  operator_idx++;
	}	

	// Solve for u_k
	//
	{

	  boost::shared_ptr<ARRAY_TYPE_ELEMENT> tmp = inner_solver_->solve(&rhs);

	  if( outer_iteration == outer_iterations_-1 ) 
	    break;

	  // Output change in u_k	
	  if( !solver_axpy( get_zero<ELEMENT_TYPE>()-get_one<ELEMENT_TYPE>(), tmp.get(), u_k.get() )){
	    this->solver_error( "sbSolver::solve : error computing u_k delta" );
	    return boost::shared_ptr<ARRAY_TYPE_ELEMENT>();
	  }	  
	  std::cout << std::endl << "u_k delta: " << solver_asum(u_k.get()) << std::endl;      	
	  
	  // Update u_k
	  *u_k = *(tmp.get());
	}

	operator_idx = 0;
	
	// Update d_k and b_k
	//
	for( unsigned int i=0; i<regularization_operators_.size(); i++ ){ // For every operator
	  
	  ARRAY_TYPE_ELEMENT tmp_sum, reg_out;
	  if( tmp_sum.create( image_dims_.get() ) < 0 || reg_out.create( image_dims_.get() ) < 0 ){
	    this->solver_error( "sbSolver::solve : memory allocation for regularization operator failed in {d_k,b_k} update" );
	    return boost::shared_ptr<ARRAY_TYPE_ELEMENT>();
	  }
	  
	  tmp_sum = *b_k[operator_idx];
	  
	  if( regularization_operators_.at(i)->mult_M( u_k.get(), &reg_out ) < 0 ){
	    this->solver_error( "sbSolver::solve : application of regularization operator failed in {d_k,b_k} update" );
	    return boost::shared_ptr<ARRAY_TYPE_ELEMENT>();
	  }
	  
	  if( !solver_axpy( get_one<ELEMENT_TYPE>(), &reg_out, &tmp_sum )){
	    this->solver_error( "sbSolver::solve : computation of shrinkage argument for d_k failed" );
	    return boost::shared_ptr<ARRAY_TYPE_ELEMENT>();
	  }
	  
	  // Update of d_k
	  if( !solver_shrink1( reciprocal<REAL>(regularization_operators_.at(i)->get_weight()), &tmp_sum, d_k[operator_idx].get() )){
	    this->solver_error( "sbSolver::solve : shrinkage of d_k failed" );
	    return boost::shared_ptr<ARRAY_TYPE_ELEMENT>();
	  }
	  
	  if( !solver_axpy( get_zero<ELEMENT_TYPE>()-get_one<ELEMENT_TYPE>(), d_k[operator_idx].get(), &reg_out )){
	    this->solver_error( "sbSolver::solve : computation of update argument to b_k failed" );
	    return boost::shared_ptr<ARRAY_TYPE_ELEMENT>();
	  }
	  
	  // Update of b_k
	  if( !solver_axpy( get_one<ELEMENT_TYPE>(), &reg_out, b_k[operator_idx].get() )){
	    this->solver_error( "sbSolver::solve : update of b_k failed" );
	    return boost::shared_ptr<ARRAY_TYPE_ELEMENT>();
	  }
	  operator_idx++;
	}
	
	unsigned int k = regularization_group_operators_.size();
	ARRAY_TYPE_ELEMENT sums[k], reg_out[k];
	
	ARRAY_TYPE_REAL s_k;
	if( s_k.create(image_dims_.get()) < 0 ){
	  this->solver_error( "sbSolver::solve : memory allocation for s_k failed" );
	  return boost::shared_ptr<ARRAY_TYPE_ELEMENT>();
	}
	
	if( !solver_clear(&s_k) ){
	  this->solver_error( "sbSolver::solve : failed to clear s_k" );
	  return boost::shared_ptr<ARRAY_TYPE_ELEMENT>();
	} 	  
	
	for( unsigned int j=0; j<k; j++ ){
	  
	  if( sums[j].create( image_dims_.get() ) < 0 || reg_out[j].create( image_dims_.get() ) < 0 ){
	    this->solver_error( "sbSolver::solve : memory allocation for regularization operator failed in {d_k,b_k} update" );
	    return boost::shared_ptr<ARRAY_TYPE_ELEMENT>();
	  }
	  
	  sums[j] = *b_k[operator_idx+j];
	  
	  if( regularization_group_operators_.at(j)->mult_M( u_k.get(), &reg_out[j] ) < 0 ){
	    this->solver_error( "sbSolver::solve : application of regularization operator failed in {d_k,b_k} update" );
	    return boost::shared_ptr<ARRAY_TYPE_ELEMENT>();
	  }
	  
	  if( !solver_axpy( get_one<ELEMENT_TYPE>(), &reg_out[j], &sums[j] )){
	    this->solver_error( "sbSolver::solve : computation of shrinkage argument for d_k failed" );
	    return boost::shared_ptr<ARRAY_TYPE_ELEMENT>();
	  }
	  
	  boost::shared_ptr<ARRAY_TYPE_REAL> tmp_s_k = solver_norm_squared(&sums[j]);
	  if( !solver_axpy( get_one<REAL>(), tmp_s_k.get(), &s_k )){
	    this->solver_error( "sbSolver::solve : accumulation of s_k failed" );
	    return boost::shared_ptr<ARRAY_TYPE_ELEMENT>();
	  }
	}
	
	if( !solver_sqrt(&s_k) ){
	  this->solver_error( "sbSolver::solve : sqrt of s_k failed" );
	  return boost::shared_ptr<ARRAY_TYPE_ELEMENT>();
	}
	
	for( unsigned int j=0; j<k; j++ ){
	  
	  // Update of d_k
	  if( !solver_shrinkd( reciprocal<REAL>(regularization_group_operators_.at(j)->get_weight()), &s_k, &sums[j], d_k[operator_idx+j].get() )){
	    this->solver_error( "sbSolver::solve : shrinkage of d_k failed" );
	    return boost::shared_ptr<ARRAY_TYPE_ELEMENT>();
	  }
	  
	  if( !solver_axpy( get_zero<ELEMENT_TYPE>()-get_one<ELEMENT_TYPE>(), d_k[operator_idx+j].get(), &reg_out[j] )){
	    this->solver_error( "sbSolver::solve : computation of update argument to b_k failed" );
	    return boost::shared_ptr<ARRAY_TYPE_ELEMENT>();
	  }
	  
	  // Update of b_k
	  if( !solver_axpy( get_one<ELEMENT_TYPE>(), &reg_out[j], b_k[operator_idx+j].get() )){
	    this->solver_error( "sbSolver::solve : update of b_k failed" );
	    return boost::shared_ptr<ARRAY_TYPE_ELEMENT>();
	  }
	}
	operator_idx += k;
      } // end of inner loop
      
      // Update f_k
      //
      ARRAY_TYPE_ELEMENT encoded_image;
      if( encoded_image.create( f.get_dimensions().get() ) < 0 ){
	this->solver_error( "sbSolver::solve : memory allocation for encoded image failed" );
	return boost::shared_ptr<ARRAY_TYPE_ELEMENT>();
      }
      
      if( encoding_operator_->mult_M( u_k.get(), &encoded_image ) < 0 ){
	this->solver_error( "sbSolver::solve : computation of encoded image failed" );
	return boost::shared_ptr<ARRAY_TYPE_ELEMENT>();
      }

      if( !solver_axpy( get_zero<ELEMENT_TYPE>()-get_one<ELEMENT_TYPE>(), &f, &encoded_image )){ // notice the deliberate sign manipulation
	this->solver_error( "sbSolver::solve : computation of update argument to f_k failed" );
	return boost::shared_ptr<ARRAY_TYPE_ELEMENT>();
      }

      // Output residual
      std::cout << std::endl << "Residual: " << cuNDA_asum<REAL>(&encoded_image) << std::endl;
      std::cout << "---------"  << std::endl;
      
      // and update f_k with residual
      if( !solver_axpy( get_zero<ELEMENT_TYPE>()-get_one<ELEMENT_TYPE>(), &encoded_image, &f_k )){ // notice the deliberate sign manipulation
	this->solver_error( "sbSolver::solve : computation of update argument to f_k failed" );
	return boost::shared_ptr<ARRAY_TYPE_ELEMENT>();
      }
    } // end of outer loop
        
    // Undo intermediate scaling of u_k ...
    solver_scal( reciprocal(image_scale), u_k.get() );

    // ... and return the result
    //
    
    return u_k;
  }

protected:
  unsigned int outer_iterations_, inner_iterations_;
  boost::shared_ptr< std::vector<unsigned int> > image_dims_;
  boost::shared_ptr< solver<ARRAY_TYPE_ELEMENT> > inner_solver_;
  boost::shared_ptr< matrixOperator<REAL, ARRAY_TYPE_ELEMENT> > encoding_operator_;
  std::vector< boost::shared_ptr< matrixOperator<REAL, ARRAY_TYPE_ELEMENT> > > regularization_operators_;
  std::vector< boost::shared_ptr< matrixOperator<REAL, ARRAY_TYPE_ELEMENT> > > regularization_group_operators_;
};
