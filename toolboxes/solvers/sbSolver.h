/*
  An implementation of the "Generalized Split Bregman Algorithm" - sec. 3.2. of the paper
  "The Split Bregman Method for L1-Regularized Problems" by Tom Goldstein and Stanley Osher. 
  Siam J. Imaging Sciences. Vol. 2, No. 2, pp. 323-343.
*/

#pragma once

#include "linearSolver.h"
#include "vector_td_utilities.h"

#include <vector>
#include <iostream>

template< class REAL, 
	  class ELEMENT_TYPE, 
	  class ARRAY_TYPE_REAL, 
	  class ARRAY_TYPE_ELEMENT, 
	  class INNER_SOLVER,
	  class OPERATOR_CONTAINER > 
class sbSolver : public linearSolver<REAL, ELEMENT_TYPE, ARRAY_TYPE_ELEMENT>
{
public:

  // Constructor
  //

  sbSolver() : linearSolver<REAL, ELEMENT_TYPE, ARRAY_TYPE_ELEMENT>() 
  { 
    tolerance_ = REAL(0);
    outer_iterations_ = 10;
    inner_iterations_ = 1;
    num_reg_operators_ = 0;
    inner_solver_ = boost::shared_ptr<INNER_SOLVER>( new INNER_SOLVER() );
  }
  

  // Destructor
  //

  virtual ~sbSolver() {}
   

  // Add regularization operator (isotropic, multiple operators per group allowed)
  //

  virtual bool add_regularization_group_operator ( boost::shared_ptr< linearOperator<REAL, ARRAY_TYPE_ELEMENT> > op ) 
  {
    if( !op.get() ){
      this->solver_error( "Error: sbSolver::add_regularization_group_operator : NULL operator provided" );
      return false;
    }    

    _regularization_group_operators_.push_back( op );

    return true;
  }
  

  // Add isotroic regularization group (multiple groups allowed)
  //

  virtual bool add_group()
  {
    if( _regularization_group_operators_.size() == 0 ){
      this->solver_error( "Error: sbSolver::add_group : no regularization group operators added" );
      return false;
    }

    regularization_group_operators_.push_back( _regularization_group_operators_ );    
    _regularization_group_operators_.clear();
    
    return true;
  }
  

  // Get regularization group operator
  //

  virtual boost::shared_ptr< linearOperator<REAL, ARRAY_TYPE_ELEMENT> > 
  get_regularization_group_operator( unsigned int group_number, unsigned int idx_in_group )
  {
    if( group_number >= regularization_group_operators_.size() ){
      this->solver_error( "Error: sbSolver::get_regularization_group_operator : group number out of range" );
      return boost::shared_ptr< linearOperator<REAL, ARRAY_TYPE_ELEMENT> >();
    }
    
    if( idx_in_group >= regularization_group_operators_.at(group_number).size() ){
      this->solver_error( "Error: sbSolver::get_regularization_group_operator : idx in group out of range" );
      return boost::shared_ptr< linearOperator<REAL, ARRAY_TYPE_ELEMENT> >();
    }
    
    return regularization_group_operators_.at(group_number).at(idx_in_group);
  }  

  
  // Set/get prior image (PICCS style). 
  // I.e. for every regularization operator (group) R that is added we minimize:
  // alpha|R(x-prior)|_{l1} + (1-alpha)|R(x)|_{l1}
  //

  virtual bool set_prior_image( boost::shared_ptr<ARRAY_TYPE_ELEMENT> prior, REAL alpha )
  {
    prior_ = prior;
    alpha_ = alpha;

    if( alpha_ < REAL(0) ){
      this->solver_warning( "Warning: sbSolver::set_prior_image : alpha cannot be negative. Clamped to 0." );
      alpha_ = REAL(0);
    }
    if( alpha_ > REAL(1) ){
      this->solver_warning( "Warning: sbSolver::set_prior_image : alpha cannot exceed 1. Clamped." );
      alpha_ = REAL(1);
    }
  }
 

  // Get the prior image and corresponding weighing factor
  //
  
  virtual boost::shared_ptr<ARRAY_TYPE_ELEMENT> get_prior_image() { return prior_; }
  virtual REAL get_prior_alpha() { return alpha_; }


  // Set termination criterium tolerance
  //

  virtual void set_tc_tolerance( REAL tolerance ) 
  {
    if( tolerance < REAL(0) ) 
      this->solver_warning( "Warning: sbSolver::set_tc_tolerence : tolerance cannot be negative. Ignored." );
    else tolerance_ = tolerance;
  }


  // Set/get maximum number of outer Split-Bregman iterations
  //

  virtual void set_max_outer_iterations( unsigned int iterations ) { outer_iterations_ = iterations; }
  virtual unsigned int get_max_outer_iterations() { return outer_iterations_; }


  // Set/get maximum number of inner Split-Bregman iterations
  //

  virtual void set_max_inner_iterations( unsigned int iterations ) { inner_iterations_ = iterations; }
  virtual unsigned int get_max_inner_iterations() { return inner_iterations_; }


  // Get the inner solver
  //

  virtual boost::shared_ptr<INNER_SOLVER> get_inner_solver() { return inner_solver_; }

  
  // Core solver functionality to be implemented in a derived class (host/device specific implementations)
  //

  virtual bool solver_clear_real( ARRAY_TYPE_REAL* ) = 0;
  virtual bool solver_clear_element( ARRAY_TYPE_ELEMENT* ) = 0;
  virtual bool solver_sqrt( ARRAY_TYPE_REAL* ) = 0;
  virtual bool solver_scal( ELEMENT_TYPE, ARRAY_TYPE_ELEMENT* ) = 0;
  virtual bool solver_axpy_real( REAL, ARRAY_TYPE_REAL*, ARRAY_TYPE_REAL* ) = 0;
  virtual bool solver_axpy_element( ELEMENT_TYPE, ARRAY_TYPE_ELEMENT*, ARRAY_TYPE_ELEMENT* ) = 0;
  virtual REAL solver_asum( ARRAY_TYPE_ELEMENT* ) = 0;
  virtual boost::shared_ptr<ARRAY_TYPE_REAL> solver_norm( ARRAY_TYPE_ELEMENT* ) = 0;
  virtual bool solver_shrink1( REAL, ARRAY_TYPE_ELEMENT*, ARRAY_TYPE_ELEMENT* ) = 0;
  virtual bool solver_shrinkd( REAL, ARRAY_TYPE_REAL*, ARRAY_TYPE_ELEMENT*, ARRAY_TYPE_ELEMENT* ) = 0;


  //
  // Main solver interface
  //

  virtual boost::shared_ptr<ARRAY_TYPE_ELEMENT> solve( ARRAY_TYPE_ELEMENT *_f )
  {
    // Check that operators etc. have been provided and consistent in dimensionality
    //
    if( !validate_solver() ){
      this->solver_error( "Error: sbSolver::solve : solver validation failed");
      return boost::shared_ptr<ARRAY_TYPE_ELEMENT>();
    }

    // Define u_k
    //
    boost::shared_ptr<ARRAY_TYPE_ELEMENT> u_k( new ARRAY_TYPE_ELEMENT() );
    
    if( !u_k->create( this->encoding_operator_->get_domain_dimensions().get() )){
      this->solver_error( "Error: sbSolver::solve : memory allocation of u_k failed" );
      return boost::shared_ptr<ARRAY_TYPE_ELEMENT>();
    }        

    // Use x0 (if provided) as starting solution estimate
    //
    if( this->get_x0().get() )
      *u_k = *(this->get_x0());
    else if( !solver_clear_element( u_k.get() )){
      this->solver_error( "Error: sbSolver::solve : failed to clear u_k" );
      return boost::shared_ptr<ARRAY_TYPE_ELEMENT>();
    }
    get_inner_solver()->set_x0( u_k );

    // Normalize (a copy of) the input data
    //

    boost::shared_ptr<ARRAY_TYPE_ELEMENT> f(new ARRAY_TYPE_ELEMENT(*_f));

    if( !f || !f->get_data_ptr() ){
      this->solver_error( "Error: sbSolver::solve : memory allocation of f failed" );
      deinitialize();
      return boost::shared_ptr<ARRAY_TYPE_ELEMENT>();
    }
    
    REAL image_scale;
    if( !normalize( f, image_scale ) ){
      this->solver_error( "Error: sbSolver::solve : normalization failed" );
      deinitialize();
      return boost::shared_ptr<ARRAY_TYPE_ELEMENT>();
    }

    // Initialze d and b (see the Split Bregman paper references above) and 
    // p_M (the regularization operators' mult_M on the image prior)
    //
    boost::shared_array< boost::shared_ptr<ARRAY_TYPE_ELEMENT> > d_k, b_k, p_M;
    
    if( !initialize( image_scale, d_k, b_k, p_M ) ){
      this->solver_error( "Error: sbSolver::solve : solver initialization failed");
      return boost::shared_ptr<ARRAY_TYPE_ELEMENT>();
    }
    
    // Invoke the core solver
    //
    if( !core( tolerance_, outer_iterations_, inner_iterations_, f, u_k, d_k, b_k, p_M ) ){
      this->solver_error( "Error: sbSolver::solve : core solver failed" );
      deinitialize();
      return boost::shared_ptr<ARRAY_TYPE_ELEMENT>();
    } 

    // Undo the intermediate scaling of u_k ...
    //
    if( !undo_normalization( u_k, 1/image_scale )){
      this->solver_error( "Error: sbSolver::solve : unable to undo normalization" );
      deinitialize();
      return boost::shared_ptr<ARRAY_TYPE_ELEMENT>();
    } 
    
    // Clean up memory occupied by the operator container and inner solver
    if( !deinitialize() ){
      this->solver_error( "Error: sbSolver::solve : unable to free internal memory" );
      return boost::shared_ptr<ARRAY_TYPE_ELEMENT>();
    }
    
    // ... and return the result
    //    
    return u_k;
  }

protected:

  //
  // Everything beyond this point is internal to the implementation
  // and not intended to be exposed as a public interface
  //

  // Validate operator
  //

  virtual bool validate_encoding_operator()
  {
    boost::shared_ptr< linearOperator<REAL, ARRAY_TYPE_ELEMENT> > op = this->get_encoding_operator();

    if( !op.get() ){
      this->solver_error( "Error: sbSolver::validate_encoding_operator : operator not set" );
      return false;
    }
    
    if( op->get_domain_dimensions()->size() == 0 ){
      this->solver_error( "Error: sbSolver::validate_encoding_operator : encoding operator must have specified domain dimensions" );
      return false;
    }
    
    if( op->get_codomain_dimensions()->size() == 0 ){
      this->solver_error( "Error: sbSolver::validate_encoding_operator : encoding operator must have specified codomain dimensions" );
      return false;
    }
    
    return true;
  }
  

  // Validate regularization operator
  //

  virtual bool validate_regularization_operators( std::vector<unsigned int> *image_dims )
  {
    boost::shared_ptr< linearOperator<REAL, ARRAY_TYPE_ELEMENT> > op;
    
    // Validate regularization operators (non-group)
    for( unsigned int i=0; i<this->get_number_of_regularization_operators(); i++ ){
      
      op = this->get_regularization_operator(i);
      
      if( !op.get() ){
	this->solver_error( "Error: sbSolver::validate_regularization_operators : invalid operator provided" );
	return false;
      }
      
      if( *(op->get_domain_dimensions()) != *image_dims ){
	op->set_domain_dimensions(image_dims);
	this->solver_warning( "Warning: sbSolver::validate_regularization_operators : operator domain dimensions set to match the image dimensions" );
      }
      
      if( *(op->get_codomain_dimensions()) != *image_dims ){
	op->set_codomain_dimensions(image_dims);
	this->solver_warning( "Warning: sbSolver::validate_regularization_operators : operator codomain dimensions set to match the image dimensions" );
      }
    }
    
    // Validate regularization operator groups
    for( unsigned int i=0; i<regularization_group_operators_.size(); i++ ){
      for( unsigned int j=0; j<regularization_group_operators_[i].size(); j++ ){
	
	op = this->get_regularization_group_operator(i,j);
	
	if( !op.get() ){
	  this->solver_error( "Error: sbSolver::validate_regularization_operators : invalid operator provided" );
	  return false;
	}
	
	if( *(op->get_domain_dimensions()) != *image_dims ){
	  op->set_domain_dimensions(image_dims);
	  this->solver_warning( "Warning: sbSolver::validate_regularization_operators : operator domain dimensions set to match the image dimensions" );
	}
	
	if( *(op->get_codomain_dimensions()) != *image_dims ){
	  op->set_codomain_dimensions(image_dims);		
	  this->solver_warning( "Warning: sbSolver::validate_regularization_operators : operator codomain dimensions set to match the image dimensions" );
	}
      }
    }
    
    return true;
  }
  
  // Validate prior
  //

  virtual bool validate_prior( std::vector<unsigned int> *image_dims )
  {
    if( get_prior_image().get() && *image_dims != *(get_prior_image()->get_dimensions()) ){
      this->solver_error( "Error: sbSolver::validate_prior : prior image dimensions mismatch the encoding operator's domain dimensions" );
      return false;
    }
    
    return true;
  }
  
  // Check that the solver is set up properly
  virtual bool validate_solver()
  {
    // Some tests to see if we are ready to go...
    //
    
    if( !validate_encoding_operator() ){
      this->solver_error( "Error: sbSolver::validate : failed to validate encoding operator" );
      return false;
    }
    
    boost::shared_ptr< std::vector<unsigned int> > image_dims = this->encoding_operator_->get_domain_dimensions();
    
    if( !validate_regularization_operators( image_dims.get() )){
      this->solver_error( "Error: sbSolver::validate : failed to validate regularization operators" );
      return false;
    }
    
    if( !validate_prior( image_dims.get() )){
      this->solver_error( "Error: sbSolver::validate : failed to validate prior image" );
      return false;
    }
    
    return true;
  }
  

  // Initialize d_k and b_k arrays to zero images and compute operator priors
  //

  virtual bool initialize( REAL image_scale,
			   boost::shared_array< boost::shared_ptr<ARRAY_TYPE_ELEMENT> > &d_k, 
			   boost::shared_array< boost::shared_ptr<ARRAY_TYPE_ELEMENT> > &b_k,
			   boost::shared_array< boost::shared_ptr<ARRAY_TYPE_ELEMENT> > &p_M )
  {
    // Get image dimensions
    boost::shared_ptr< std::vector<unsigned int> > image_dims = this->encoding_operator_->get_domain_dimensions();
    
    // Determine length of arrays
    num_reg_operators_ = this->regularization_operators_.size();
    for( unsigned int i=0; i<regularization_group_operators_.size(); i++ ){
      num_reg_operators_ += regularization_group_operators_.at(i).size();
    }

    // PICCS style : if a prior is provided it generates two cost terms per regularization operator
    unsigned int l = ( get_prior_image().get() ) ? num_reg_operators_*2 : num_reg_operators_;
    
    // Allocate d_k, b_k
    d_k = boost::shared_array< boost::shared_ptr<ARRAY_TYPE_ELEMENT> > 
      ( new boost::shared_ptr<ARRAY_TYPE_ELEMENT>[l] );
    
    b_k = boost::shared_array< boost::shared_ptr<ARRAY_TYPE_ELEMENT> > 
      ( new boost::shared_ptr<ARRAY_TYPE_ELEMENT>[l] );

    if( !d_k.get() || !b_k.get() ){
      this->solver_error( "Error: sbSolver::initialize : memory allocation of d_k or b_k failed" );
      return false;
    }
    
    // Allocate prior
    if( get_prior_image().get() ){
      
      p_M = boost::shared_array< boost::shared_ptr<ARRAY_TYPE_ELEMENT> > 
	( new boost::shared_ptr<ARRAY_TYPE_ELEMENT>[num_reg_operators_] );
      
      if( !p_M.get() ){
	this->solver_error( "Error: sbSolver::initialize : memory allocation of prior terms failed" );
	return false;
      }
    }    

    // Clear d_k, b_k
    for( unsigned int i=0; i<l; i++ ){

      d_k[i] = boost::shared_ptr<ARRAY_TYPE_ELEMENT>(new ARRAY_TYPE_ELEMENT());
      
      if( !d_k[i]->create( image_dims.get() ) ){
	this->solver_error( "Error: sbSolver::initialize : memory allocation of d_k failed" );
	return false;
      }

      if( !solver_clear_element( d_k[i].get() )){
	this->solver_error( "Error: sbSolver::initialize : failed to clear internal memory buffer d_k" );
	return false;
      }

      b_k[i] = boost::shared_ptr<ARRAY_TYPE_ELEMENT>(new ARRAY_TYPE_ELEMENT());
      
      if( !b_k[i]->create( image_dims.get() ) ){
	this->solver_error( "Error: sbSolver::initialize : memory allocation of b_k failed" );
	return false;
      }

      if( !solver_clear_element( b_k[i].get() )){
	this->solver_error( "Error: sbSolver::initialize : failed to clear internal memory buffer b_k" );
	return false;
      } 
    }
    
    // Compute regularization operator :: mult_M on the prior image and store p_M
    if( get_prior_image().get() ){
      
      unsigned int operator_idx = 0;
      
      // Make copy of prior (we need to scale it)
      ARRAY_TYPE_ELEMENT tmp_prior;
      tmp_prior = *get_prior_image();

      if( !tmp_prior.get_data_ptr() || !solver_scal( image_scale, &tmp_prior ) ){
	this->solver_error( "Error: sbSolver::initialize : failed to scale prior" );
	return false;
      }
      
      // Non-group operators
      for( unsigned int i=0; i<this->get_number_of_regularization_operators(); i++ ){
      
	boost::shared_ptr< linearOperator<REAL, ARRAY_TYPE_ELEMENT> > op = 
	  this->get_regularization_operator(i);
	
	p_M[i] = boost::shared_ptr<ARRAY_TYPE_ELEMENT>(new ARRAY_TYPE_ELEMENT());
	
	if( !p_M[i]->create( image_dims.get() ) ){
	  this->solver_error( "Error: sbSolver::initialize : memory allocation of p_M failed (1)" );
	  return false;
	}
	
	if( op->mult_M( &tmp_prior, p_M[i].get() )){
	  this->solver_error( "Error: sbSolver::initialize : failed to apply operator to prior (1)" );
	  return false;
	}
	operator_idx++;
      }
            
      // Group operators
      for( unsigned int i=0; i<regularization_group_operators_.size(); i++ ){
	for( unsigned int j=0; j<regularization_group_operators_[i].size(); j++ ){
	  	  
	  boost::shared_ptr< linearOperator<REAL, ARRAY_TYPE_ELEMENT> > op = 
	    this->get_regularization_group_operator(i,j);
	  
	  p_M[operator_idx] = boost::shared_ptr<ARRAY_TYPE_ELEMENT>(new ARRAY_TYPE_ELEMENT());
	  
	  if( !p_M[operator_idx]->create( image_dims.get() ) ){
	    this->solver_error( "Error: sbSolver::initialize : memory allocation of p_M failed (2)" );
	    return false;
	  }
	  
	  if( op->mult_M( &tmp_prior, p_M[operator_idx].get() )){
	    this->solver_error( "Error: sbSolver::initialize : failed to apply operator to prior (2)" );
	    return false;
	  }
	  operator_idx++;
	}
      }
    }
    
    // Set up inner solver
    //

    if( inner_solver_->get_encoding_operator().get() ){
      this->solver_warning( "Warning: sbSolver::initialize : overriding inner solver's encoding operator" );
      this->solver_warning( "Warning: sbSolver::initialize : the encoding operator cannot be set from outside" );
    }
    
    enc_op_container_ = boost::shared_ptr<OPERATOR_CONTAINER>( new OPERATOR_CONTAINER() );
    inner_solver_->set_encoding_operator( enc_op_container_ );
    
    if( !enc_op_container_->add_operator( this->encoding_operator_ ) ){
      this->solver_error( "Error: sbSolver::initialize : failed to add encoding operator to inner solver" );
      return false;
    }
    
    // The domain is constant for all operators (the image dimensions)
    enc_op_container_->set_domain_dimensions( image_dims.get() );

    // Clear backup array before putting new stuff in
    weights_backup_.clear();
    
    // Add non-group operators to inner solver
    regularization_prior_operators_.clear();

    for( unsigned int i=0; i<this->regularization_operators_.size(); i++ ){
      
      boost::shared_ptr< linearOperator<REAL, ARRAY_TYPE_ELEMENT> > op = 
	this->get_regularization_operator(i);
      
      if( !enc_op_container_->add_operator( op ) ){
	this->solver_error( "Error: sbSolver::initialize : failed to add regularization operator to the inner solver" );		
	return false;
      }
      
      if( get_prior_image().get() ){
	
	// Make shallow copy of the operator
	boost::shared_ptr< linearOperator<REAL, ARRAY_TYPE_ELEMENT> > op_cpy = op->clone();
	regularization_prior_operators_.push_back(op_cpy);

	// Modify weight by alpha. Store the original values.
	REAL w = op->get_weight();
	weights_backup_.push_back(w);
	REAL w1 = w*(REAL(1)-get_prior_alpha());
	REAL w2 = w*get_prior_alpha();
	op->set_weight(w1);
	op_cpy->set_weight(w2);

	if( !enc_op_container_->add_operator( op_cpy ) ){
	  this->solver_error( "Error: sbSolver::initialize : failed to add regularization operator prior to the inner solver" );		
	  return false;
	}
      }
    }
    
    // Add group operators to inner solver
    regularization_group_prior_operators_.clear();
    for( unsigned int i=0; i<regularization_group_operators_.size(); i++ ){
      for( unsigned int p=0; p<((get_prior_image().get()) ? 2 : 1); p++ ){
	for( unsigned int j=0; j<regularization_group_operators_[i].size(); j++ ){
	  
	  boost::shared_ptr< linearOperator<REAL, ARRAY_TYPE_ELEMENT> > op = 
	    this->get_regularization_group_operator(i,j);
	  
	  if( p==0 ){
	    if( !enc_op_container_->add_operator( op ) ){
	      this->solver_error( "Error: sbSolver::initialize : failed to add regularization group operator to the inner solver" );		
	      return false;
	    }
	  }
	  
	  if( p==1 ){
	    
	    // Make shallow copy of the operator
	    boost::shared_ptr< linearOperator<REAL, ARRAY_TYPE_ELEMENT> > op_cpy = op->clone();
	    regularization_group_prior_operators_.push_back(op_cpy);
	    
	    // Modify weight by alpha. Store the original values.
	    REAL w = op->get_weight();
	    weights_backup_.push_back(w);
	    REAL w1 = w*(REAL(1)-get_prior_alpha());
	    REAL w2 = w*get_prior_alpha();
	    op->set_weight(w1);
	    op_cpy->set_weight(w2);	  
	    
	    if( !enc_op_container_->add_operator( op_cpy ) ){
	      this->solver_error( "Error: sbSolver::initialize : failed to add regularization group operator to the inner solver" );
	      return false;
	    }
	  }
	}
      }     
    }
    
    return true;
  }
  
  // Clean up operator memory in the inner solver
  // Also restore the weights we temporarily changed

  virtual bool deinitialize() 
  {
    enc_op_container_ = boost::shared_ptr<OPERATOR_CONTAINER>( new OPERATOR_CONTAINER() );
    inner_solver_->set_encoding_operator( enc_op_container_ );

    if( get_prior_image().get() ){
      
      unsigned int operator_idx = 0;
      
      for( unsigned int i=0; i<this->regularization_operators_.size(); i++ ){
	
	boost::shared_ptr< linearOperator<REAL, ARRAY_TYPE_ELEMENT> > op = 
	  this->get_regularization_operator(i);
	
	op->set_weight( weights_backup_[operator_idx] );		
	operator_idx++;
      }

      for( unsigned int i=0; i<regularization_group_operators_.size(); i++ ){
	for( unsigned int j=0; j<regularization_group_operators_[i].size(); j++ ){

	  boost::shared_ptr< linearOperator<REAL, ARRAY_TYPE_ELEMENT> > op = 
	    this->get_regularization_group_operator(i,j);
	  
	  op->set_weight( weights_backup_[operator_idx] );		
	  operator_idx++;
	}
      }
    }

    regularization_prior_operators_.clear();
    regularization_group_prior_operators_.clear();

    return true;
  }

      
  // Normalize data
  //

  virtual bool normalize( boost::shared_ptr<ARRAY_TYPE_ELEMENT> f, 
			  REAL &image_scale )    
  {    
    // Normalize in image space
    //

    ARRAY_TYPE_ELEMENT tmp;
    if( !tmp.create( this->encoding_operator_->get_domain_dimensions().get() )){
      this->solver_error( "Error: sbSolver::normalize : memory allocation failed" );
      return false;
    }    
    
    if( this->encoding_operator_->mult_MH( f.get(), &tmp ) < 0 ){
      this->solver_error( "Error: sbSolver::normalize : adjoint encoding operation failed on f" );
      return false;
    }

    // Normalize to an average energy of "one intensity unit per image element"
    //

    REAL sum = solver_asum( &tmp );
    image_scale = (REAL) (tmp.get_number_of_elements())/REAL(sum);

    if(	!solver_scal( image_scale, f.get() )){
      this->solver_error( "Error: sbSolver::normalize : unable to scale f" );
      return false;
    }

    return true;
  }


  // Undo normalization
  //

  virtual bool undo_normalization( boost::shared_ptr<ARRAY_TYPE_ELEMENT> u_k, ELEMENT_TYPE undo_scale )
  {
    if(	!solver_scal( undo_scale, u_k.get() )){
      this->solver_error( "Error: sbSolver::undo_normalization : unable to undo scaling of u_k" );
      return false;
    }
    
    return true;
  }


  // The core of the Split Bregman solver.
  //

  virtual bool core( REAL tolerance, unsigned int outer_iterations, unsigned int inner_iterations,
		     boost::shared_ptr<ARRAY_TYPE_ELEMENT> f, 
		     boost::shared_ptr<ARRAY_TYPE_ELEMENT> u_k,
		     boost::shared_array< boost::shared_ptr<ARRAY_TYPE_ELEMENT> > d_k,
		     boost::shared_array< boost::shared_ptr<ARRAY_TYPE_ELEMENT> > b_k,
		     boost::shared_array< boost::shared_ptr<ARRAY_TYPE_ELEMENT> > p_M )
  {
    // Image space dimensions
    boost::shared_ptr< std::vector<unsigned int> > image_dims = 
      this->encoding_operator_->get_domain_dimensions();
    
    // Keep a copy of the "previous" u_k to compute the outer loop change of u_k
    // 

    ARRAY_TYPE_ELEMENT u_k_prev;
    if( tolerance > REAL(0) || 
	this->output_mode_ >= solver<ARRAY_TYPE_ELEMENT, ARRAY_TYPE_ELEMENT>::OUTPUT_VERBOSE ){

      u_k_prev = *u_k;
      if( !u_k_prev.get_data_ptr() ){
	this->solver_error( "Error: sbSolver::core : memory allocation of u_k_prev failed" );
	return false;
      }
    }

    //
    // Outer loop
    //

    for( unsigned int outer_iteration=0; outer_iteration<outer_iterations; outer_iteration++ ) {

      if( this->output_mode_ >= solver<ARRAY_TYPE_ELEMENT, ARRAY_TYPE_ELEMENT>::OUTPUT_MAX )
	std::cout << std::endl << "SB outer loop iteration " << outer_iteration << std::endl << std::endl;

      //
      // Inner loop
      //

      for( unsigned int inner_iteration=0; inner_iteration<inner_iterations; inner_iteration++ ) {

	if( this->output_mode_ >= solver<ARRAY_TYPE_ELEMENT, ARRAY_TYPE_ELEMENT>::OUTPUT_MAX )
	  std::cout << std::endl << "SB inner loop iteration " << inner_iteration << std::endl << std::endl;
	
	// Setup input vector to the encoding operator container (argument to the inner solver's solve)
	//	

	{ // Brackets used to free 'data' below as soon as it goes out of scope

	  ARRAY_TYPE_ELEMENT data, tmp;
	  if( !data.create(enc_op_container_->get_codomain_dimensions().get()) ){
	    this->solver_error( "Error: sbSolver::core : memory allocation for container operator data failed (1)" );
	    return false;
	  }
	  
	  // First add the encoding operator data, f
	  //
	  
	  if( !tmp.create( f->get_dimensions().get(), data.get_data_ptr() )){
	    this->solver_error( "Error: sbSolver::core : memory allocation for container operator data failed (2)" );
	    return false;
	  }
	  
	  tmp = *f;	  
	  
	  // Next add the regularization operators' data, d_k - b_k
	  //

	  unsigned int operator_idx = 0;

	  for( unsigned int i=0; i<this->regularization_operators_.size(); i++ ){
	    for( unsigned int p=0; p<((get_prior_image().get()) ? 2 : 1); p++ ){
	      
	      if( !tmp.create( image_dims.get(), data.get_data_ptr()+enc_op_container_->get_offset(operator_idx+1) )){
		this->solver_error( "Error: sbSolver::core : memory allocation for container operator data failed (3)" );
		return false;
	      }
	      
	      tmp = *d_k[operator_idx];
	      
	      if( !solver_axpy_element( ELEMENT_TYPE(-1), b_k[operator_idx].get(), &tmp )){
		this->solver_error( "Error: sbSolver::core : axpy on difference image failed (1)" );
		return false;
	      }
	      
	      if( p == 1 ){
		if( !solver_axpy_element( ELEMENT_TYPE(1), p_M[operator_idx/2].get(), &tmp )){
		  this->solver_error( "Error: sbSolver::core : axpy on difference image failed (2)" );
		  return false;
		}
	      }
	      operator_idx++;
	    }
	  }

	  for( unsigned int i=0; i<this->regularization_group_operators_.size(); i++ ){
	    for( unsigned int p=0; p<((get_prior_image().get()) ? 2 : 1); p++ ){
	      for( unsigned int j=0; j<regularization_group_operators_[i].size(); j++ ){
		
		if( !tmp.create( image_dims.get(), data.get_data_ptr()+enc_op_container_->get_offset(operator_idx+1) )){
		  this->solver_error( "Error: sbSolver::core : memory allocation for container operator data failed (3)" );
		  return false;
		}
		
		tmp = *d_k[operator_idx];
		
		if( !solver_axpy_element( ELEMENT_TYPE(-1), b_k[operator_idx].get(), &tmp )){
		  this->solver_error( "Error: sbSolver::core : axpy on difference image failed (1)" );
		  return false;
		}
		
		if( p == 1 ){
		  if( !solver_axpy_element( ELEMENT_TYPE(1), p_M[operator_idx/2].get(), &tmp )){
		    this->solver_error( "Error: sbSolver::core : axpy on difference image failed (2)" );
		    return false;
		  }
		}
		operator_idx++;
	      }
	    }
	  }
	  
	  // Solve for u_k
	  //
	  
	  {
	    boost::shared_ptr<ARRAY_TYPE_ELEMENT> tmp_u_k = get_inner_solver()->solve( &data );

	    // Compute change in u_k
	    if( this->output_mode_ >= solver<ARRAY_TYPE_ELEMENT, ARRAY_TYPE_ELEMENT>::OUTPUT_VERBOSE ){
	      if( !solver_axpy_element( ELEMENT_TYPE(-1), tmp_u_k.get(), u_k.get() )){
		this->solver_error( "Error: sbSolver::core : error computing inner loop u_k delta" );
		return false;
	      }
	      std::cout << std::endl << "u_k delta (inner loop): " << solver_asum(u_k.get()) << std::endl;
	    }
	    
	    // Update u_k 
	    *u_k = *tmp_u_k;
	  }
	}

	unsigned int operator_idx = 0;

	// Update d_k (and in the final inner iteration b_k also)
	//

	for( unsigned int i=0; i<this->regularization_operators_.size(); i++ ){
	  for( unsigned int p=0; p<((get_prior_image().get()) ? 2 : 1); p++ ){
	    
	    boost::shared_ptr< linearOperator<REAL, ARRAY_TYPE_ELEMENT> > op = 
	      (p==0 ) ? this->regularization_operators_.at(i) : this->regularization_prior_operators_.at(i);
	    
	    ARRAY_TYPE_ELEMENT tmp_sum, reg_out;
	    if( !tmp_sum.create( image_dims.get() ) || !reg_out.create( image_dims.get() ) ){
	      this->solver_error( "Error: sbSolver::core : memory allocation for regularization operator failed in {d_k,b_k} update" );
	      return false;
	    }
	    
	    tmp_sum = *b_k[operator_idx];
	    
	    if( op->mult_M( u_k.get(), &reg_out ) < 0 ){
	      this->solver_error( "Error: sbSolver::core : application of regularization operator failed in {d_k,b_k} update" );
	      return false;
	    }
	    
	    if( p==1 ){
	      if( !solver_axpy_element( ELEMENT_TYPE(-1), p_M[operator_idx/2].get(), &reg_out )){
		this->solver_error( "Error: sbSolver::core : computation of shrinkage_1 argument for d_k failed (2)" );
		return false;
	      }
	    }	  
	    
	    if( !solver_axpy_element( ELEMENT_TYPE(1), &reg_out, &tmp_sum )){
	      this->solver_error( "Error: sbSolver::core : computation of shrinkage_1 argument for d_k failed (1)" );
	      return false;
	    }
	    
	    // Update of d_k

	    if( !solver_shrink1( REAL(1)/op->get_weight(), &tmp_sum, d_k[operator_idx].get() )){
	      this->solver_error( "Error: sbSolver::core : shrinkage_1 of d_k failed" );
	      return false;
	    }
	    
	    // Update of b_k (only in the last inner iteration)
	    if( inner_iteration == inner_iterations-1 ){
	      
	      if( !solver_axpy_element( ELEMENT_TYPE(-1), d_k[operator_idx].get(), &reg_out )){
		this->solver_error( "Error: sbSolver::core : computation of update argument to b_k failed" );
		return false;
	      }
	      
	      // Update of b_k
	      if( !solver_axpy_element( ELEMENT_TYPE(1), &reg_out, b_k[operator_idx].get() )){
		this->solver_error( "Error: sbSolver::core : update of b_k failed" );
		return false;
	      }
	    }
	    operator_idx++;
	  }
	}
	
	unsigned int group_idx = 0;
	  
	for( unsigned int i=0; i<regularization_group_operators_.size(); i++ ){

	  unsigned int k = regularization_group_operators_.at(i).size();

	  ARRAY_TYPE_ELEMENT *sums    = new ARRAY_TYPE_ELEMENT[k*((get_prior_image().get()) ? 2 : 1)];
	  ARRAY_TYPE_ELEMENT *reg_out = new ARRAY_TYPE_ELEMENT[k*((get_prior_image().get()) ? 2 : 1)];
	  ARRAY_TYPE_REAL    *s_k     = new ARRAY_TYPE_REAL   [  ((get_prior_image().get()) ? 2 : 1)];

	  if( !sums || !reg_out || !s_k ){
	    this->solver_error( "Error: sbSolver::core : host memory allocation for temporary arrays failed" );
	    return false;
	  }

	  for( unsigned int p=0; p<((get_prior_image().get()) ? 2 : 1); p++ ){
	    if( s_k[p].create(image_dims.get()) < 0 ){
	      this->solver_error( "Error: sbSolver::core : memory allocation for s_k failed" );
	      return false;
	    }
	    
	    if( !solver_clear_real(&s_k[p]) ){
	      this->solver_error( "Error: sbSolver::core : failed to clear s_k" );
	      return false;
	    } 	  
	  }
	  
	  for( unsigned int p=0; p<((get_prior_image().get()) ? 2 : 1); p++ ){
	    for( unsigned int j=0; j<k; j++ ){
	      
	      unsigned int idx = p*k+j;

	      boost::shared_ptr< linearOperator<REAL, ARRAY_TYPE_ELEMENT> > op = 
		(p==0 ) ? this->regularization_group_operators_.at(i).at(j) : 
		this->regularization_group_prior_operators_.at(group_idx+j);	      
	      
	      if( sums[idx].create( image_dims.get() ) < 0 || reg_out[idx].create( image_dims.get() ) < 0 ){
		this->solver_error( "Error: sbSolver::core : memory allocation for regularization operator failed in {d_k,b_k} update" );
		return false;
	      }
	      
	      sums[idx] = *b_k[operator_idx+idx];
	      
	      if( op->mult_M( u_k.get(), &reg_out[idx] ) < 0 ){
		this->solver_error( "Error: sbSolver::core : application of regularization operator failed in {d_k,b_k} update" );
		return false;
	      }
	      
	      if( p==1 ){
		if( !solver_axpy_element( ELEMENT_TYPE(-1), p_M[operator_idx/2+j].get(), &reg_out[idx] )){
		  this->solver_error( "Error: sbSolver::core : computation of shrinkage_d argument for d_k failed (2)" );
		  return false;
		}
	      }	  
	      
	      if( !solver_axpy_element( ELEMENT_TYPE(1), &reg_out[idx], &sums[idx] )){
		this->solver_error( "Error: sbSolver::core : computation of shrinkage_d argument for d_k failed" );
		return false;
	      }
	      
	      boost::shared_ptr<ARRAY_TYPE_REAL> tmp_s_k = solver_norm( &sums[idx] );
	      if( !solver_axpy_real( REAL(1), tmp_s_k.get(), &s_k[p] )){
		this->solver_error( "Error: sbSolver::core : accumulation of s_k failed" );
		return false;
	      }
	    }
	    
	    if( !solver_sqrt( &s_k[p] )){
	      this->solver_error( "Error: sbSolver::core : sqrt of s_k failed" );
	      return false;
	    }
	  }

	  for( unsigned int p=0; p<((get_prior_image().get()) ? 2 : 1); p++ ){
	    for( unsigned int j=0; j<k; j++ ){
	      
	      unsigned int idx = p*k+j;

	      boost::shared_ptr< linearOperator<REAL, ARRAY_TYPE_ELEMENT> > op = 
		(p==0 ) ? this->regularization_group_operators_.at(i).at(j) : 
		this->regularization_group_prior_operators_.at(group_idx+j);
	      
	      // Update of d_k
	      if( !solver_shrinkd( REAL(1)/op->get_weight(), &s_k[p], &sums[idx], d_k[operator_idx+idx].get() )){
		this->solver_error( "Error: sbSolver::core : shrinkage_d of d_k failed" );
		return false;
	      }
	      
	      // Update of b_k (only in the last inner iteration)
	      if( inner_iteration == inner_iterations-1 ){
		
		if( !solver_axpy_element( ELEMENT_TYPE(-1), d_k[operator_idx+idx].get(), &reg_out[idx] )){
		  this->solver_error( "Error: sbSolver::core : computation of update argument to b_k failed" );
		  return false;
		}
		
		if( !solver_axpy_element( ELEMENT_TYPE(1), &reg_out[idx], b_k[operator_idx+idx].get() )){
		  this->solver_error( "Error: sbSolver::core : update of b_k failed" );
		  return false;
		}
	      }
	    }
	  }
	  operator_idx += (k*((get_prior_image().get()) ? 2 : 1)); group_idx += k;
	  delete[] sums; delete[] reg_out; delete[] s_k;
	}

      } // end of inner loop

      // Output change in u_k
      if( tolerance > REAL(0) || this->output_mode_ >= solver<ARRAY_TYPE_ELEMENT, ARRAY_TYPE_ELEMENT>::OUTPUT_VERBOSE ){

	if( !solver_scal( ELEMENT_TYPE(-1), &u_k_prev ) ){
	  this->solver_error( "Error: sbSolver::core : error computing inner loop u_k delta (scale)" );
	  return false;
	}

	if( !solver_axpy_element( ELEMENT_TYPE(1), u_k.get(), &u_k_prev )){
	  this->solver_error( "Error: sbSolver::core : error computing inner loop u_k delta (axpy)" );
	  return false;
	}

	REAL delta = solver_asum(&u_k_prev);

	if( this->output_mode_ >= solver<ARRAY_TYPE_ELEMENT, ARRAY_TYPE_ELEMENT>::OUTPUT_VERBOSE )
	  std::cout << std::endl << "u_k delta (outer loop): " << delta << std::endl << std::endl;

	if( delta < tolerance )
	  break;

	u_k_prev = *u_k;
      }

    } // end of outer loop

    return true;
  }

protected:
  REAL tolerance_;
  unsigned int outer_iterations_, inner_iterations_;
  unsigned int num_reg_operators_;
  std::vector< boost::shared_ptr< linearOperator<REAL, ARRAY_TYPE_ELEMENT> > > _regularization_group_operators_;
  std::vector< std::vector< boost::shared_ptr< linearOperator<REAL, ARRAY_TYPE_ELEMENT> > > > regularization_group_operators_;
  std::vector< boost::shared_ptr< linearOperator<REAL, ARRAY_TYPE_ELEMENT> > > regularization_prior_operators_;
  std::vector< boost::shared_ptr< linearOperator<REAL, ARRAY_TYPE_ELEMENT> > > regularization_group_prior_operators_;
  boost::shared_ptr<ARRAY_TYPE_ELEMENT> prior_;
  REAL alpha_;
  boost::shared_ptr<INNER_SOLVER> inner_solver_;
  boost::shared_ptr<OPERATOR_CONTAINER> enc_op_container_;
  std::vector<unsigned int> weights_backup_;
};
