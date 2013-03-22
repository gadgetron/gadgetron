/** \file sbSolver.h
    \brief Implementation of an unconstraint Split Bregman solver

    The file sbSolver.h is a device independent implementation of the unconstraint Split Bregman solver. 
    For a constraint solver version see instead the file sbcSolver.h.
    Detail on the algorithm can be found in sec. 3.2. of the paper
    "The Split Bregman Method for L1-Regularized Problems" 
    by Tom Goldstein and Stanley Osher. Siam J. Imaging Sciences. Vol. 2, No. 2, pp. 323-343.

    To simplify the actual solver instantiation we refer to the files
    - the class(/file) hoSbSolver(/.h) for a cpu instantiated solver using the hoNDArray class
    - the class(/file) cuSbSolver(/.h) for a gpu instantiated solver using the cuNDArray class
    - the class(/file) hoCuSbSolver(/.h) for a gpu based solver using a host memory interface. 

    The latter version is intended for large reconstructions in which device memory cannot hold 
    the entire data from the image and encoded image domains. 
    In the "hoCu" scenario, suitable encoding and regularization operators
    capable of batching their mult_M and mult_MHM functions should be chosen.

    In all cases, the encoding and regularization operators added to the solver 
    must adhere to the underlying instantiation of the NDArray data type.
*/

#pragma once

#include "linearSolver.h"
#include "vector_td_utilities.h"
#include "encodingOperatorContainer.h"

#include <vector>
#include <iostream>
#include <algorithm>

namespace Gadgetron{

template<
	  class ARRAY_TYPE_REAL, 
	  class ARRAY_TYPE_ELEMENT, 
	  class INNER_SOLVER>

class sbSolver : public linearSolver<ARRAY_TYPE_ELEMENT>
{
	 typedef typename ARRAY_TYPE_ELEMENT::element_type ELEMENT_TYPE;
	 typedef typename realType<ELEMENT_TYPE>::type REAL;
	 typedef encodingOperatorContainer<ARRAY_TYPE_ELEMENT> OPERATOR_CONTAINER;
public:

  // Constructor
  //

  sbSolver() : linearSolver<ARRAY_TYPE_ELEMENT>()
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

  virtual bool add_regularization_group_operator ( boost::shared_ptr< linearOperator< ARRAY_TYPE_ELEMENT> > op )
  {
    if( !op.get() ){
      throw std::runtime_error( "Error: sbSolver::add_regularization_group_operator : NULL operator provided" );
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
      throw std::runtime_error( "Error: sbSolver::add_group : no regularization group operators added" );
      return false;
    }

    regularization_group_operators_.push_back( _regularization_group_operators_ );    
    _regularization_group_operators_.clear();
    
    return true;
  }
  

  // Get regularization group operator
  //

  virtual boost::shared_ptr< linearOperator<ARRAY_TYPE_ELEMENT> >
  get_regularization_group_operator( unsigned int group_number, unsigned int idx_in_group )
  {
    if( group_number >= regularization_group_operators_.size() ){
      throw std::runtime_error( "Error: sbSolver::get_regularization_group_operator : group number out of range" );
      return boost::shared_ptr< linearOperator<ARRAY_TYPE_ELEMENT> >();
    }
    
    if( idx_in_group >= regularization_group_operators_.at(group_number).size() ){
      throw std::runtime_error( "Error: sbSolver::get_regularization_group_operator : idx in group out of range" );
      return boost::shared_ptr< linearOperator<ARRAY_TYPE_ELEMENT> >();
    }
    
    return regularization_group_operators_.at(group_number).at(idx_in_group);
  }  

  
  // Set/get prior image (PICCS style). 
  // I.e. for every regularization operator (group) R that is added we minimize:
  // alpha|R(x-prior)|_{l1} + (1-alpha)|R(x)|_{l1}
  //

  virtual void set_prior_image( boost::shared_ptr<ARRAY_TYPE_ELEMENT> prior, REAL alpha )
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

  


  // Provide the user an option to access u_k right after its update. 
  //

  virtual bool post_linear_solver_callback( ARRAY_TYPE_ELEMENT* ) { return true; }

  //
  // Main solver interface
  //

  virtual boost::shared_ptr<ARRAY_TYPE_ELEMENT> solve( ARRAY_TYPE_ELEMENT *_f )
  {
    // Check that operators etc. have been provided and consistent in dimensionality
    //
    validate_solver();
    // Define u_k
    //
    boost::shared_ptr<ARRAY_TYPE_ELEMENT> u_k( new ARRAY_TYPE_ELEMENT(this->encoding_operator_->get_domain_dimensions()) );
    


    // Use x0 (if provided) as starting solution estimate
    //
    if( this->get_x0().get() )
      *u_k = *(this->get_x0());
    u_k->clear();

    //    get_inner_solver()->set_x0( u_k );

    // Normalize (a copy of) the input data
    //

    boost::shared_ptr<ARRAY_TYPE_ELEMENT> f(new ARRAY_TYPE_ELEMENT(*_f));

    if( !f || !f->get_data_ptr() ){
    	deinitialize();
      throw std::runtime_error( "Error: sbSolver::solve : memory allocation of f failed" );


    }
    
    REAL image_scale;
    normalize( f, image_scale );

    // Initialze d and b (see the Split Bregman paper references above) and 
    // p_M (the regularization operators' mult_M on the image prior)
    //
    boost::shared_array< boost::shared_ptr<ARRAY_TYPE_ELEMENT> > d_k, b_k, p_M;
    
    initialize( image_scale, d_k, b_k, p_M );
    
    // Invoke the core solver
    //
    core( tolerance_, outer_iterations_, inner_iterations_, f, u_k, d_k, b_k, p_M );

    // Undo the intermediate scaling of u_k ...
    //
    undo_normalization( u_k, 1/image_scale );
    
    // Clean up memory occupied by the operator container and inner solver
    deinitialize();
    
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

  virtual void validate_encoding_operator()
  {
    boost::shared_ptr< linearOperator<ARRAY_TYPE_ELEMENT> > op = this->get_encoding_operator();

    if( !op.get() ){
      throw std::runtime_error( "Error: sbSolver::validate_encoding_operator : operator not set" );

    }
    
    if( op->get_domain_dimensions()->size() == 0 ){
      throw std::runtime_error( "Error: sbSolver::validate_encoding_operator : encoding operator must have specified domain dimensions" );

    }
    
    if( op->get_codomain_dimensions()->size() == 0 ){
      throw std::runtime_error( "Error: sbSolver::validate_encoding_operator : encoding operator must have specified codomain dimensions" );

    }
    
  }
  

  // Validate regularization operator
  //

  virtual bool validate_regularization_operators( std::vector<unsigned int> *image_dims )
  {
    boost::shared_ptr< linearOperator<ARRAY_TYPE_ELEMENT> > op;
    
    // Validate regularization operators (non-group)
    for( unsigned int i=0; i<this->get_number_of_regularization_operators(); i++ ){
      
      op = this->get_regularization_operator(i);
      
      if( !op.get() ){
      	throw std::runtime_error( "Error: sbSolver::validate_regularization_operators : invalid operator provided" );

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
	  throw std::runtime_error( "Error: sbSolver::validate_regularization_operators : invalid operator provided" );

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

  virtual void validate_prior( std::vector<unsigned int> *image_dims )
  {
    if( get_prior_image().get() && *image_dims != *(get_prior_image()->get_dimensions()) ){
      throw std::runtime_error( "Error: sbSolver::validate_prior : prior image dimensions mismatch the encoding operator's domain dimensions" );
    }
    
  }
  
  // Check that the solver is set up properly
  virtual void validate_solver()
  {
    // Some tests to see if we are ready to go...
    //
    
    validate_encoding_operator();
    
    boost::shared_ptr< std::vector<unsigned int> > image_dims = this->encoding_operator_->get_domain_dimensions();
    validate_regularization_operators( image_dims.get() );
    validate_prior( image_dims.get() );
    
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
      throw std::runtime_error( "Error: sbSolver::initialize : memory allocation of d_k or b_k failed" );
      return false;
    }
    
    // Allocate prior
    if( get_prior_image().get() ){
      
      p_M = boost::shared_array< boost::shared_ptr<ARRAY_TYPE_ELEMENT> > 
	( new boost::shared_ptr<ARRAY_TYPE_ELEMENT>[num_reg_operators_] );
      
      if( !p_M.get() ){
	throw std::runtime_error( "Error: sbSolver::initialize : memory allocation of prior terms failed" );
	return false;
      }
    }    

    // Clear d_k, b_k
    for( unsigned int i=0; i<l; i++ ){

      d_k[i] = boost::shared_ptr<ARRAY_TYPE_ELEMENT>(new ARRAY_TYPE_ELEMENT(image_dims));
      d_k[i]->clear();

      b_k[i] = boost::shared_ptr<ARRAY_TYPE_ELEMENT>(new ARRAY_TYPE_ELEMENT(image_dims));
      b_k[i]->clear();

    }
    
    // Compute regularization operator :: mult_M on the prior image and store p_M
    if( get_prior_image().get() ){
      
      unsigned int operator_idx = 0;
      
      // Make copy of prior (we need to scale it)
      ARRAY_TYPE_ELEMENT tmp_prior;
      tmp_prior = *get_prior_image();
      tmp_prior *= image_scale;
      // Non-group operators
      for( unsigned int i=0; i<this->get_number_of_regularization_operators(); i++ ){
      
	boost::shared_ptr< linearOperator<ARRAY_TYPE_ELEMENT> > op =
	  this->get_regularization_operator(i);
	
	p_M[i] = boost::shared_ptr<ARRAY_TYPE_ELEMENT>(new ARRAY_TYPE_ELEMENT(image_dims));
	
	op->mult_M( &tmp_prior, p_M[i].get() );
	operator_idx++;
      }
            
      // Group operators
      for( unsigned int i=0; i<regularization_group_operators_.size(); i++ ){
	for( unsigned int j=0; j<regularization_group_operators_[i].size(); j++ ){
	  	  
	  boost::shared_ptr< linearOperator<ARRAY_TYPE_ELEMENT> > op =
	    this->get_regularization_group_operator(i,j);
	  
	  p_M[operator_idx] = boost::shared_ptr<ARRAY_TYPE_ELEMENT>(new ARRAY_TYPE_ELEMENT(image_dims));
	  
	  
	   op->mult_M( &tmp_prior, p_M[operator_idx].get());
	  operator_idx++;
	}
      }
    }
    
    // Set up inner solver
    //

    enc_op_container_ = boost::shared_ptr<OPERATOR_CONTAINER>( new OPERATOR_CONTAINER() );
    inner_solver_->set_encoding_operator( enc_op_container_ );
    
    enc_op_container_->add_operator( this->encoding_operator_ );
    
    // The domain is constant for all operators (the image dimensions)
    enc_op_container_->set_domain_dimensions( image_dims.get() );

    // Clear backup array before putting new stuff in
    weights_backup_.clear();
    
    // Add non-group operators to inner solver
    regularization_prior_operators_.clear();

    for( unsigned int i=0; i<this->regularization_operators_.size(); i++ ){
      
      boost::shared_ptr< linearOperator<ARRAY_TYPE_ELEMENT> > op =
	this->get_regularization_operator(i);
      
      enc_op_container_->add_operator( op );
      
      if( get_prior_image().get() ){
	
				// Make shallow copy of the operator
				boost::shared_ptr< linearOperator<ARRAY_TYPE_ELEMENT> > op_cpy = op->clone();
				regularization_prior_operators_.push_back(op_cpy);

				// Modify weight by alpha. Store the original values.
				REAL w = op->get_weight();
				weights_backup_.push_back(w);
				REAL w1 = w*(REAL(1)-get_prior_alpha());
				REAL w2 = w*get_prior_alpha();
				op->set_weight(w1);
				op_cpy->set_weight(w2);

				enc_op_container_->add_operator( op_cpy );
      }
    }
    
    // Add group operators to inner solver
    regularization_group_prior_operators_.clear();
    for( unsigned int i=0; i<regularization_group_operators_.size(); i++ ){
      for( unsigned int p=0; p<((get_prior_image().get()) ? 2 : 1); p++ ){
	for( unsigned int j=0; j<regularization_group_operators_[i].size(); j++ ){
	  
	  boost::shared_ptr< linearOperator<ARRAY_TYPE_ELEMENT> > op =
	    this->get_regularization_group_operator(i,j);
	  
	  if( p==0 ){
	    enc_op_container_->add_operator( op );
	  }
	  
	  if( p==1 ){
	    
	    // Make shallow copy of the operator
	    boost::shared_ptr< linearOperator<ARRAY_TYPE_ELEMENT> > op_cpy = op->clone();
	    regularization_group_prior_operators_.push_back(op_cpy);
	    
	    // Modify weight by alpha. Store the original values.
	    REAL w = op->get_weight();
	    weights_backup_.push_back(w);
	    REAL w1 = w*(REAL(1)-get_prior_alpha());
	    REAL w2 = w*get_prior_alpha();
	    op->set_weight(w1);
	    op_cpy->set_weight(w2);	  
	    
	    enc_op_container_->add_operator( op_cpy );
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
	
	boost::shared_ptr< linearOperator<ARRAY_TYPE_ELEMENT> > op =
	  this->get_regularization_operator(i);
	
	op->set_weight( weights_backup_[operator_idx] );		
	operator_idx++;
      }

      for( unsigned int i=0; i<regularization_group_operators_.size(); i++ ){
	for( unsigned int j=0; j<regularization_group_operators_[i].size(); j++ ){

	  boost::shared_ptr< linearOperator<ARRAY_TYPE_ELEMENT> > op =
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

  virtual void normalize( boost::shared_ptr<ARRAY_TYPE_ELEMENT> f,
			  REAL &image_scale )    
  {    
    // Normalize in image space
    //

    ARRAY_TYPE_ELEMENT tmp(this->encoding_operator_->get_domain_dimensions());

    
     this->encoding_operator_->mult_MH( f.get(), &tmp );

    // Normalize to an average energy of "one intensity unit per image element"
    //

    REAL sum = asum( &tmp );
    //image_scale = (REAL) (tmp.get_number_of_elements())/REAL(sum);

    // image_scale = std::max( REAL(tmp.get_number_of_elements()), REAL(1e6) )/REAL(sum);
    if ( REAL(tmp.get_number_of_elements()) > REAL(1e6) )
    {
        image_scale = REAL(tmp.get_number_of_elements())/REAL(sum);
    }
    else
    {
        image_scale = REAL(1e6)/REAL(sum);
    }

    *f *= image_scale;



  }


  // Undo normalization
  //

  virtual void undo_normalization( boost::shared_ptr<ARRAY_TYPE_ELEMENT> u_k, ELEMENT_TYPE undo_scale )
  {
  	*u_k *= undo_scale;
    

  }


  // The core of the Split Bregman solver.
  //

  virtual void core( REAL tolerance, unsigned int outer_iterations, unsigned int inner_iterations,
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

	  ARRAY_TYPE_ELEMENT data(enc_op_container_->get_codomain_dimensions());
	  ARRAY_TYPE_ELEMENT tmp(f->get_dimensions().get(), data.get_data_ptr());
	  
	  // First add the encoding operator data, f
	  //
	  
	  
	  tmp = *f;	  
	  
	  // Next add the regularization operators' data, d_k - b_k
	  //

	  unsigned int operator_idx = 0;

	  for( unsigned int i=0; i<this->regularization_operators_.size(); i++ ){
	    for( unsigned int p=0; p<((get_prior_image().get()) ? 2 : 1); p++ ){
	      
	      tmp.create( image_dims.get(), data.get_data_ptr()+enc_op_container_->get_offset(operator_idx+1));
	      
	      tmp = *d_k[operator_idx];
	      tmp -= *b_k[operator_idx];
	      
	      if( p == 1 ){
	      	tmp +=  *p_M[operator_idx/2];
	      }
	      operator_idx++;
	    }
	  }

	  for( unsigned int i=0; i<this->regularization_group_operators_.size(); i++ ){
	    for( unsigned int p=0; p<((get_prior_image().get()) ? 2 : 1); p++ ){
	      for( unsigned int j=0; j<regularization_group_operators_[i].size(); j++ ){
		
		tmp.create( image_dims.get(), data.get_data_ptr()+enc_op_container_->get_offset(operator_idx+1) );
		
		tmp = *d_k[operator_idx];
		tmp -= *b_k[operator_idx];
		
		if( p == 1 ){
			tmp += *p_M[operator_idx/2];
		}
		operator_idx++;
	      }
	    }
	  }
	  
	  // Solve for u_k
	  //
	  
	  {
	    boost::shared_ptr<ARRAY_TYPE_ELEMENT> tmp_u_k = get_inner_solver()->solve( &data );

	    // Invoke the post inner solver callback
	    if( !post_linear_solver_callback( tmp_u_k.get() ) ){
	      throw std::runtime_error( "Error: sbSolver::core : error computing inner loop u_k delta" );

	    }
	    
	    // Compute change in u_k
	    if( this->output_mode_ >= solver<ARRAY_TYPE_ELEMENT, ARRAY_TYPE_ELEMENT>::OUTPUT_VERBOSE ){
	    	*u_k -= *tmp_u_k;
	      std::cout << "u_k delta (inner loop): " << asum(u_k.get()) << std::endl;
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
	    
	    boost::shared_ptr< linearOperator<ARRAY_TYPE_ELEMENT> > op =
	      (p==0 ) ? this->regularization_operators_.at(i) : this->regularization_prior_operators_.at(i);
	    
	    ARRAY_TYPE_ELEMENT tmp_sum(image_dims);
	    ARRAY_TYPE_ELEMENT reg_out(image_dims);
	    
	    tmp_sum = *b_k[operator_idx];
	    
	    op->mult_M( u_k.get(), &reg_out );
	    
	    if( p==1 ){
	    	reg_out -= *p_M[operator_idx/2];
	    }	  
	    reg_out += tmp_sum;
	    
	    // Update of d_k

	    shrink1( REAL(1)/op->get_weight(), &tmp_sum, d_k[operator_idx].get() );
	    
	    // Update of b_k (only in the last inner iteration)
	    if( inner_iteration == inner_iterations-1 ){
	      reg_out -= *d_k[operator_idx];
	      
	      // Update of b_k
	      reg_out += *b_k[operator_idx];
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
	    throw std::runtime_error( "Error: sbSolver::core : host memory allocation for temporary arrays failed" );

	  }

	  for( unsigned int p=0; p<((get_prior_image().get()) ? 2 : 1); p++ ){
	    s_k[p].create(image_dims.get());
	    s_k[p].clear();
	  }
	  
	  for( unsigned int p=0; p<((get_prior_image().get()) ? 2 : 1); p++ ){
	    for( unsigned int j=0; j<k; j++ ){
	      
	      unsigned int idx = p*k+j;

	      boost::shared_ptr< linearOperator<ARRAY_TYPE_ELEMENT> > op =
		(p==0 ) ? this->regularization_group_operators_.at(i).at(j) : 
		this->regularization_group_prior_operators_.at(group_idx+j);	      
	      sums[idx].create( image_dims.get() );
	      reg_out[idx].create( image_dims.get() );
	      
	      sums[idx] = *b_k[operator_idx+idx];
	      
	      op->mult_M( u_k.get(), &reg_out[idx] );
	      
	      if( p==1 ){
	      	reg_out[idx] -= *p_M[operator_idx/2+j];
	      }	  
	      reg_out[idx] += sums[idx];
	      
	      boost::shared_ptr<ARRAY_TYPE_REAL> tmp_s_k = abs( &sums[idx] );
	      s_k[p] += *tmp_s_k;
	    }
	    s_k[p].sqrt();
	  }

	  for( unsigned int p=0; p<((get_prior_image().get()) ? 2 : 1); p++ ){
	    for( unsigned int j=0; j<k; j++ ){
	      
	      unsigned int idx = p*k+j;

	      boost::shared_ptr< linearOperator<ARRAY_TYPE_ELEMENT> > op =
		(p==0 ) ? this->regularization_group_operators_.at(i).at(j) : 
		this->regularization_group_prior_operators_.at(group_idx+j);
	      
	      // Update of d_k
	      shrinkd( REAL(1)/op->get_weight(), &s_k[p], &sums[idx], d_k[operator_idx+idx].get() );
	      
	      // Update of b_k (only in the last inner iteration)
	      if( inner_iteration == inner_iterations-1 ){
	      	reg_out[idx] -= *d_k[operator_idx+idx];
	      	*b_k[operator_idx+idx] += reg_out[idx];
	      }
	    }
	  }
	  operator_idx += (k*((get_prior_image().get()) ? 2 : 1)); group_idx += k;
	  delete[] sums; delete[] reg_out; delete[] s_k;
	}

      } // end of inner loop

      // Output change in u_k
      if( tolerance > REAL(0) || this->output_mode_ >= solver<ARRAY_TYPE_ELEMENT, ARRAY_TYPE_ELEMENT>::OUTPUT_VERBOSE ){
      	u_k_prev *= ELEMENT_TYPE(-1);
      	u_k_prev += *u_k;

      	REAL delta = asum(&u_k_prev);

				if( this->output_mode_ >= solver<ARRAY_TYPE_ELEMENT, ARRAY_TYPE_ELEMENT>::OUTPUT_VERBOSE )
					std::cout << "u_k delta (outer loop): " << delta << std::endl << std::endl;

				if( delta < tolerance )
					break;

				u_k_prev = *u_k;
      }

    } // end of outer loop


  }

protected:
  REAL tolerance_;
  unsigned int outer_iterations_, inner_iterations_;
  unsigned int num_reg_operators_;
  std::vector< boost::shared_ptr< linearOperator<ARRAY_TYPE_ELEMENT> > > _regularization_group_operators_;
  std::vector< std::vector< boost::shared_ptr< linearOperator<ARRAY_TYPE_ELEMENT> > > > regularization_group_operators_;
  std::vector< boost::shared_ptr< linearOperator<ARRAY_TYPE_ELEMENT> > > regularization_prior_operators_;
  std::vector< boost::shared_ptr< linearOperator<ARRAY_TYPE_ELEMENT> > > regularization_group_prior_operators_;
  boost::shared_ptr<ARRAY_TYPE_ELEMENT> prior_;
  REAL alpha_;
  boost::shared_ptr<INNER_SOLVER> inner_solver_;
  boost::shared_ptr<encodingOperatorContainer<ARRAY_TYPE_ELEMENT> > enc_op_container_;
  std::vector<unsigned int> weights_backup_;

};
}
