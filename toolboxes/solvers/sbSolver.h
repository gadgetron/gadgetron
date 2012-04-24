/*
  An implementation of the "Generalized Split Bregman Algorithm" - sec. 3.2. of the paper
  "The Split Bregman Method for L1-Regularized Problems" by Tom Goldstein and Stanley Osher. 
  Siam J. Imaging Sciences. Vol. 2, No. 2, pp. 323-343.
*/

#pragma once

#include "linearSolver.h"
#include "encodingOperatorContainer.h"
#include "vector_td_utilities.h"

#include <vector>
#include <iostream>

template <class REAL, class ELEMENT_TYPE, class ARRAY_TYPE_REAL, class ARRAY_TYPE_ELEMENT, class INNER_SOLVER> 
class sbSolver : public linearSolver<REAL, ELEMENT_TYPE, ARRAY_TYPE_ELEMENT>
{
public:

  // Constructor
  sbSolver() : linearSolver<REAL, ELEMENT_TYPE, ARRAY_TYPE_ELEMENT>() 
  { 
    tolerance_ = REAL(0);
    outer_iterations_ = 10;
    inner_iterations_ = 1;
    num_reg_operators_ = 0;
    inner_solver_ = boost::shared_ptr<INNER_SOLVER>( new INNER_SOLVER() );
    enc_op_container_ = boost::shared_ptr< encodingOperatorContainer<REAL, ARRAY_TYPE_ELEMENT> >
      ( new encodingOperatorContainer<REAL, ARRAY_TYPE_ELEMENT>() );
    inner_solver_->set_encoding_operator( enc_op_container_ );
  }
  
  // Destructor
  virtual ~sbSolver() {}

  // Set encoding operator (only one allowed)
  virtual bool set_encoding_operator( boost::shared_ptr< linearOperator<REAL, ARRAY_TYPE_ELEMENT> > op )
  {
    if( !validate_encoding_operator( op ) ){
      this->solver_error( "Error: sbSolver::add_encoding_operator : operator validation failed" );
      return false;
    }
    
    if( regularization_operators_.size() || _regularization_group_operators_.size() || regularization_group_operators_.size() ){
      this->solver_error( "Error: sbSolver::add_encoding_operator : the encoding operator must be set prior to adding any regularization operators" );
      return false;
    }
    
    if( !enc_op_container_->add_operator( op ) ){
      this->solver_error( "Error: sbSolver::add_encoding_operator : failed to add encoding operator" );
      return false;
    }
    
    // The domain is constant for all operators (the image dimensions)
    enc_op_container_->set_domain_dimensions( op->get_domain_dimensions().get() );

    encoding_operator_ = op;
    return true;
  }

  // Add regularization operator (anisotropic, multiple allowed)
  virtual bool add_regularization_operator( boost::shared_ptr< linearOperator<REAL, ARRAY_TYPE_ELEMENT> > op )
  {
    if( !validate_regularization_operator( op ) ){
      this->solver_error( "Error: sbSolver::add_regularization_operator : operator validation failed" );
      return false;
    }
    
    if( !enc_op_container_->add_operator( op ) ){
      this->solver_error( "Error: sbSolver::add_regularization_operator : failed to add regularization operator" );
      return false;
    }
    
    regularization_operators_.push_back(op);
    
    return true;
  }

  // Add regularization operator (isotropic, multiple operators per group allowed)
  virtual bool add_regularization_group_operator ( boost::shared_ptr< linearOperator<REAL, ARRAY_TYPE_ELEMENT> > op ) 
  {
    if( !validate_regularization_operator( op ) ){
      this->solver_error( "Error: sbSolver::add_regularization_group_operator : encoding operator must be set before adding regularization operators" );
      return false;
    }
    
    _regularization_group_operators_.push_back( op );

    return true;
  }
  
  // Add isotroic regularization group (multiple groups allowed)
  virtual bool add_group()
  {
    if( _regularization_group_operators_.size() == 0 ){
      this->solver_error( "Error: sbSolver::add_group : no regularization group operators added" );
      return false;
    }
    
    for( unsigned int i=0; i<_regularization_group_operators_.size(); i++ ){
      if( !enc_op_container_->add_operator( _regularization_group_operators_.at(i) ) ){
	this->solver_error( "Error: sbSolver::add_group : failed to add regularization operator" );
	return false;
      }
    }
    
    regularization_group_operators_.push_back( _regularization_group_operators_ );    
    _regularization_group_operators_.clear();
    
    return true;
  }

  // Set termination criterium tolerance
  virtual void set_tc_tolerance( REAL tolerance ) 
  {
    if( tolerance < REAL(0) ) 
      this->solver_warning( "Warning: sbSolver::set_tc_tolerence : tolerance cannot be negative. Ignored." );
    else tolerance_ = tolerance;
  }

  // Set/get maximum number of outer Split-Bregman iterations
  virtual void set_max_outer_iterations( unsigned int iterations ) { outer_iterations_ = iterations; }
  virtual unsigned int get_max_outer_iterations() { return outer_iterations_; }

  // Set/get maximum number of inner Split-Bregman iterations
  virtual void set_max_inner_iterations( unsigned int iterations ) { inner_iterations_ = iterations; }
  virtual unsigned int get_max_inner_iterations() { return inner_iterations_; }

  // Get the inner solver
  virtual boost::shared_ptr<INNER_SOLVER> get_inner_solver() { return inner_solver_; }
  
  // Core solver functionality to be implemented in a derived class (host/device specific implementations)
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

  // Main solver interface
  virtual boost::shared_ptr<ARRAY_TYPE_ELEMENT> solve( ARRAY_TYPE_ELEMENT *_f )
  {
    // Check if everything is set up right
    //
    if( !validate_solve() ){
      this->solver_error( "Error: sbSolver::solve : solver validation failed");
      return boost::shared_ptr<ARRAY_TYPE_ELEMENT>();
    }
    
    // Define u_k
    //
    boost::shared_ptr<ARRAY_TYPE_ELEMENT> u_k( new ARRAY_TYPE_ELEMENT() );
    
    if( !u_k->create( encoding_operator_->get_domain_dimensions().get() )){
      this->solver_error( "Error: sbSolver::solve : memory allocation of u_k failed" );
      return boost::shared_ptr<ARRAY_TYPE_ELEMENT>();
    }        

    // Use x0 (if provided) as starting estimate
    if( this->get_x0().get() )
      *u_k = *(this->get_x0());
    else if( !solver_clear_element( u_k.get() )){
      this->solver_error( "Error: sbSolver::solve : failed to clear u_k" );
      return boost::shared_ptr<ARRAY_TYPE_ELEMENT>();
    }
    get_inner_solver()->set_x0( u_k );

    // Initialze d and b
    //
    boost::shared_array< boost::shared_ptr<ARRAY_TYPE_ELEMENT> > d_k, b_k;
    
    if( !initialize( d_k, b_k ) ){
      this->solver_error( "Error: sbSolver::solve : solver initialization failed");
      return boost::shared_ptr<ARRAY_TYPE_ELEMENT>();
    }

    // Normalize (a copy of) the input data
    //

    boost::shared_ptr<ARRAY_TYPE_ELEMENT> f(new ARRAY_TYPE_ELEMENT(*_f));

    if( !f || !f->get_data_ptr() ){
      this->solver_error( "Error: sbSolver::solve : memory allocation of f failed" );
      return boost::shared_ptr<ARRAY_TYPE_ELEMENT>();
    }

    ELEMENT_TYPE image_scale;
    if( !normalize( f, image_scale ) ){
      this->solver_error( "Error: sbSolver::solve : normalization failed" );
      return boost::shared_ptr<ARRAY_TYPE_ELEMENT>();
    }

    // Invoke the core solver
    //
    if( !core( tolerance_, outer_iterations_, inner_iterations_, f, u_k, d_k, b_k ) ){
      this->solver_error( "Error: sbSolver::solve : core solver failed" );
      return boost::shared_ptr<ARRAY_TYPE_ELEMENT>();
    } 

    // Undo the intermediate scaling of u_k ...
    //
    if( !undo_normalization( u_k, 1/image_scale )){
      this->solver_error( "Error: sbSolver::solve : unable to undo normalization" );
      return boost::shared_ptr<ARRAY_TYPE_ELEMENT>();
    } 

    // ... and return the result
    //    
    return u_k;
  }

protected:

  // Validate operator
  virtual bool validate_encoding_operator( boost::shared_ptr< linearOperator<REAL, ARRAY_TYPE_ELEMENT> > op )
  {
    if( !op.get() ){
      this->solver_error( "Error: sbSolver::validate_encoding_operator : NULL operator provided" );
      return false;
    }
    
    if( op->get_domain_dimensions()->size() == 0 ){
      this->solver_error( "Error: sbSolver::validate_encoding_operator : operator must have specified domain dimensions" );
      return false;
    }

    if( op->get_codomain_dimensions()->size() == 0 ){
      this->solver_error( "Error: sbSolver::validate_encoding_operator : operator must have specified codomain dimensions" );
      return false;
    }

    return true;
  }

  // Validate regularization operator
  virtual bool validate_regularization_operator( boost::shared_ptr< linearOperator<REAL, ARRAY_TYPE_ELEMENT> > op )
  {
    if( !op.get() ){
      this->solver_error( "Error: sbSolver::validate_regularization_operator : NULL operator provided" );
      return false;
    }
    
    if( !encoding_operator_.get() ){
      this->solver_error( "Error: sbSolver::validate_regularization_operator : encoding operator must be set before adding regularizers" );
      return false;
    }
    
    boost::shared_ptr< std::vector<unsigned int> > image_dims = encoding_operator_->get_domain_dimensions();
    boost::shared_ptr< std::vector<unsigned int> > domain_dims = op->get_domain_dimensions();
    boost::shared_ptr< std::vector<unsigned int> > codomain_dims = op->get_codomain_dimensions();

    if( domain_dims->size() > 0 ){
      if( *domain_dims != *image_dims ){
	this->solver_error( "Error: sbSolver::validate_regularization_operator : operator domain dimensions mismatch the encoding operator" );
	return false;
      }
    }
    else{
      this->solver_warning( "Warning: sbSolver::validate_regularization_operator : input operator modified -- domain dimensions set" );
      op->set_domain_dimensions( image_dims.get() );
    }
    
    if( codomain_dims->size() > 0 ){
      if( *codomain_dims != *image_dims ){
	this->solver_error( "Error: sbSolver::validate_regularization_operator : operator codomain dimensions mismatch the image dimensions" );
	return false;
      }
    }
    else{
      this->solver_warning( "Warning: sbSolver::validate_regularization_operator : input operator modified -- codomain dimensions set" );
      op->set_codomain_dimensions( image_dims.get() );
    }

    return true;
  }
  
  // Check that the solver is set up properly
  virtual bool validate_solve()
  {
    // Some tests to see if we are ready to go...
    //

    if( !encoding_operator_.get() ){
      this->solver_error( "Error: sbSolver::validate : encoding operator has not been set" );
      return false;
    }

    if( regularization_operators_.size() == 0 && regularization_group_operators_.size() == 0 ){
      this->solver_error( "Error: sbSolver::validate : at least one matrix sparsifying regularizer must be added" );
      return false;
    }

    return true;
  }

  // Initialize d_k and b_k arrays to zero images
  //
  virtual bool initialize( boost::shared_array< boost::shared_ptr<ARRAY_TYPE_ELEMENT> > &d_k, 
			   boost::shared_array< boost::shared_ptr<ARRAY_TYPE_ELEMENT> > &b_k )
  {
    // Determine length of arrays
    num_reg_operators_ = regularization_operators_.size();
    for( unsigned int i=0; i<regularization_group_operators_.size(); i++ ){
      num_reg_operators_ += regularization_group_operators_.at(i).size();
    }

    d_k = boost::shared_array< boost::shared_ptr<ARRAY_TYPE_ELEMENT> > ( new boost::shared_ptr<ARRAY_TYPE_ELEMENT>[num_reg_operators_] );
    b_k = boost::shared_array< boost::shared_ptr<ARRAY_TYPE_ELEMENT> > ( new boost::shared_ptr<ARRAY_TYPE_ELEMENT>[num_reg_operators_] );

    if( !d_k.get() || !b_k.get() ){
      this->solver_error( "Error: sbSolver::initialize : memory allocation of d_k or b_k failed" );
      return false;
    }

    boost::shared_ptr< std::vector<unsigned int> > image_dims = encoding_operator_->get_domain_dimensions();

    for( unsigned int i=0; i<num_reg_operators_; i++ ){

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

    return true;
  }

  // Normalize data
  //
  virtual bool normalize( boost::shared_ptr<ARRAY_TYPE_ELEMENT> f, 
			  ELEMENT_TYPE &image_scale )    
  {    
    // Normalize in image space
    //

    ARRAY_TYPE_ELEMENT tmp;
    if( !tmp.create( encoding_operator_->get_domain_dimensions().get() )){
      this->solver_error( "Error: sbSolver::normalize : memory allocation failed" );
      return false;
    }    
    
    if( encoding_operator_->mult_MH( f.get(), &tmp ) < 0 ){
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
    // Undo normalization of u_k
    //
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
		     boost::shared_array< boost::shared_ptr<ARRAY_TYPE_ELEMENT> > b_k )
  {
    // Image space dimensions
    //

    boost::shared_ptr< std::vector<unsigned int> > image_dims = encoding_operator_->get_domain_dimensions();
    
    // Keep a copy of the "previous" u_k to compute the outer loop change of u_k
    // 

    ARRAY_TYPE_ELEMENT u_k_prev;
    if( tolerance > REAL(0) || this->output_mode_ >= solver<ARRAY_TYPE_ELEMENT, ARRAY_TYPE_ELEMENT>::OUTPUT_VERBOSE ){
      u_k_prev = *u_k;
      if( !u_k_prev.get_data_ptr() ){
	this->solver_error( "Error: sbSolver::core : memory allocation of u_k_prev failed" );
	return false;
      }
    }
    
    // Outer loop
    //

    for( unsigned int outer_iteration=0; outer_iteration<outer_iterations; outer_iteration++ ) {

      if( this->output_mode_ >= solver<ARRAY_TYPE_ELEMENT, ARRAY_TYPE_ELEMENT>::OUTPUT_MAX )
	std::cout << std::endl << "SB outer loop iteration " << outer_iteration << std::endl << std::endl;

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
	  
	  for( unsigned int i=0; i<num_reg_operators_; i++ ){
	    
	    if( !tmp.create( image_dims.get(), data.get_data_ptr()+enc_op_container_->get_offset(i+1) )){
	      this->solver_error( "Error: sbSolver::core : memory allocation for container operator data failed (3)" );
	      return false;
	    }
	    
	    tmp = *d_k[i];
	    
	    if( !solver_axpy_element( ELEMENT_TYPE(-1), b_k[i].get(), &tmp )){
	      this->solver_error( "Error: sbSolver::core : axpy on difference image failed" );
	      return false;
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

	for( unsigned int i=0; i<regularization_operators_.size(); i++ ){

	  ARRAY_TYPE_ELEMENT tmp_sum, reg_out;
	  if( !tmp_sum.create( image_dims.get() ) || !reg_out.create( image_dims.get() ) ){
	    this->solver_error( "Error: sbSolver::core : memory allocation for regularization operator failed in {d_k,b_k} update" );
	    return false;
	  }

	  tmp_sum = *b_k[operator_idx];

	  if( regularization_operators_.at(i)->mult_M( u_k.get(), &reg_out ) < 0 ){
	    this->solver_error( "Error: sbSolver::core : application of regularization operator failed in {d_k,b_k} update" );
	    return false;
	  }

	  if( !solver_axpy_element( ELEMENT_TYPE(1), &reg_out, &tmp_sum )){
	    this->solver_error( "Error: sbSolver::core : computation of shrinkage_1 argument for d_k failed" );
	    return false;
	  }

	  // Update of d_k
	  if( !solver_shrink1( 1/regularization_operators_.at(i)->get_weight(), &tmp_sum, d_k[operator_idx].get() )){
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

	for( unsigned int i=0; i<regularization_group_operators_.size(); i++ ){

	  unsigned int k = regularization_group_operators_.at(i).size();
	  ARRAY_TYPE_ELEMENT *sums = new ARRAY_TYPE_ELEMENT[k], *reg_out = new ARRAY_TYPE_ELEMENT[k];

	  if( !sums || !reg_out ){
	    this->solver_error( "Error: sbSolver::core : host memory allocation for temporary arrays failed" );
	    return false;
	  }

	  ARRAY_TYPE_REAL s_k;
	  if( s_k.create(image_dims.get()) < 0 ){
	    this->solver_error( "Error: sbSolver::core : memory allocation for s_k failed" );
	    return false;
	  }

	  if( !solver_clear_real(&s_k) ){
	    this->solver_error( "Error: sbSolver::core : failed to clear s_k" );
	    return false;
	  } 	  

	  for( unsigned int j=0; j<k; j++ ){

	    if( sums[j].create( image_dims.get() ) < 0 || reg_out[j].create( image_dims.get() ) < 0 ){
	      this->solver_error( "Error: sbSolver::core : memory allocation for regularization operator failed in {d_k,b_k} update" );
	      return false;
	    }

	    sums[j] = *b_k[operator_idx+j];

	    if( regularization_group_operators_.at(i).at(j)->mult_M( u_k.get(), &reg_out[j] ) < 0 ){
	      this->solver_error( "Error: sbSolver::core : application of regularization operator failed in {d_k,b_k} update" );
	      return false;
	    }

	    if( !solver_axpy_element( ELEMENT_TYPE(1), &reg_out[j], &sums[j] )){
	      this->solver_error( "Error: sbSolver::core : computation of shrinkage_d argument for d_k failed" );
	      return false;
	    }

	    boost::shared_ptr<ARRAY_TYPE_REAL> tmp_s_k = solver_norm( &sums[j] );
	    if( !solver_axpy_real( REAL(1), tmp_s_k.get(), &s_k )){
	      this->solver_error( "Error: sbSolver::core : accumulation of s_k failed" );
	      return false;
	    }
	  }

	  if( !solver_sqrt( &s_k )){
	    this->solver_error( "Error: sbSolver::core : sqrt of s_k failed" );
	    return false;
	  }

	  for( unsigned int j=0; j<k; j++ ){

	    // Update of d_k
	    if( !solver_shrinkd( 1/(regularization_group_operators_.at(i).at(j)->get_weight()), &s_k, &sums[j], d_k[operator_idx+j].get() )){
	      this->solver_error( "Error: sbSolver::core : shrinkage_d of d_k failed" );
	      return false;
	    }

	    // Update of b_k (only in the last inner iteration)
	    if( inner_iteration == inner_iterations-1 ){

	      if( !solver_axpy_element( ELEMENT_TYPE(-1), d_k[operator_idx+j].get(), &reg_out[j] )){
		this->solver_error( "Error: sbSolver::core : computation of update argument to b_k failed" );
		return false;
	      }

	      if( !solver_axpy_element( ELEMENT_TYPE(1), &reg_out[j], b_k[operator_idx+j].get() )){
		this->solver_error( "Error: sbSolver::core : update of b_k failed" );
		return false;
	      }
	    }
	  }
	  operator_idx += k;
	  delete[] sums; delete[] reg_out;
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
  boost::shared_ptr< linearOperator<REAL, ARRAY_TYPE_ELEMENT> > encoding_operator_;
  std::vector< boost::shared_ptr< linearOperator<REAL, ARRAY_TYPE_ELEMENT> > > regularization_operators_;
  std::vector< boost::shared_ptr< linearOperator<REAL, ARRAY_TYPE_ELEMENT> > > _regularization_group_operators_;
  std::vector< std::vector< boost::shared_ptr< linearOperator<REAL, ARRAY_TYPE_ELEMENT> > > > regularization_group_operators_;
  boost::shared_ptr<INNER_SOLVER> inner_solver_;
  boost::shared_ptr<encodingOperatorContainer<REAL, ARRAY_TYPE_ELEMENT> > enc_op_container_;
};
