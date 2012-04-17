/*
  An implementation of the "Generalized Split Bregman Algorithm" - sec. 3.2. of the paper
  "The Split Bregman Method for L1-Regularized Problems" by Tom Goldstein and Stanley Osher. 
  Siam J. Imaging Sciences. Vol. 2, No. 2, pp. 323-343.
*/

#pragma once

#include "solver.h"
#include "linearOperator.h"
#include "vector_td_utilities.h"

#include <vector>
#include <iostream>

template <class REAL, class ELEMENT_TYPE, class ARRAY_TYPE_REAL, class ARRAY_TYPE_ELEMENT> class sbSolver 
  : public solver<ARRAY_TYPE_ELEMENT, ARRAY_TYPE_ELEMENT>
{
public:

  sbSolver() : solver<ARRAY_TYPE_ELEMENT, ARRAY_TYPE_ELEMENT>() 
  { 
    tolerance_ = REAL(0);
    outer_iterations_ = 10;
    inner_iterations_ = 1;
  }

  virtual ~sbSolver() {}

  virtual int set_inner_solver( boost::shared_ptr< solver<ARRAY_TYPE_ELEMENT,ARRAY_TYPE_ELEMENT> > solver ) {
    inner_solver_ = solver;
    return 0;
  }

  virtual bool set_encoding_operator( boost::shared_ptr< linearOperator<REAL, ARRAY_TYPE_ELEMENT> > op ) {
    encoding_operator_ = op;
    return true;
  }

  virtual bool add_regularization_operator( boost::shared_ptr< linearOperator<REAL, ARRAY_TYPE_ELEMENT> > op, ARRAY_TYPE_ELEMENT *prior = 0x0  )
  {
    regularization_operators_.push_back(op);
    boost::shared_ptr<ARRAY_TYPE_ELEMENT> opp;
    if( prior ){
      opp = boost::shared_ptr<ARRAY_TYPE_ELEMENT>( new ARRAY_TYPE_ELEMENT() );
      if( !opp.get() || opp->create( prior->get_dimensions().get() ) < 0 ){
	this->solver_error( "sbSolver::add_regularization_operator : allocation failed for prior" );
	regularization_priors_.push_back( boost::shared_ptr<ARRAY_TYPE_ELEMENT>() );
	return false;
      }
      if( op->mult_M( prior, opp.get() ) < 0 ){
	this->solver_error( "sbSolver::add_regularization_operator : could not apply operator to prior" );
	regularization_priors_.push_back( boost::shared_ptr<ARRAY_TYPE_ELEMENT>() );
	return false;
      }      
    }
    regularization_priors_.push_back( opp );
    return true;
  }

  virtual bool add_regularization_group_operator( boost::shared_ptr< linearOperator<REAL, ARRAY_TYPE_ELEMENT> > op ) {
    _regularization_group_operators_.push_back(op);
    return true;
  }

  virtual bool add_group( ARRAY_TYPE_ELEMENT *prior = 0x0 )
  {
    regularization_group_operators_.push_back(_regularization_group_operators_);
    _regularization_group_operators_.clear();

    boost::shared_ptr<ARRAY_TYPE_ELEMENT> opp;
    std::vector< boost::shared_ptr<ARRAY_TYPE_ELEMENT> > _regularization_group_priors_;

    if( prior ){
      for( unsigned int i=0; i<regularization_group_operators_.back().size(); i++ ){

	opp = boost::shared_ptr<ARRAY_TYPE_ELEMENT>( new ARRAY_TYPE_ELEMENT() );

	if( !opp.get() || opp->create( prior->get_dimensions().get() ) < 0 ){
	  this->solver_error( "sbSolver::add_group : allocation failed for prior" );
	  regularization_group_priors_.push_back( std::vector< boost::shared_ptr< ARRAY_TYPE_ELEMENT > >() );
	  return false;
	}

	if( regularization_group_operators_.back().at(i)->mult_M( prior, opp.get() ) < 0 ){
	  this->solver_error( "sbSolver::add_group : could not apply operator to prior" );
	  regularization_group_priors_.push_back( std::vector< boost::shared_ptr< ARRAY_TYPE_ELEMENT > >() );
	  return false;
	}      

	_regularization_group_priors_.push_back( opp );
      }
    }
    regularization_group_priors_.push_back( _regularization_group_priors_ );

    return true;
  }

  virtual void set_tolerance( REAL tolerance ) {
    if( tolerance < REAL(0) ) this->solver_error( "sbSolver::set_tolerence : tolerance cannot be negative" );
    else tolerance_ = tolerance;
  }

  virtual void set_outer_iterations( unsigned int iterations ) {
    outer_iterations_ = iterations;
  }

  virtual void set_inner_iterations( unsigned int iterations ) {
    inner_iterations_ = iterations;
  }

  virtual void set_image_dimensions( boost::shared_ptr< std::vector<unsigned int> > dims ){
    image_dims_ = dims;
  }

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
    if( !validate() ){
      this->solver_error( "sbSolver::solve : setup failed validation");
      return boost::shared_ptr<ARRAY_TYPE_ELEMENT>();
    }

    // Initialze d and b
    //
    boost::shared_array< boost::shared_ptr<ARRAY_TYPE_ELEMENT> > d_k;
    boost::shared_array< boost::shared_ptr<ARRAY_TYPE_ELEMENT> > b_k;
    if( !initialize( d_k, b_k ) ){
      this->solver_error( "sbSolver::solve : failed to initialize");
      return boost::shared_ptr<ARRAY_TYPE_ELEMENT>();
    }

    // Make a copy of _f before normalization
    //
    boost::shared_ptr<ARRAY_TYPE_ELEMENT> f( new ARRAY_TYPE_ELEMENT(*_f) );
    if( !f->get_data_ptr() ){
      this->solver_error( "sbSolver::solve : memory allocation of f failed" );
      return boost::shared_ptr<ARRAY_TYPE_ELEMENT>();
    }

    // Define u_k
    //
    boost::shared_ptr<ARRAY_TYPE_ELEMENT> u_k;

    // Normalize the data
    //
    ELEMENT_TYPE image_scale = ELEMENT_TYPE(0);
    if( !normalize( f, u_k, image_scale ) ){
      this->solver_error( "sbSolver::solve : normalization failed" );
      return boost::shared_ptr<ARRAY_TYPE_ELEMENT>();
    }

    // Keep a copy of E^H f for subsequent rhs computations in the inner solver
    // 
    boost::shared_ptr<ARRAY_TYPE_ELEMENT> muEHf( new ARRAY_TYPE_ELEMENT(*u_k) );
    if( !muEHf->get_data_ptr() ){
      this->solver_error( "sbSolver::solve : memory allocation of muEHf failed" );
      return boost::shared_ptr<ARRAY_TYPE_ELEMENT>();
    }

    // Scale EHf with mu
    //
    if( !solver_scal( REAL(encoding_operator_->get_weight()), muEHf.get() ) ){
      this->solver_error( "sbSolver::solve : error scaling EHf with mu" );
      return boost::shared_ptr<ARRAY_TYPE_ELEMENT>();
    }

    // Invoke the core solver
    //
    if( !core( tolerance_, outer_iterations_, inner_iterations_, f, muEHf, u_k, d_k, b_k ) ){
      this->solver_error( "sbSolver::solve : core solver failed" );
      return boost::shared_ptr<ARRAY_TYPE_ELEMENT>();
    } 

    // Undo the intermediate scaling of u_k ...
    //
    if( !undo_normalization( u_k, 1/image_scale )){
      this->solver_error( "sbSolver::solve : unable to undo normalization" );
      return boost::shared_ptr<ARRAY_TYPE_ELEMENT>();
    } 

    // ... and return the result
    //    
    return u_k;
  }

protected:

  // Check that the solver is set up properly
  virtual bool validate()
  {
    // Some tests to see if we are ready to go...
    //
    if( !inner_solver_.get() ){
      this->solver_error( "sbSolver::validate : inner solver has not been set" );
      return false;
    }

    if( !encoding_operator_.get() ){
      this->solver_error( "sbSolver::validate : encoding operator has not been set" );
      return false;
    }

    if( regularization_operators_.size() == 0 && regularization_group_operators_.size() == 0 ){
      this->solver_error( "sbSolver::validate : at least one matrix regularizer must be added" );
      return false;
    }

    if( !image_dims_.get() ){
      this->solver_error( "sbSolver::validate : image dimensions have not been set" );
      return false;
    }

    for( unsigned i=0; i<regularization_priors_.size(); i++ ){
      if( regularization_priors_.at(i) && !regularization_priors_.at(i)->dimensions_equal( image_dims_.get()) ){
	this->solver_error( "sbSolver::validate : Regularization prior does not match specified image dimensions" );
	return false;
      }
    }

    for( unsigned i=0; i<regularization_group_priors_.size(); i++ ){
      for( unsigned j=0; j<regularization_group_priors_.at(i).size(); j++ ){
	if( regularization_group_priors_.at(i).size() > 0 && !regularization_group_priors_.at(i).at(j)->dimensions_equal( image_dims_.get()) ){
	  this->solver_error( "sbSolver::validate : Regularization group prior does not match specified image dimensions" );
	  return false;
	}
      }
    }
    return true;
  }

  // Initialize d_k and b_k arrays to zero images
  //
  virtual bool initialize( boost::shared_array< boost::shared_ptr<ARRAY_TYPE_ELEMENT> > &d_k, 
			   boost::shared_array< boost::shared_ptr<ARRAY_TYPE_ELEMENT> > &b_k )
  {
    // Determine length of arrays
    unsigned int num_reg_operators = regularization_operators_.size();
    for( unsigned int i=0; i<regularization_group_operators_.size(); i++ )
      num_reg_operators += regularization_group_operators_.at(i).size();

    d_k = boost::shared_array< boost::shared_ptr<ARRAY_TYPE_ELEMENT> > ( new boost::shared_ptr<ARRAY_TYPE_ELEMENT>[num_reg_operators] );
    b_k = boost::shared_array< boost::shared_ptr<ARRAY_TYPE_ELEMENT> > ( new boost::shared_ptr<ARRAY_TYPE_ELEMENT>[num_reg_operators] );

    if( !d_k.get() || !b_k.get() ){
      this->solver_error( "sbSolver::initialize : memory allocation of d_k or b_k failed" );
      return false;
    }

    for( unsigned int i=0; i<num_reg_operators; i++ ){

      d_k[i] = boost::shared_ptr<ARRAY_TYPE_ELEMENT>(new ARRAY_TYPE_ELEMENT());

      if( d_k[i] ) d_k[i]->create( image_dims_.get() );

      if( !d_k[i]->get_data_ptr() ){
	this->solver_error( "sbSolver::initialize : memory allocation of d_k failed" );
	return false;
      }

      if( !solver_clear_element( d_k[i].get() )){
	this->solver_error( "sbSolver::initialize : failed to clear internal memory buffer d_k" );
	return false;
      }

      b_k[i] = boost::shared_ptr<ARRAY_TYPE_ELEMENT>(new ARRAY_TYPE_ELEMENT());

      if( b_k[i] ) b_k[i]->create( image_dims_.get() );

      if( !b_k[i]->get_data_ptr() ){
	this->solver_error( "sbSolver::initialize : memory allocation of b_k failed" );
	return false;
      }

      if( !solver_clear_element( b_k[i].get() )){
	this->solver_error( "sbSolver::initialize : failed to clear internal memory buffer b_k" );
	return false;
      } 
    }
    return true;
  }

  // Normalize data
  //
  virtual bool normalize( boost::shared_ptr<ARRAY_TYPE_ELEMENT> &f, 
			  boost::shared_ptr<ARRAY_TYPE_ELEMENT> &u_k, 
			  ELEMENT_TYPE &image_scale )    
  {    
    // Initialize u_k to E^H f 
    //
    u_k = boost::shared_ptr<ARRAY_TYPE_ELEMENT>(new ARRAY_TYPE_ELEMENT());
    if( u_k.get() ) u_k->create( image_dims_.get() );

    if( !u_k.get() || !u_k->get_data_ptr() ){
      this->solver_error( "sbSolver::normalize : memory allocation of u_k failed" );
      return false;
    }    

    if( encoding_operator_->mult_MH( f.get(), u_k.get() ) < 0 ){
      this->solver_error( "sbSolver::normalize : adjoint encoding operation failed on f" );
      return false;
    }

    // Normalize to an average energy of "one intensity unit per image element"
    //
    REAL sum = solver_asum( u_k.get() );
    image_scale =  (REAL) (u_k->get_number_of_elements())/REAL(sum);

    // Normalize u_k and f
    //
    if( !solver_scal( image_scale, u_k.get() )){
      this->solver_error( "sbSolver::normalize : unable to scale u_k" );
      return false;
    }

    if(	!solver_scal( image_scale, f.get() )){
      this->solver_error( "sbSolver::normalize : unable to scale f" );
      return false;
    }

    // Normalize regularization priors if available
    //
    for( unsigned i=0; i<regularization_priors_.size(); i++ ){
      if( regularization_priors_.at(i) ){
	if( !solver_scal( image_scale, regularization_priors_.at(i).get() )){
	  this->solver_error( "sbSolver::normalize : unable to scale prior" );
	  return false;
	}
      }
    }

    for( unsigned i=0; i<regularization_group_priors_.size(); i++ ){
      for( unsigned j=0; j<regularization_group_priors_.at(i).size(); j++ ){
	if( !solver_scal( image_scale, regularization_group_priors_.at(i).at(j).get() )){
	  this->solver_error( "sbSolver::normalize : unable to scale group prior" );
	  return false;
	}
      }
    }

    return true;
  }

  // Undo normalization
  //
  virtual bool undo_normalization( boost::shared_ptr<ARRAY_TYPE_ELEMENT> &u_k, ELEMENT_TYPE undo_scale )
  {
    // Undo normalization of u_k
    //
    if(	!solver_scal( undo_scale, u_k.get() )){
      this->solver_error( "sbSolver::undo_normalization : unable to undo scaling of u_k" );
      return false;
    }

    // Undo normalization of regularization priors
    //
    for( unsigned i=0; i<regularization_priors_.size(); i++ ){
      if( regularization_priors_.at(i) ){
	if( !solver_scal( undo_scale, regularization_priors_.at(i).get() )){
	  this->solver_error( "sbSolver::undo_normalization : unable to undo scaling of prior" );
	  return false;
	}
      }
    }

    for( unsigned i=0; i<regularization_group_priors_.size(); i++ ){
      for( unsigned j=0; j<regularization_group_priors_.at(i).size(); j++ ){
	if( !solver_scal( undo_scale, regularization_group_priors_.at(i).at(j).get() )){
	  this->solver_error( "sbSolver::normalization : unable to undo scaling of group prior" );
	  return false;
	}
      }
    }

    return true;
  }

  // The core of the Split Bregman solver.
  //
  virtual bool core( REAL tolerance, unsigned int outer_iterations, unsigned int inner_iterations,
		     boost::shared_ptr<ARRAY_TYPE_ELEMENT> &f, 
		     boost::shared_ptr<ARRAY_TYPE_ELEMENT> &muEHf,
		     boost::shared_ptr<ARRAY_TYPE_ELEMENT> &u_k,
		     boost::shared_array< boost::shared_ptr<ARRAY_TYPE_ELEMENT> > &d_k,
		     boost::shared_array< boost::shared_ptr<ARRAY_TYPE_ELEMENT> > &b_k )
  {

    // Keep a copy of the "previous" u_k to compute the outer loop change of u_k
    // 
    ARRAY_TYPE_ELEMENT u_k_prev;
    if( tolerance > REAL(0) || this->output_mode_ >= solver<ARRAY_TYPE_ELEMENT, ARRAY_TYPE_ELEMENT>::OUTPUT_VERBOSE ){
      u_k_prev = *u_k;
      if( !u_k_prev.get_data_ptr() ){
	this->solver_error( "sbSolver::core : memory allocation of u_k_prev failed" );
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

	// Form rhs for inner loop solver (initializes to adjoint encoding operator)
	// 
	ARRAY_TYPE_ELEMENT rhs(*muEHf);
	if( !rhs.get_data_ptr() ){
	  this->solver_error( "sbSolver::core : memory allocation of rhs failed" );
	  return false;
	}

	unsigned int operator_idx = 0;

	// Add regularization operators to rhs
	//
	for( unsigned int i=0; i<regularization_operators_.size(); i++ ){

	  ARRAY_TYPE_ELEMENT tmp_diff, reg_out;
	  if( tmp_diff.create( image_dims_.get() ) < 0 || reg_out.create( image_dims_.get() ) < 0 ){
	    this->solver_error( "sbSolver::core : memory allocation for regularization operator failed in rhs computation" );
	    return false;
	  }    

	  tmp_diff = *d_k[operator_idx];

	  if( regularization_priors_.at(i).get() ){
	    if( !solver_axpy_element( ELEMENT_TYPE(1), regularization_priors_.at(i).get(), &tmp_diff )){
	      this->solver_error( "sbSolver::core : could not add regularization prior in rhs computation" );
	      return false;
	    }
	  }

	  if( !solver_axpy_element( ELEMENT_TYPE(0)-ELEMENT_TYPE(1), b_k[operator_idx].get(), &tmp_diff )){
	    this->solver_error( "sbSolver::core : computation of regularization argument failed in rhs computation" );
	    return false;
	  }    

	  if( regularization_operators_.at(i)->mult_MH( &tmp_diff, &reg_out ) < 0 ){
	    this->solver_error( "sbSolver::core : application of regularization operator failed in rhs computation" );
	    return false;
	  }    

	  if( !solver_axpy_element( REAL(regularization_operators_.at(i)->get_weight()), &reg_out, &rhs )){
	    this->solver_error( "sbSolver::core : accumulation in rhs computation failed (1)" );
	    return false;
	  }
	  operator_idx++;
	}	

	// Add regularization group operators to rhs
	//
	for( unsigned int i=0; i<regularization_group_operators_.size(); i++ ){
	  for( unsigned int j=0; j<regularization_group_operators_.at(i).size(); j++ ){

	    ARRAY_TYPE_ELEMENT tmp_diff, reg_out;
	    if( tmp_diff.create( image_dims_.get() ) < 0 || reg_out.create( image_dims_.get() ) < 0 ){
	      this->solver_error( "sbSolver::core : memory allocation for group regularization operator failed in rhs computation" );
	      return false;
	    }    

	    tmp_diff = *d_k[operator_idx];

	    if( regularization_group_priors_.at(i).size() > 0 && regularization_group_priors_.at(i).at(j).get() ){
	      if( !solver_axpy_element( ELEMENT_TYPE(1), regularization_group_priors_.at(i).at(j).get(), &tmp_diff )){
		this->solver_error( "sbSolver::core : could not add regularization group prior in rhs computation" );
		return false;
	      }
	    }

	    if( !solver_axpy_element( ELEMENT_TYPE(0)-ELEMENT_TYPE(1), b_k[operator_idx].get(), &tmp_diff )){
	      this->solver_error( "sbSolver::core : computation of group regularization argument failed in rhs computation" );
	      return false;
	    }    

	    if( regularization_group_operators_.at(i).at(j)->mult_MH( &tmp_diff, &reg_out ) < 0 ){
	      this->solver_error( "sbSolver::core : application of group regularization operator failed in rhs computation" );
	      return false;
	    }    

	    if( !solver_axpy_element( REAL(regularization_group_operators_.at(i).at(j)->get_weight()), &reg_out, &rhs )){
	      this->solver_error( "sbSolver::core : accumulation in rhs computation failed (2)" );
	      return false;
	    }
	    operator_idx++;
	  }
	}	

	// Solve for u_k
	//
	{
	  boost::shared_ptr<ARRAY_TYPE_ELEMENT> tmp;//TODO = inner_solver_->solve(&rhs);

	  // Compute change in u_k
	  //
	  if( this->output_mode_ >= solver<ARRAY_TYPE_ELEMENT, ARRAY_TYPE_ELEMENT>::OUTPUT_VERBOSE ){
	    if( !solver_axpy_element( ELEMENT_TYPE(0)-ELEMENT_TYPE(1), tmp.get(), u_k.get() )){
	      this->solver_error( "sbSolver::core : error computing inner loop u_k delta" );
	      return false;
	    }
	    std::cout << std::endl << "u_k delta (inner loop): " << solver_asum(u_k.get()) << std::endl;
	  }

	  // Update u_k
	  u_k = tmp;
	}

	operator_idx = 0;

	// Update d_k (and in the final inner iteration b_k also)
	//
	for( unsigned int i=0; i<regularization_operators_.size(); i++ ){

	  ARRAY_TYPE_ELEMENT tmp_sum, reg_out;
	  if( tmp_sum.create( image_dims_.get() ) < 0 || reg_out.create( image_dims_.get() ) < 0 ){
	    this->solver_error( "sbSolver::core : memory allocation for regularization operator failed in {d_k,b_k} update" );
	    return false;
	  }

	  tmp_sum = *b_k[operator_idx];

	  if( regularization_operators_.at(i)->mult_M( u_k.get(), &reg_out ) < 0 ){
	    this->solver_error( "sbSolver::core : application of regularization operator failed in {d_k,b_k} update" );
	    return false;
	  }

	  if( regularization_priors_.at(i).get() ) {
	    if( !solver_axpy_element( ELEMENT_TYPE(0)-ELEMENT_TYPE(1), regularization_priors_.at(i).get(), &reg_out )){
	      this->solver_error( "sbSolver::core : application of regularization prior failed in {d_k,b_k} update" );
	      return false;
	    }
	  }

	  if( !solver_axpy_element( ELEMENT_TYPE(1), &reg_out, &tmp_sum )){
	    this->solver_error( "sbSolver::core : computation of shrinkage_1 argument for d_k failed" );
	    return false;
	  }

	  // Update of d_k
	  if( !solver_shrink1( 1/regularization_operators_.at(i)->get_weight(), &tmp_sum, d_k[operator_idx].get() )){
	    this->solver_error( "sbSolver::core : shrinkage_1 of d_k failed" );
	    return false;
	  }

	  // Update of b_k (only in the last inner iteration)
	  if( inner_iteration == inner_iterations-1 ){

	    if( !solver_axpy_element( ELEMENT_TYPE(0)-ELEMENT_TYPE(1), d_k[operator_idx].get(), &reg_out )){
	      this->solver_error( "sbSolver::core : computation of update argument to b_k failed" );
	      return false;
	    }

	    // Update of b_k
	    if( !solver_axpy_element( ELEMENT_TYPE(1), &reg_out, b_k[operator_idx].get() )){
	      this->solver_error( "sbSolver::core : update of b_k failed" );
	      return false;
	    }
	  }
	  operator_idx++;
	}

	for( unsigned int i=0; i<regularization_group_operators_.size(); i++ ){

	  unsigned int k = regularization_group_operators_.at(i).size();
	  ARRAY_TYPE_ELEMENT *sums = new ARRAY_TYPE_ELEMENT[k], *reg_out = new ARRAY_TYPE_ELEMENT[k];

	  if( !sums || !reg_out ){
	    this->solver_error( "sbSolver::core : host memory allocation for temporary arrays failed" );
	    return false;
	  }

	  ARRAY_TYPE_REAL s_k;
	  if( s_k.create(image_dims_.get()) < 0 ){
	    this->solver_error( "sbSolver::core : memory allocation for s_k failed" );
	    return false;
	  }

	  if( !solver_clear_real(&s_k) ){
	    this->solver_error( "sbSolver::core : failed to clear s_k" );
	    return false;
	  } 	  

	  for( unsigned int j=0; j<k; j++ ){

	    if( sums[j].create( image_dims_.get() ) < 0 || reg_out[j].create( image_dims_.get() ) < 0 ){
	      this->solver_error( "sbSolver::core : memory allocation for regularization operator failed in {d_k,b_k} update" );
	      return false;
	    }

	    sums[j] = *b_k[operator_idx+j];

	    if( regularization_group_operators_.at(i).at(j)->mult_M( u_k.get(), &reg_out[j] ) < 0 ){
	      this->solver_error( "sbSolver::core : application of regularization operator failed in {d_k,b_k} update" );
	      return false;
	    }

	    if( regularization_group_priors_.at(i).size() > 0 && regularization_group_priors_.at(i).at(j).get() ) {
	      if( !solver_axpy_element( ELEMENT_TYPE(0)-ELEMENT_TYPE(1), regularization_group_priors_.at(i).at(j).get(), &reg_out[j] )){
		this->solver_error( "sbSolver::core : application of regularization group prior failed in {d_k,b_k} update" );
		return false;
	      }
	    }

	    if( !solver_axpy_element( ELEMENT_TYPE(1), &reg_out[j], &sums[j] )){
	      this->solver_error( "sbSolver::core : computation of shrinkage_d argument for d_k failed" );
	      return false;
	    }

	    boost::shared_ptr<ARRAY_TYPE_REAL> tmp_s_k = solver_norm(&sums[j]);
	    if( !solver_axpy_real( REAL(1), tmp_s_k.get(), &s_k )){
	      this->solver_error( "sbSolver::core : accumulation of s_k failed" );
	      return false;
	    }
	  }

	  if( !solver_sqrt(&s_k) ){
	    this->solver_error( "sbSolver::core : sqrt of s_k failed" );
	    return false;
	  }

	  for( unsigned int j=0; j<k; j++ ){

	    // Update of d_k
	    if( !solver_shrinkd( 1/(regularization_group_operators_.at(i).at(j)->get_weight()), &s_k, &sums[j], d_k[operator_idx+j].get() )){
	      this->solver_error( "sbSolver::core : shrinkage_d of d_k failed" );
	      return false;
	    }

	    // Update of b_k (only in the last inner iteration)
	    if( inner_iteration == inner_iterations-1 ){

	      if( !solver_axpy_element( ELEMENT_TYPE(0)-ELEMENT_TYPE(1), d_k[operator_idx+j].get(), &reg_out[j] )){
		this->solver_error( "sbSolver::core : computation of update argument to b_k failed" );
		return false;
	      }

	      if( !solver_axpy_element( ELEMENT_TYPE(1), &reg_out[j], b_k[operator_idx+j].get() )){
		this->solver_error( "sbSolver::core : update of b_k failed" );
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

	if( !solver_scal( ELEMENT_TYPE(0)-ELEMENT_TYPE(1), &u_k_prev ) ){
	  this->solver_error( "sbSolver::core : error computing inner loop u_k delta (scale)" );
	  return false;
	}

	if( !solver_axpy_element( ELEMENT_TYPE(1), u_k.get(), &u_k_prev )){
	  this->solver_error( "sbSolver::core : error computing inner loop u_k delta (axpy)" );
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
  boost::shared_ptr< std::vector<unsigned int> > image_dims_;
  boost::shared_ptr< solver<ARRAY_TYPE_ELEMENT, ARRAY_TYPE_ELEMENT> > inner_solver_;
  boost::shared_ptr< linearOperator<REAL, ARRAY_TYPE_ELEMENT> > encoding_operator_;
  std::vector< boost::shared_ptr< linearOperator<REAL, ARRAY_TYPE_ELEMENT> > > regularization_operators_;
  std::vector< boost::shared_ptr< linearOperator<REAL, ARRAY_TYPE_ELEMENT> > > _regularization_group_operators_;
  std::vector< std::vector< boost::shared_ptr< linearOperator<REAL, ARRAY_TYPE_ELEMENT> > > > regularization_group_operators_;
  std::vector< boost::shared_ptr<ARRAY_TYPE_ELEMENT> > regularization_priors_;
  std::vector< std::vector< boost::shared_ptr<ARRAY_TYPE_ELEMENT> > > regularization_group_priors_;
};
