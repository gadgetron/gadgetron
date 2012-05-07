#pragma once

#include "linearSolver.h"
#include "cgCallback.h"
#include "cgPreconditioner.h"
#include "real_utilities.h"
#include "cgCallback.h"
#include <vector>
#include <iostream>

template <class REAL, class ELEMENT_TYPE, class ARRAY_TYPE> class cgSolver 
  : public linearSolver<REAL, ELEMENT_TYPE, ARRAY_TYPE>
{
  
public:

  friend class cgTerminationCallback<REAL,ELEMENT_TYPE,ARRAY_TYPE>;

  // Constructor
  //

  cgSolver() : linearSolver<REAL, ELEMENT_TYPE, ARRAY_TYPE>() {
    alpha_ = ELEMENT_TYPE(NAN);
    iterations_ = 10;
    tc_tolerance_ = (REAL)1e-3;
    cb_ = boost::shared_ptr< relativeResidualTCB<REAL, ELEMENT_TYPE, ARRAY_TYPE> >
      ( new relativeResidualTCB<REAL, ELEMENT_TYPE, ARRAY_TYPE>() );
  }
  

  // Destructor
  //

  virtual ~cgSolver() {}


  // Set preconditioner
  //

  virtual void set_preconditioner( boost::shared_ptr< cgPreconditioner<ARRAY_TYPE> > precond ) {
    precond_ = precond;
  }
  

  // Set termination callback
  //

  virtual void set_termination_callback( boost::shared_ptr< cgTerminationCallback<REAL, ELEMENT_TYPE, ARRAY_TYPE> > cb ){
    cb_ = cb;
  }
  

  // Set/get maximally allowed number of iterations
  //

  virtual void set_max_iterations( unsigned int iterations ) { iterations_ = iterations; }
  virtual unsigned int get_max_iterations() { return iterations_; }  


  // Set/get tolerance threshold for termination criterium
  //

  virtual void set_tc_tolerance( REAL tolerance ) { tc_tolerance_ = tolerance; }
  virtual REAL get_tc_tolerance() { return tc_tolerance_; }
  

  // Pre/post solver callbacks
  //
  // 'pre_solve' is invoked _both_ by 'solve' and 'solve_from_rhs'. 
  // - and since 'solve' itself invokes 'solve_from_rhs' the callback is triggered twice in 'solve'
  //
  // 'post_solve' is invoked on the resulting image right before it is returned
  //

  virtual bool pre_solve( ARRAY_TYPE** ) { return true; }
  virtual bool post_solve( boost::shared_ptr<ARRAY_TYPE>& ) { return true; }


  // Pure virtual functions defining core solver functionality
  // Implemented on the host/device respectively in a derived class
  //

  virtual ELEMENT_TYPE solver_dot( ARRAY_TYPE*, ARRAY_TYPE* ) = 0;
  virtual bool solver_clear( ARRAY_TYPE* ) = 0;
  virtual bool solver_scal( ELEMENT_TYPE, ARRAY_TYPE* ) = 0;
  virtual bool solver_axpy( ELEMENT_TYPE, ARRAY_TYPE*, ARRAY_TYPE* ) = 0;
  virtual bool solver_dump( ARRAY_TYPE* ) { return true; }


  //
  // Main solver interface
  //

  virtual boost::shared_ptr<ARRAY_TYPE> solve( ARRAY_TYPE *_d )
  {
    
    // Make copy of the input pointer for the pre_solve callback
    //
    
    ARRAY_TYPE *d = _d;
    if( !d->get_data_ptr() ){
      this->solver_error( "Error: cgSolver::solve : memory allocation failed" );
      return boost::shared_ptr<ARRAY_TYPE>();
    }
    
    // Custom initialization
    //

    if( !pre_solve( &d ))
      this->solver_warning( "Warning: cgSolver::solve : pre_solve callback failed" );
    
    // Compute right hand side
    //
    
    boost::shared_ptr<ARRAY_TYPE> rhs = compute_rhs( d );
    if( !rhs.get() ){
      this->solver_error( "Error: cgSolver::solve : compute_rhs failed" );
      return boost::shared_ptr<ARRAY_TYPE>();
    }
    
    // Solve for the given rhs
    //

    boost::shared_ptr<ARRAY_TYPE> result = solve_from_rhs( rhs.get() );

    if( !result.get() ){
      this->solver_error( "Error: cgSolver::solve : solve_from_rhs failed" );
      return boost::shared_ptr<ARRAY_TYPE>();
    }
    
    return result;
  }


  //
  // Alternative solver interface when given the right hand side
  //

  virtual boost::shared_ptr<ARRAY_TYPE> solve_from_rhs( ARRAY_TYPE *_rhs ) 
  {
    // Make copy of the input pointer for the pre_solve callback
    //
    
    ARRAY_TYPE *rhs = _rhs;
    if( !rhs->get_data_ptr() ){
      this->solver_error( "Error: cgSolver::solve_from_rhs : memory allocation failed" );
      return boost::shared_ptr<ARRAY_TYPE>();
    }
    
    // Custom initialization
    //

    if( !pre_solve( &rhs ))
      this->solver_warning( "Warning: cgSolver::solve_from_rhs : pre_solve callback failed" );
    
    // Initialize
    //

    if( !this->initialize( rhs )){
      this->solver_error( "Error: cgSolver::solve_from_rhs : initialization failed" );
      return boost::shared_ptr<ARRAY_TYPE>();
    }
    
    // Iterate
    //

    if( this->output_mode_ >= solver<ARRAY_TYPE,ARRAY_TYPE>::OUTPUT_VERBOSE ){
      std::cout << "Iterating..." << std::endl;
    }
    
    for( unsigned int it=0; it<iterations_; it++ ){

      REAL tc_metric;
      bool tc_terminate;
      
      if( !this->iterate( it, &tc_metric, &tc_terminate )){
	this->solver_error( "Error: cgSolver::solve_from_rhs : iteration failed" );
	deinitialize();
	return boost::shared_ptr<ARRAY_TYPE>();
      }

      if( !solver_dump( x_.get()) )
	this->solver_warning( "Warning: cgSolver::solve_from_rhs : image dump callback failed" );
      
      if( tc_terminate )
	break;
    }

    // Invoke post solve callback
    //

    if( !post_solve( x_ ) )
      this->solver_warning( "Warning : cgSolver::solve_from_rhs : post_solve callback failed" );
    
    // Clean up and we are done
    //

    boost::shared_ptr<ARRAY_TYPE> tmpx = x_;
    deinitialize();
    return tmpx;
  }


  // Compute right hand side
  //

  virtual boost::shared_ptr<ARRAY_TYPE> compute_rhs( ARRAY_TYPE *d )
  {
    
    if( this->encoding_operator_.get() == 0 ){
      this->solver_error( "Error: cgSolver::compute_rhs : no encoding operator is set" );
      return boost::shared_ptr<ARRAY_TYPE>();
    } 
        
    // Get image space dimensions from the encoding operator
    //

    boost::shared_ptr< std::vector<unsigned int> > image_dims = this->encoding_operator_->get_domain_dimensions();
    if( image_dims->size() == 0 ){
      this->solver_error( "Error: cgSolver::compute_rhs : encoding operator has not set domain dimension" );
      return boost::shared_ptr<ARRAY_TYPE>();
    }

    // Create result array
    //

    boost::shared_ptr<ARRAY_TYPE> result = boost::shared_ptr<ARRAY_TYPE>(new ARRAY_TYPE());
    if( !result->create( image_dims.get() )) {
      this->solver_error( "Error: cgSolver::compute_rhs : memory allocation failed (1)" );
      return boost::shared_ptr<ARRAY_TYPE>();
    }
    
    // Clear result
    //

    if( !solver_clear( result.get() )){
      this->solver_error( "Error: cgSolver::compute_rhs : failed to clear result" );
      return boost::shared_ptr<ARRAY_TYPE>();
    }
    
    // Create temporary array
    //

    ARRAY_TYPE tmp;
    if( !tmp.create( image_dims.get() )) {
      this->solver_error( "Error: cgSolver::compute_rhs : memory allocation failed (2)" );
      return boost::shared_ptr<ARRAY_TYPE>();
    }
    
    // Compute operator adjoint
    //

    if( this->encoding_operator_->mult_MH( d, &tmp ) < 0 ) {
      this->solver_error( "Error: cgSolver::compute_rhs : failed to apply matrix operator" );
      return boost::shared_ptr<ARRAY_TYPE>();
    }
    
    // Apply weight
    //

    if( !solver_axpy(this->encoding_operator_->get_weight(), &tmp, result.get() )) {
      this->solver_error( "Error: cgSolver::compute_rhs : failed to accumulate result" );
      return boost::shared_ptr<ARRAY_TYPE>();
    }
    
    return result;
  }

protected:
  
  //
  // Everything beyond this point is internal to the implementation
  // and not intended to be exposed as a public interface
  //

  // Initialize solver
  //

  virtual bool initialize( ARRAY_TYPE *rhs ) 
  {
    // Input validity test
    //

    if( !rhs || rhs->get_number_of_elements() == 0 ){
      this->solver_error( "Error: cgSolver::initialize : empty or NULL rhs provided" );
      return false;
    }
    
    // Result, x
    //

    x_ = boost::shared_ptr<ARRAY_TYPE>( new ARRAY_TYPE() );
    
    if( !x_->create( rhs->get_dimensions().get() )) {
      this->solver_error( "Error: cgSolver::initialize : Unable to allocate temporary storage (1)" );
      return false;
    }
    
    // Initialize r
    //

    r_ = boost::shared_ptr<ARRAY_TYPE>( new ARRAY_TYPE(*rhs) );
    
    // Initialize x,p
    //

    if( !this->get_x0().get() ){ // no starting image provided      
      
      if( !solver_clear( x_.get()) ){ // set x to zero
	this->solver_error( "Error: cgSolver::initialize : Unable to clear result" );
	return false;
      }
    }
    p_ = boost::shared_ptr<ARRAY_TYPE>( new ARRAY_TYPE(*r_) );

    // Apply preconditioning, twice (should change preconditioners to do this)
    //
      
    if( precond_.get() ) {	
      if( precond_->apply( p_.get(), p_.get() ) < 0 ) {
	this->solver_error( "Error: cgSolver::initialize : failed to apply preconditioner to p (1)" );
	return false;
      }
      if( precond_->apply( p_.get(), p_.get() ) < 0 ) {
	this->solver_error( "Error: cgSolver::initialize : failed to apply preconditioner to p (2)" );
	return false;
      }
    }
    rq0_ = real( solver_dot( r_.get(), p_.get() ));

    if (this->get_x0().get()){
	
      if( !this->get_x0()->dimensions_equal( rhs )){
	this->solver_error( "Error: cgSolver::initialize : RHS and initial guess must have same dimensions" );
	return false;
      }
	
      *x_ = *(this->get_x0());
	
      ARRAY_TYPE mhmX;
      if( !mhmX.create( rhs->get_dimensions().get() )) {
	this->solver_error( "Error: cgSolver::initialize : Unable to allocate temporary storage (2)" );
	return false;
      }

      if( this->output_mode_ >= solver<ARRAY_TYPE,ARRAY_TYPE>::OUTPUT_VERBOSE ) {
	std::cout << "Preparing guess..." << std::endl;
      }

      if( !mult_MH_M( this->get_x0().get(), &mhmX )){
	this->solver_error( "Error: cgSolver::initialize : Error in performing mult_MH_M for initial guess" );
	return false;
      }

      if( !solver_axpy( -ELEMENT_TYPE(1), &mhmX, r_.get() )){
	this->solver_error( "Error: cgSolver::initialize : Error in performing axpy for initial guess" );
	return false;
      }

      p_ = boost::shared_ptr<ARRAY_TYPE>( new ARRAY_TYPE(*r_) );

      // Apply preconditioning, twice (should change preconditioners to do this)
      //

      if( precond_.get() ) {
	if( precond_->apply( p_.get(), p_.get() ) < 0 ) {
	  this->solver_error( "Error: cgSolver::initialize : failed to apply preconditioner to p (3)" );
	  return false;
	}
	if( precond_->apply( p_.get(), p_.get() ) < 0 ) {
	  this->solver_error( "Error: cgSolver::initialize : failed to apply preconditioner to p (4)" );
	  return false;
	}
      }
    }

    rq_ = real( solver_dot( r_.get(), p_.get() ));

    // Invoke termination callback initialization
    //
    
    if( !cb_->initialize(this) ){
      this->solver_error( "Error: cgSolver::initialize : termination callback initialization failed" );
      return false;
    }    

    return true;
  }
  

  // Clean up
  //

  virtual bool deinitialize() 
  {
    p_ = boost::shared_ptr<ARRAY_TYPE>();
    r_ = boost::shared_ptr<ARRAY_TYPE>();
    x_ = boost::shared_ptr<ARRAY_TYPE>();
    return true;
  }


  // Perform full cg iteration
  //

  virtual bool iterate( unsigned int iteration, REAL *tc_metric, bool *tc_terminate )
  {
    ARRAY_TYPE q;
    if( !q.create( x_->get_dimensions().get() )) {
      this->solver_error( "Error: cgSolver::iterate : memory allocation failed for temporary image" );
      return false;
    }

    // Perform one iteration of the solver
    //

    if( !mult_MH_M( p_.get(), &q )){
      this->solver_error( "Error: cgSolver::iterate : mult_MH_M failed" );
      return false;
    }
    
    alpha_ = rq_/solver_dot( p_.get(), &q );

    // Update solution
    //

    if( !solver_axpy( alpha_, p_.get(), x_.get()) ) {
      this->solver_error( "Error: cgSolver::iterate : failed to update solution" );
      return false;
    }

    // Update residual
    //

    if( !solver_axpy( -alpha_, &q, r_.get() )) {
      this->solver_error( "Error: cgSolver::iterate : failed to update residual" );
      return false;
    }

    // Apply preconditioning
    //

    if( precond_.get() ){

      if( precond_->apply( r_.get(), &q ) < 0 ) {
	this->solver_error( "Error: cgSolver::iterate : failed to apply preconditioner to q (1)" );
	return false;
      }
      if( precond_->apply( &q, &q ) < 0 ) {
	this->solver_error( "Error: cgSolver::iterate : failed to apply preconditioner to q (2)" );
	return false;
      }
      
      REAL tmp_rq = real(solver_dot( r_.get(), &q ));
      
      if( !solver_scal( (tmp_rq/rq_)*ELEMENT_TYPE(1), p_.get() )){
	this->solver_error( "Error: cgSolver::iterate : failed to scale p (1)" );
	return false;
      }

      if( !solver_axpy( ELEMENT_TYPE(1), &q, p_.get() )) {
	this->solver_error( "Error: cgSolver::iterate : failed to update p (1)" );
	return false;
      }

      rq_ = tmp_rq;
    } 
    else{

      REAL tmp_rq = real(solver_dot( r_.get(), r_.get()) );
      
      if( !solver_scal( (tmp_rq/rq_)*ELEMENT_TYPE(1), p_.get() )){
	this->solver_error( "cgSolver::iterate : failed to scale p (2)" );
	return false;
      }
      
      if( !solver_axpy( ELEMENT_TYPE(1), r_.get(), p_.get() )) {
	this->solver_error( "cgSolver::iterate : failed to update p (2)" );
	return false;
      }

      rq_ = tmp_rq;      
    }
    
    // Invoke termination callback iteration
    //

    if( !cb_->iterate( iteration, tc_metric, tc_terminate ) ){
      this->solver_error( "Error: cgSolver::iterate : termination callback iteration failed" );
      return false;
    }    

    return true;
  }


  // Perform mult_MH_M of the encoding and regularization matrices
  //

  bool mult_MH_M( ARRAY_TYPE *in, ARRAY_TYPE *out )
  {
    // Basic validity checks
    //

    if( !in || !out ){
      this->solver_error( "Error: cgSolver::mult_MH_M : invalid input pointer(s)" );
      return false;
    }

    if( in->get_number_of_elements() != out->get_number_of_elements() ){
      this->solver_error( "Error: cgSolver::mult_MH_M : array dimensionality mismatch" );
      return false;
    }
    
    // Intermediate storage
    //

    ARRAY_TYPE q;
    if( !q.create( in->get_dimensions().get() )) {
      this->solver_error( "Error: cgSolver::mult_MH_M : memory allocation failed for temporary image" );
      return false;
    }

    // Start by clearing the output
    //

    if( !solver_clear( out )){
      this->solver_error( "Error: cgSolver::mult_MH_M : failed to clear output" );
      return false;
    }

    // Apply encoding operator
    //

    if( this->encoding_operator_->mult_MH_M( in, &q, false ) < 0 ) {
      this->solver_error( "Error: cgSolver::mult_MH_M : failed to apply encoding operator" );
      return false;
    }
    
    if( !solver_axpy( this->encoding_operator_->get_weight(), &q, out )) {
      this->solver_error( "Error: cgSolver::mult_MH_M : failed to add result from  encoding operator" );
      return false;
    }

    // Iterate over regularization operators
    //

    for( unsigned int i=0; i<this->regularization_operators_.size(); i++ ){
      
      if( this->regularization_operators_[i]->mult_MH_M( in, &q, false ) < 0 ) {
	this->solver_error( "Error: cgSolver::mult_MH_M : failed to apply linear operator" );
	return false;
      }
      
      if( !solver_axpy( this->regularization_operators_[i]->get_weight(), &q, out )) {
	this->solver_error( "Error: cgSolver::mult_MH_M : failed to add result from linear operator" );
	return false;
      }
    }

    return true;
  }

protected:

  // Preconditioner
  boost::shared_ptr< cgPreconditioner<ARRAY_TYPE> > precond_;

  // Termination criterium callback
  boost::shared_ptr< cgTerminationCallback<REAL, ELEMENT_TYPE, ARRAY_TYPE> > cb_;

  // Termination criterium threshold
  REAL tc_tolerance_;

  // Maximum number of iterations
  unsigned int iterations_;

  // Internal variables. 
  // Exposed by member functions to allow access for the termination criterium callbacks.

  REAL rq_;
  REAL rq0_;
  ELEMENT_TYPE alpha_;
  boost::shared_ptr<ARRAY_TYPE> x_, p_, r_;
};
