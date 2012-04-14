#pragma once

#include "solver.h"
#include "matrixOperator.h"
#include "cgPreconditioner.h"
#include "real_utilities.h"

#include <vector>
#include <iostream>

template <class REAL, class ELEMENT_TYPE, class ARRAY_TYPE> class cgSolver 
  : public solver<ARRAY_TYPE, ARRAY_TYPE>
{
public:

  // Constructor
  cgSolver() : solver<ARRAY_TYPE,ARRAY_TYPE>() { 
    iterations_ = 10;
    limit_ = (REAL)1e-3;
    tc_history_enabled_ = false;
  }

  // Destructor
  virtual ~cgSolver() {}

  // Add matrix operator to the solver
  // ---------------------------------
  // The latter two arguments are used during subsequent calls to 'solve' and 'solve_from_rhs'
  //
  // When using the 'solve_from_rhs' interface, 'contributes_to_rhs' and 'rhs_prior' are ignored.
  // When using the 'solve' interface
  // - 'contributions_to_rhs' indicates if this operator contributes to the right hand side (rhs):
  // - if true, the adjoint matrix operator (op.mult_MH) is computed on 'rhs_prior' during the rhs computation
  // - if true and 'rhs_prior' is 0x0, (op.mult_MH) is computed on the input data to the 'solve' method 

  inline int add_matrix_operator( boost::shared_ptr< matrixOperator<REAL, ARRAY_TYPE> > op,
				  bool contributes_to_rhs = false, 
				  boost::shared_ptr<ARRAY_TYPE> rhs_prior = boost::shared_ptr<ARRAY_TYPE>() ) 
  {
    operators_.push_back(op);
    indicators_.push_back(contributes_to_rhs);
    rhs_priors_.push_back(rhs_prior);
    return 0;
  }

  // Set preconditioner (optional)
  inline int set_preconditioner( boost::shared_ptr< cgPreconditioner<ARRAY_TYPE> > precond ) {
    precond_ = precond;
    return 0;
  }

  // Set termination threshold for convergence criterium
  inline void set_limit( REAL limit ) {
    limit_ = limit;
  }
  
  // Set maximally allowed number of iterations
  inline void set_max_iterations( unsigned int iterations ) {
    iterations_ = iterations;
  }

  // Toggle (on/off) record keeping of termination criterium (tc) evaluations
  void enable_tc_history( bool b ) {
    tc_history_enabled_ = b;
  }

  // Return history of the termination criterium evaluations
  inline std::vector<REAL> get_tc_history() {
    return tc_values_;
  }

  // Pre/post solver callbacks
  //
  // 'pre_solve' is invoked _both_ by 'solve' and 'solve_from_rhs'. 
  // - and since 'solve' itself invokes 'solve_from_rhs' the callback is triggered twice in 'solve'
  //
  // 'post_solve' is invoked right before the image is returned

  virtual bool pre_solve( ARRAY_TYPE** ) { return true; }
  virtual bool post_solve( ARRAY_TYPE** ) { return true; }

  // Pure virtual functions defining core solver functionality
  // Implemented on the host/device respectively in a derived class

  virtual ELEMENT_TYPE solver_dot( ARRAY_TYPE*, ARRAY_TYPE* ) = 0;
  virtual bool solver_clear( ARRAY_TYPE* ) = 0;
  virtual bool solver_scal( ELEMENT_TYPE, ARRAY_TYPE* ) = 0;
  virtual bool solver_axpy( ELEMENT_TYPE, ARRAY_TYPE*, ARRAY_TYPE* ) = 0;
  virtual bool solver_dump( ARRAY_TYPE *x ) { return true; }

  // Inherited solver interface
  virtual boost::shared_ptr<ARRAY_TYPE> solve( ARRAY_TYPE *_d )
  {
    // Make copy of the input pointer for the pre_solve callback
    ARRAY_TYPE *d = _d;
    
    // Custom initialization
    if( !pre_solve(&d) ){
      this->solver_error( "Warning: cgSolver::solve : pre_solve callback failed" );
    }
    
    // Compute right hand side
    boost::shared_ptr<ARRAY_TYPE> rhs = compute_rhs(d);
    if( !rhs.get() ){
      this->solver_error( "Error: cgSolver::solve : compute_rhs failed" );
      return boost::shared_ptr<ARRAY_TYPE>();
    }
    
    // Solve for the given rhs
    boost::shared_ptr<ARRAY_TYPE> result = solve_from_rhs( rhs.get() );

    if( !result.get() ){
      this->solver_error( "Error: cgSolver::solve : solve_from_rhs failed" );
      return boost::shared_ptr<ARRAY_TYPE>();
    }
    
    return result;
  }
  
  // Solver interface given the right hand side
  virtual boost::shared_ptr<ARRAY_TYPE> solve_from_rhs( ARRAY_TYPE *_rhs ) 
  {
    // Make copy of the input pointer for the pre_solve callback
    ARRAY_TYPE *rhs = _rhs;
    
    // Custom initialization
    if( !pre_solve(&rhs) ){
      this->solver_error( "Warning: cgSolver::solve_from_rhs : pre_solve callback failed" );
    }

    // Initialize
    if( !this->initialize(rhs, this->x0_.get()) ){
      this->solver_error( "Error: cgSolver::solve_from_rhs : initialization failed" );
      return boost::shared_ptr<ARRAY_TYPE>();
    }

    // Iterate
    //

    if( this->output_mode_ >= solver<ARRAY_TYPE,ARRAY_TYPE>::OUTPUT_VERBOSE ){
      std::cout << "Iterating..." << std::endl;
    }
    
    REAL cc_last = get_max<REAL>();
    for (unsigned int it=0; it < iterations_; it++ ){
      REAL convergence_criteria;

      if( !this->iterate(&convergence_criteria) ){
	boost::shared_ptr<ARRAY_TYPE> tmpx = x;
	deinitialize();
	return tmpx;
      }
      tc_values_.push_back(convergence_criteria);

      if( this->output_mode_ >= solver<ARRAY_TYPE,ARRAY_TYPE>::OUTPUT_WARNINGS ) {
	if( this->output_mode_ >= solver<ARRAY_TYPE,ARRAY_TYPE>::OUTPUT_VERBOSE ) {
	  std::cout << "Iteration " << it+1 << ". rq/rq_0 = " << convergence_criteria << std::endl;
	}
	if( cc_last-convergence_criteria < REAL(0) ) {
	  std::cout << "----- Warning: CG residual increase. Stability problem! -----" << std::endl;
	}
      }

      if( !solver_dump(x.get()) ) {
	this->solver_error( "Warning: cgSolver::solve_from_rhs : image dump callback failed" );
	boost::shared_ptr<ARRAY_TYPE> tmpx = x;
	deinitialize();
	return tmpx;
      }

      if( convergence_criteria < limit_ ) {
	break;
      } else {
	cc_last = convergence_criteria;
      }
    }

    ARRAY_TYPE* xp = x.get();
    if( !post_solve( &xp ) ){
      this->solver_error( "Warning : cgSolver::solve_from_rhs : post_solve callback failed" );
    }

    boost::shared_ptr<ARRAY_TYPE> tmpx = x;
    deinitialize();
    return tmpx;
  }

protected:

  virtual boost::shared_ptr<ARRAY_TYPE> compute_rhs( ARRAY_TYPE *d )
  {

    // Make sure we have at least one operator on which to compute the adjoint
    // We also fetch the image dimensions from a suitable operator

    int first_idx = -1;
    for( int i=0; i<operators_.size(); i++){
      if( indicators_[i] ){ 
	if( operators_[i]->get_domain_dimensions().size() == 0 ){
	  first_idx--;
	  continue;	  
	}
	first_idx = i;
	break;
      }
    }
    
    if( first_idx == -1 ){
      this->solver_error( "Error: cgSolver::compute_rhs : no matrix operators are configured to contribute to rhs" );
      return boost::shared_ptr<ARRAY_TYPE>();
    }
    
    if( first_idx < -1 ){
      this->solver_error( "Error: cgSolver::compute_rhs : no matrix operators contributing to rhs has defined domain dimensions" );
      return boost::shared_ptr<ARRAY_TYPE>();
    }
    
    // We can take the image space dimensions from the first operator
    std::vector<unsigned int> image_dims = operators_[first_idx]->get_domain_dimensions();
    if( image_dims.size() == 0 ){
      this->solver_error( "Error: cgSolver::compute_rhs : " );
      return boost::shared_ptr<ARRAY_TYPE>();
    }

    // Create result array and clear to zero
    boost::shared_ptr<ARRAY_TYPE> result = boost::shared_ptr<ARRAY_TYPE>(new ARRAY_TYPE());
    if( !result->create( &image_dims )) {
      this->solver_error( "Error: cgSolver::compute_rhs : memory allocation failed (1)" );
      return boost::shared_ptr<ARRAY_TYPE>();
    }
    solver_clear( result.get() );
    
    // Create temporary array
    boost::shared_ptr<ARRAY_TYPE> tmp = boost::shared_ptr<ARRAY_TYPE>(new ARRAY_TYPE());
    if( !tmp->create( &image_dims )) {
      this->solver_error( "Error: cgSolver::compute_rhs : memory allocation failed (2)" );
      return boost::shared_ptr<ARRAY_TYPE>();
    }
    
    // Iterate over operators to compute rhs
    //

    for( unsigned int i=0; i<operators_.size(); i++){
      if( indicators_[i] ) {

	// Compute operator adjoint
	if( operators_[i]->mult_MH( (rhs_priors_[i].get()) ? rhs_priors_[i].get() : d, tmp.get() ) < 0 ) {
	  this->solver_error( "Error: cgSolver::compute_rhs : failed to apply matrix operator" );
	  return boost::shared_ptr<ARRAY_TYPE>();
	}

	// Accumulate
	if( !solver_axpy(operators_[i]->get_weight(), tmp.get(), result.get() )) {
	  this->solver_error( "Error: cgSolver::compute_rhs : failed to accumulate result" );
	  return boost::shared_ptr<ARRAY_TYPE>();
	}
      }
    }
    return result;
  }

  virtual bool initialize( ARRAY_TYPE *rhs, ARRAY_TYPE *x0 ) 
  {
    // Input validity test
    if( !rhs || rhs->get_number_of_elements() == 0 ){
      this->solver_error( "cgSolver::solve : empty or NULL rhs provided" );
      return false;
    }
  
    q2 = new ARRAY_TYPE();
    if( !q2->create( rhs->get_dimensions().get() )) {
      this->solver_error( "cgSolver::initialize : Unable to allocate temp storage (q2)" );
      return false;
    }

    // Result, rho
    x = boost::shared_ptr<ARRAY_TYPE>( new ARRAY_TYPE() );

    if( !x->create(rhs->get_dimensions().get() )) {
      this->solver_error( "cgSolver::solve : Unable to allocate temp storage (x)" );
      return false;
    }

    if (x0){
      *(x.get()) = *x0;
    } else{
      // Clear rho
      solver_clear(x.get());
    }

    // Calculate residual r
    r = new ARRAY_TYPE(*rhs);
    //r =  *rhs;

    q = new ARRAY_TYPE();
    if( !q->create(rhs->get_dimensions().get() )) {
      this->solver_error( "cgSolver::solve : Unable to allocate temp storage (q)" );
      return false;
    }

    *q = *r;

    if( precond_.get() ) {
      // Apply preconditioning, twice. Should change preconditioners to do this
      if( precond_->apply(q,q) < 0 ) {
	this->solver_error( "cgSolver::solve : failed to apply preconditioner to q" );
	return false;
      }

      if( precond_->apply(q,q) < 0 ) {
	this->solver_error( "cgSolver::solve : failed to apply preconditioner to q" );
	return false;
      }
    }

    rq_0 = real(solver_dot(r, q));
    rq = rq_0;

    if (x0) {
      if (!x0->dimensions_equal(rhs)){
	this->solver_error( "cgSolver::solve : RHS and initial guess must have same dimensions" );
	return false;
      }

      ARRAY_TYPE mhmX;
      if( !mhmX.create(rhs->get_dimensions().get() )) {
	this->solver_error( "cgSolver::solve : Unable to allocate temp storage (mhmX)" );
	return false;
      }

      if( this->output_mode_ >= solver<ARRAY_TYPE,ARRAY_TYPE>::OUTPUT_VERBOSE ) {
	std::cout << "Preparing guess..." << std::endl;
      }

      if(!mult_MH_M(x0,&mhmX)){
	this->solver_error( "cgSolver::solve : Error in performing mult_MH_M for initial guess" );
	return false;
      }

      if (!solver_axpy(-ELEMENT_TYPE(1),&mhmX,r)) {
	this->solver_error( "cgSolver::solve : Error in performing axpy for initial guess" );
	return false;
      }

      *q = *r;

      if( precond_.get() ) {
	// Apply preconditioning, twice. Should change preconditioners to do this
	if( precond_->apply(q,q) < 0 ) {
	  this->solver_error( "cgSolver::initialize : failed to apply preconditioner to q" );
	  return false;
	}

	if( precond_->apply(q,q) < 0 ) {
	  this->solver_error( "cgSolver::initialize : failed to apply preconditioner to q" );
	  return false;
	}
      }

      rq = real(solver_dot(r, q));
    }

    p = new ARRAY_TYPE();
    if( !p->create( rhs->get_dimensions().get() )) {
      this->solver_error( "cgSolver::initialize : Unable to allocate temp storage (p)" );
      return false;
    }

    *p = *q;
    rq_new = rq;
    return true;
  }

  virtual bool deinitialize() {
    delete p;
    delete r;
    delete q;
    delete q2;
    x = boost::shared_ptr<ARRAY_TYPE>();
    return true;
  }

  virtual bool iterate(REAL* convergence_criteria) {
    if (!mult_MH_M(p,q)){
      this->solver_error( "cgSolver::solve : error in performing mult_MH_M" );
      return false;
    }
    alpha = rq/solver_dot(p,q);

    // Update solution
    if( !solver_axpy(alpha,p,x.get()) ) {
      this->solver_error( "cgSolver::solve : failed to update solution" );
      return false;
    }

    // Update residual
    if( !solver_axpy(-alpha,q,r) ) {
      this->solver_error( "cgSolver::solve : failed to update residual" );
      return false;
    }

    if (precond_.get()){
      if( precond_->apply(r,q) < 0 ) {
	this->solver_error( "cgSolver::solve : failed to apply preconditioner to q" );
	return false;
      }
      if( precond_->apply(q,q) < 0 ) {
	this->solver_error( "cgSolver::solve : failed to apply preconditioner to q" );
	return false;
      }

      rq_new = real(solver_dot(r, q));
      solver_scal((rq_new/rq)*ELEMENT_TYPE(1),p);

      if( !solver_axpy(ELEMENT_TYPE(1),q,p) ) {
	this->solver_error( "cgSolver::solve : failed to update solution" );
	return false;
      }
    } else{
      rq_new = real(solver_dot(r, r));
      solver_scal((rq_new/rq)*ELEMENT_TYPE(1),p);
      if( !solver_axpy(ELEMENT_TYPE(1),r,p) ) {
	this->solver_error( "cgSolver::solve : failed to update solution" );
	return false;
      }
    }

    rq = rq_new;

    // Calculate relative residual norm
    *convergence_criteria = rq/rq_0;

    return true;
  }

  bool mult_MH_M( ARRAY_TYPE *in, ARRAY_TYPE *out ){

    solver_clear(out);


    for (unsigned int i = 0; i < operators_.size(); i++) {
      if( operators_[i]->mult_MH_M(in, q2, false) < 0 ) {
	this->solver_error( "cgSolver::solve : failed to apply matrix operator" );
	return false;
      }
      if( !solver_axpy(operators_[i]->get_weight(), q2, out) ) {
	this->solver_error( "cgSolver::solve : failed to add result from operator" );
	return false;
      }
    }
    return true;
  }

protected:

  // Vector of matrix operators
  std::vector< boost::shared_ptr< matrixOperator<REAL, ARRAY_TYPE> > > operators_;
  
  // Vector of boolean rhs indicators for the operators
  std::vector<bool> indicators_;
  
  // Vector of rhs priors for the operators
  std::vector< boost::shared_ptr<ARRAY_TYPE> > rhs_priors_;
  
  // Preconditioner
  boost::shared_ptr< cgPreconditioner<ARRAY_TYPE> > precond_;

  // Maximum number of iterations
  unsigned int iterations_;

  // Termination criterium threshold
  REAL limit_;

  // Toogle record keeping of termination criterium
  bool tc_history_enabled_;

  // History of termination criterium values
  std::vector<REAL> tc_values_;

  REAL rq_new, rq, rq_0;
  ELEMENT_TYPE alpha;
  ARRAY_TYPE *p, *r, *q, *q2;
  boost::shared_ptr<ARRAY_TYPE> x;

};
