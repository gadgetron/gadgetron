#pragma once

#include "solver.h"
#include "matrixOperator.h"
#include "cgPreconditioner.h"
#include "real_utilities.h"
#include "solvers_export.h"

#include <vector>
#include <iostream>

template <class REAL, class ELEMENT_TYPE, class ARRAY_TYPE> class cgSolver : public solver<ARRAY_TYPE, ARRAY_TYPE>
{
public:

  cgSolver( int output_mode = solver<ARRAY_TYPE, ARRAY_TYPE>::OUTPUT_SILENT ) : solver<ARRAY_TYPE, ARRAY_TYPE>( output_mode ) { 
    iterations_ = 10;
    limit_ = (REAL)1e-3;
    operators_ = boost::shared_ptr< std::vector< boost::shared_ptr< matrixOperator<REAL, ARRAY_TYPE> > > >
      (new std::vector< boost::shared_ptr< matrixOperator<REAL, ARRAY_TYPE> > >);
  }

  virtual ~cgSolver() {}

  virtual int add_matrix_operator( boost::shared_ptr< matrixOperator<REAL, ARRAY_TYPE> > op ) {
    operators_->push_back(op);
    return 0;
  }

  virtual int set_preconditioner( boost::shared_ptr< cgPreconditioner<ARRAY_TYPE> > precond ) {
    precond_ = precond;
    return 0;
  }

  virtual void set_limit( REAL limit ) {
    limit_ = limit;
  }

  virtual void set_iterations( unsigned int iterations ) {
    iterations_ = iterations;
  }

  virtual bool pre_solve( ARRAY_TYPE **rhs ) { return true; }
  virtual bool post_solve( ARRAY_TYPE **rho ) { return true; }

  virtual ELEMENT_TYPE solver_dot( ARRAY_TYPE*, ARRAY_TYPE* ) = 0;
  virtual bool solver_clear( ARRAY_TYPE* ) = 0;
  virtual bool solver_scal( ELEMENT_TYPE, ARRAY_TYPE* ) = 0;
  virtual bool solver_axpy( ELEMENT_TYPE, ARRAY_TYPE*, ARRAY_TYPE* ) = 0;
  virtual bool solver_dump( ARRAY_TYPE *rho ) { return true; }

  virtual boost::shared_ptr<ARRAY_TYPE> solve( ARRAY_TYPE *_rhs )
  {
    // Input validity test
    if( !_rhs || _rhs->get_number_of_elements() == 0 ){
      this->solver_error( "cgSolver::solve : empty or NULL rhs provided" );
      return boost::shared_ptr<ARRAY_TYPE>();
    }

    // Make copy of the input pointer for the pre_solve callback
    ARRAY_TYPE *rhs = _rhs;

    // Custom initialization
    if( !pre_solve(&rhs) ){
      this->solver_error( "cgSolver::solve : error in pre_solve" );
      return boost::shared_ptr<ARRAY_TYPE>();
    }

    // Result, rho
    ARRAY_TYPE *rho = new ARRAY_TYPE();

    if( !rho->create(rhs->get_dimensions().get() )) {
      this->solver_error( "cgSolver::solve : Unable to allocate temp storage (rho)" );
      return boost::shared_ptr<ARRAY_TYPE>(rho);
    }

    // Clear rho
    solver_clear(rho);

    // Calculate residual r
    ARRAY_TYPE r;
    if( precond_.get() ) {
      if( !r.create( rhs->get_dimensions().get() )) {
	this->solver_error( "cgSolver::solve : Unable to allocate storage (r)" );
	return boost::shared_ptr<ARRAY_TYPE>(rho);
      }
      if( precond_->apply( rhs, &r ) < 0 ) {
	this->solver_error( "cgSolver::solve : Unable to apply preconditioning to rhs" );
	return boost::shared_ptr<ARRAY_TYPE>(rho);
      }
    } else {
      r =  *rhs;
    }

    REAL rr_0    = real<REAL>(solver_dot(&r, &r));
    REAL rr_1    = rr_0;
    REAL rr      = get_zero<REAL>();
    REAL rr_last = get_max<REAL>();

    ARRAY_TYPE p;
    if( !p.create( rhs->get_dimensions().get() )) {
      this->solver_error( "cgSolver::solve : Unable to allocate temp storage (p)" );
      return boost::shared_ptr<ARRAY_TYPE>(rho);
    }

    ARRAY_TYPE p_precond;
    if( precond_.get() ) { // We only need this additional storage if we are using a preconditioner
      if( !p_precond.create( rhs->get_dimensions().get() )) {
	this->solver_error( "cgSolver::solve : Unable to allocate temp storage (p_precond)" );
	return boost::shared_ptr<ARRAY_TYPE>(rho);
      }
    }

    ARRAY_TYPE q;
    if( !q.create( rhs->get_dimensions().get() )) {
      this->solver_error( "cgSolver::solve : Unable to allocate temp storage (q)" );
      return boost::shared_ptr<ARRAY_TYPE>(rho);
    }

    ARRAY_TYPE q2;
    if( !q2.create( rhs->get_dimensions().get() )) {
      this->solver_error( "cgSolver::solve : Unable to allocate temp storage (q2)" );
      return boost::shared_ptr<ARRAY_TYPE>(rho);
    }

    REAL rel_res;

    if( this->output_mode_ >= solver<ARRAY_TYPE,ARRAY_TYPE>::OUTPUT_VERBOSE ) {
      std::cout << "Iterating..." << std::endl;
    }

    for( unsigned int it = 0; it < iterations_; it++ ) {

      rr_1 = rr;
      rr = real<REAL,ELEMENT_TYPE>( solver_dot(&r, &r) );

      // Update p
      if( it == 0 ){
	p = r;
      } else {        
	ELEMENT_TYPE beta = mul<REAL>(rr/rr_1, get_one<ELEMENT_TYPE>());
	if( !solver_scal(beta,&p) ) {
	  this->solver_error( "cgSolver::solve : failed to scale p" );
	  return boost::shared_ptr<ARRAY_TYPE>(rho);
	}
	if( !solver_axpy(get_one<ELEMENT_TYPE>(),&r,&p) ) {
	  this->solver_error( "cgSolver::solve : failed to add r to scaled p" );
	  return boost::shared_ptr<ARRAY_TYPE>(rho);
	}
      }

      // Now we need to multiply with the system matrix
      solver_clear(&q);

      // Take care of preconditioning
      ARRAY_TYPE* cur_p = &p;
      if( precond_.get() ) {
	if( precond_->apply(&p,&p_precond) < 0 ) {
	  this->solver_error( "cgSolver::solve : failed to apply preconditioner to p" );
	  return boost::shared_ptr<ARRAY_TYPE>(rho);
	}
	cur_p = &p_precond;
      }

      for (unsigned int i = 0; i < operators_->size(); i++) {

	if( (*operators_)[i]->mult_MH_M(cur_p, &q2, false) < 0 ) {
	  this->solver_error( "cgSolver::solve : failed to apply matrix operator" );
	  return boost::shared_ptr<ARRAY_TYPE>(rho);
	}

	if( !solver_axpy(mul<REAL>((*operators_)[i]->get_weight(), get_one<ELEMENT_TYPE>()), &q2, &q) ) {
	  this->solver_error( "cgSolver::solve : failed to add result from operator" );
	  return boost::shared_ptr<ARRAY_TYPE>(rho);
	}
      }

      if( precond_.get() ) {
	if( precond_->apply(&q,&q) < 0 ) {
	  this->solver_error( "cgSolver::solve : failed to apply preconditioner to q" );
	  return boost::shared_ptr<ARRAY_TYPE>(rho);
	}
      }

      ELEMENT_TYPE alpha = mul<REAL>(rr, reciprocal<ELEMENT_TYPE>(solver_dot(&p,&q)));

      // Update solution
      if( !solver_axpy(alpha,&p,rho) ) {
	this->solver_error( "cgSolver::solve : failed to update solution" );
	return boost::shared_ptr<ARRAY_TYPE>(rho);
      }

    if( !solver_dump(rho) ) {
        this->solver_error( "cgSolver::solve : failed to dump" );
        return boost::shared_ptr<ARRAY_TYPE>(rho);
    }

      // Update residual
      if( !solver_axpy(mul<REAL>(-get_one<REAL>(),alpha),&q,&r) ) {
	this->solver_error( "cgSolver::solve : failed to update residual" );
	return boost::shared_ptr<ARRAY_TYPE>(rho);
      }

      // Calculate relative residual norm
      rel_res = rr/rr_0;

      if( this->output_mode_ >= solver<ARRAY_TYPE,ARRAY_TYPE>::OUTPUT_WARNINGS ) {
	if( this->output_mode_ >= solver<ARRAY_TYPE,ARRAY_TYPE>::OUTPUT_VERBOSE ) {
	  std::cout << "Iteration " << it+1 << ". rr/rr_0 = " << rel_res << std::endl;
	}
	if( rr_last-rel_res < get_zero<REAL>() ) {
	  std::cout << "----- Warning: CG residual increase. Stability problem! -----" << std::endl;
	}
      }

      if( rel_res < limit_ ) {
	break;
      } else {
	rr_last = rel_res;
      }
    }

    if( precond_.get() ) {
      if( precond_->apply(rho,rho) < 0 ) {
	this->solver_error( "cgSolver::solve : failed to apply preconditioner to rho" );
	return boost::shared_ptr<ARRAY_TYPE>(rho);
      }
    }

    if( !post_solve(&rho) ){
      this->solver_error( "cgSolver::solve : error in post_solve" );
      return boost::shared_ptr<ARRAY_TYPE>(rho);
    }

    return boost::shared_ptr<ARRAY_TYPE>(rho);
  }

protected:
  boost::shared_ptr< std::vector< boost::shared_ptr< matrixOperator<REAL, ARRAY_TYPE> > > > operators_;
  boost::shared_ptr< cgPreconditioner<ARRAY_TYPE> > precond_;
  unsigned int iterations_;
  REAL limit_;
};
