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

  cgSolver() : solver<ARRAY_TYPE,ARRAY_TYPE>() { 
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

  virtual boost::shared_ptr<ARRAY_TYPE> solve( ARRAY_TYPE *_rhs ){
    return solve(_rhs,0);
  }

  virtual std::vector<REAL> get_convergence_list() {
      return ccList;
  }

  virtual bool prepare_for_iterating( ARRAY_TYPE *_rhs, ARRAY_TYPE * x0 ) {
    // Input validity test
    if( !_rhs || _rhs->get_number_of_elements() == 0 ){
      this->solver_error( "cgSolver::solve : empty or NULL rhs provided" );
      return false;
    }

    // Make copy of the input pointer for the pre_solve callback
    ARRAY_TYPE *rhs = _rhs;

    // Custom initialization
    if( !pre_solve(&rhs) ){
      this->solver_error( "cgSolver::solve : error in pre_solve" );
      return false;
    }

    // Result, rho
    x = new ARRAY_TYPE();

    if( !x->create(rhs->get_dimensions().get() )) {
      this->solver_error( "cgSolver::solve : Unable to allocate temp storage (x)" );
      return false;
    }

    if (x0){
      *x = *x0;
    } else{
      // Clear rho
      solver_clear(x);
    }

    // Calculate residual r
    r =  *rhs;

    if( !q.create(rhs->get_dimensions().get() )) {
      this->solver_error( "cgSolver::solve : Unable to allocate temp storage (q)" );
      return false;
    }

    q = r;

    if( precond_.get() ) {
        // Apply preconditioning, twice. Should change preconditioners to do this
        if( precond_->apply(&q,&q) < 0 ) {
            this->solver_error( "cgSolver::solve : failed to apply preconditioner to q" );
            return false;
        }

        if( precond_->apply(&q,&q) < 0 ) {
            this->solver_error( "cgSolver::solve : failed to apply preconditioner to q" );
            return false;
        }
    }

    rq_0 = real(solver_dot(&r, &q));
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

        if (!solver_axpy(-ELEMENT_TYPE(1),&mhmX,&r)) {
            this->solver_error( "cgSolver::solve : Error in performing axpy for initial guess" );
            return false;
        }

        q = r;

        if( precond_.get() ) {
            // Apply preconditioning, twice. Should change preconditioners to do this
            if( precond_->apply(&q,&q) < 0 ) {
                this->solver_error( "cgSolver::solve : failed to apply preconditioner to q" );
                return false;
            }

            if( precond_->apply(&q,&q) < 0 ) {
                this->solver_error( "cgSolver::solve : failed to apply preconditioner to q" );
                return false;
            }
        }

        rq = real(solver_dot(&r, &q));
    }

    if( !p.create( rhs->get_dimensions().get() )) {
        this->solver_error( "cgSolver::solve : Unable to allocate temp storage (p)" );
        return false;
    }

    p = q;

    rq_new = rq;

    return true;
  }

  virtual bool iterate(REAL* convergence_criteria) {
      if (!mult_MH_M(&p,&q)){
          this->solver_error( "cgSolver::solve : error in performing mult_MH_M" );
          return false;
      }
      alpha = rq/solver_dot(&p,&q);

      // Update solution
      if( !solver_axpy(alpha,&p,x) ) {
          this->solver_error( "cgSolver::solve : failed to update solution" );
          return false;
      }

      if( !solver_dump(x) ) {
          this->solver_error( "cgSolver::solve : failed to dump" );
          return false;
      }

      // Update residual
      if( !solver_axpy(-alpha,&q,&r) ) {
          this->solver_error( "cgSolver::solve : failed to update residual" );
          return false;
      }

      if (precond_.get()){
          if( precond_->apply(&r,&q) < 0 ) {
              this->solver_error( "cgSolver::solve : failed to apply preconditioner to q" );
              return false;
          }
          if( precond_->apply(&q,&q) < 0 ) {
              this->solver_error( "cgSolver::solve : failed to apply preconditioner to q" );
              return false;
          }

          rq_new = real(solver_dot(&r, &q));
          solver_scal((rq_new/rq)*ELEMENT_TYPE(1),&p);

          if( !solver_axpy(ELEMENT_TYPE(1),&q,&p) ) {
              this->solver_error( "cgSolver::solve : failed to update solution" );
              return false;
          }
      } else{
          rq_new = real(solver_dot(&r, &r));
          solver_scal((rq_new/rq)*ELEMENT_TYPE(1),&p);
          if( !solver_axpy(ELEMENT_TYPE(1),&r,&p) ) {
              this->solver_error( "cgSolver::solve : failed to update solution" );
              return false;
          }
      }

      rq = rq_new;

      // Calculate relative residual norm
      *convergence_criteria = rq/rq_0;

      return true;
  }

  virtual boost::shared_ptr<ARRAY_TYPE> solve( ARRAY_TYPE *_rhs, ARRAY_TYPE * x0 ) {
      bool status = this->prepare_for_iterating(_rhs, x0);
      if (status == false)
          return boost::shared_ptr<ARRAY_TYPE>();

      if( this->output_mode_ >= solver<ARRAY_TYPE,ARRAY_TYPE>::OUTPUT_VERBOSE ) {
          std::cout << "Iterating..." << std::endl;
      }

      REAL cc_last =get_max<REAL>();
      for (unsigned int it=0; it< iterations_; it++ ) {
          REAL convergence_criteria;
          status = this->iterate(&convergence_criteria);
          if (status == false)
              return boost::shared_ptr<ARRAY_TYPE>(x);
          
          
          if( this->output_mode_ >= solver<ARRAY_TYPE,ARRAY_TYPE>::OUTPUT_WARNINGS ) {
              if( this->output_mode_ >= solver<ARRAY_TYPE,ARRAY_TYPE>::OUTPUT_VERBOSE ) {
                  std::cout << "Iteration " << it+1 << ". rq/rq_0 = " << convergence_criteria << std::endl;
              }
              ccList.push_back(convergence_criteria);
              if( cc_last-convergence_criteria < REAL(0) ) {
                  std::cout << "----- Warning: CG residual increase. Stability problem! -----" << std::endl;
              }
          }

          if( convergence_criteria < limit_ ) {
              break;
          } else {
              cc_last = convergence_criteria;
          }
      }

      if( !post_solve(&x) ){
          this->solver_error( "cgSolver::solve : error in post_solve" );
      }

      return boost::shared_ptr<ARRAY_TYPE>(x);
  }

protected:
  bool mult_MH_M( ARRAY_TYPE *in, ARRAY_TYPE *out ){

    solver_clear(out);
    ARRAY_TYPE q2;

    if( !q2.create( in->get_dimensions().get() )) {
      this->solver_error( "cgSolver::solve : Unable to allocate temp storage (q2)" );
      return false;
    }

    for (unsigned int i = 0; i < operators_->size(); i++) {
      if( (*operators_)[i]->mult_MH_M(in, &q2, false) < 0 ) {
          this->solver_error( "cgSolver::solve : failed to apply matrix operator" );
          return false;
      }
      if( !solver_axpy((*operators_)[i]->get_weight(), &q2, out) ) {
          this->solver_error( "cgSolver::solve : failed to add result from operator" );
          return false;
      }
    }
    return true;
  }

protected:
  boost::shared_ptr< std::vector< boost::shared_ptr< matrixOperator<REAL, ARRAY_TYPE> > > > operators_;
  boost::shared_ptr< cgPreconditioner<ARRAY_TYPE> > precond_;
  unsigned int iterations_;
  REAL limit_;
  std::vector<REAL> ccList;

  REAL rq_new, rq, rq_0;
  ELEMENT_TYPE alpha;
  ARRAY_TYPE p, r, q;
  ARRAY_TYPE *x;
};
