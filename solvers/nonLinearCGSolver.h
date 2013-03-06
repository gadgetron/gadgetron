/*
 * Nonlinear Conjugate gradient using the Polak Ribiere Polyak (PRP) method.
 * Implementation based on
 * Y. H. DAI and Y. YUAN,SIAM J. OPTIM. Vol. 10, No. 1, pp. 177-182, 1999
 */
#pragma once

#include "linearSolver.h"

#include "real_utilities.h"
#include "complext.h"
#include "sparsifyingOperator.h"
#include <vector>
#include <iostream>


template <class REAL, class ELEMENT_TYPE, class ARRAY_TYPE> class nonLinearCGSolver
: public linearSolver<REAL, ELEMENT_TYPE, ARRAY_TYPE> {
	public:
	nonLinearCGSolver(): linearSolver<REAL, ELEMENT_TYPE, ARRAY_TYPE>() {

	    iterations_ = 10;
	    tc_tolerance_ = (REAL)1e-6;
	    non_negativity_constraint_=false;
	    dump_residual = false;
	    restart_interval=20;
	  }
	void add_csOperator(boost::shared_ptr<sparsifyingOperator<REAL,ARRAY_TYPE> > cs){
		csOperators_.push_back(cs);
	}

	boost::shared_ptr<ARRAY_TYPE> solve(ARRAY_TYPE* in){

		if( this->encoding_operator_.get() == 0 ){
		  this->solver_error( "Error: nonLinearCGSolver::compute_rhs : no encoding operator is set" );
		  return boost::shared_ptr<ARRAY_TYPE>();
		}

		// Get image space dimensions from the encoding operator
		//

		boost::shared_ptr< std::vector<unsigned int> > image_dims = this->encoding_operator_->get_domain_dimensions();
		if( image_dims->size() == 0 ){
		  this->solver_error( "Error: nonLinearCGSolver::compute_rhs : encoding operator has not set domain dimension" );
		  return boost::shared_ptr<ARRAY_TYPE>();
		}

		ARRAY_TYPE * x = new ARRAY_TYPE;
		x->create(image_dims.get());


		ARRAY_TYPE * g = new ARRAY_TYPE;
		g->create(image_dims.get());
		ARRAY_TYPE *  g_old = new ARRAY_TYPE;
		g_old->create(image_dims.get());

		ARRAY_TYPE * d = new ARRAY_TYPE;
		d->create(image_dims.get());
		ARRAY_TYPE *  d_old = new ARRAY_TYPE;
		d_old->create(image_dims.get());

		if (this->x0_.get()){
			*x = *(this->x0_.get());
		} else  {
			solver_clear(x);
		}
		//REAL r0 = solver_dot(rhs.get(),rhs.get());
		REAL r0;
		REAL reg_res,data_res;
		ARRAY_TYPE encoding_space;
		encoding_space.create(in->get_dimensions().get());
		ARRAY_TYPE encoding_space2;
		encoding_space2.create(in->get_dimensions().get());
		ARRAY_TYPE encoding_space3;
		encoding_space3.create(in->get_dimensions().get());
	  if( this->output_mode_ >= solver<ARRAY_TYPE,ARRAY_TYPE>::OUTPUT_VERBOSE ){
		  std::cout << "Iterating..." << std::endl;
		}
		for (int i = 0; i < iterations_; i++){
			if ((i==0) && (!this->x0_.get())){
				solver_clear(d);
				solver_clear(g);
				encoding_space = *in;

				this->encoding_operator_->mult_MH(in,x);

				solver_axpy(REAL(-1),x,g);
				solver_clear(x);

			} else {
				this->encoding_operator_->mult_M(x,&encoding_space);
				solver_axpy(REAL(-1),in,&encoding_space);
				 if( this->output_mode_ >= solver<ARRAY_TYPE,ARRAY_TYPE>::OUTPUT_VERBOSE ){
					 std::cout << "Data residual: " << solver_nrm2(&encoding_space) << std::endl;
				}
				 if (dump_residual){
					 data_res = solver_nrm2(&encoding_space);
				 }

				this->encoding_operator_->mult_MH(&encoding_space,g);
				reg_res = sqrt(mult_MH_M(x,g,true));
			}
			if (i==0){
				r0=solver_dot(g,g);
			}


			for (int n  = 0; n < csOperators_.size(); n++){
				csOperators_[n]->gradient(x,g,true);
			}

			if (non_negativity_constraint_) solver_non_negativity_filter(x,g);
			if (i%restart_interval==0){
				solver_axpy(REAL(-1),g,d);
			} else {
				 //Polak Ribiere Polyak (PRP) value of Beta
				REAL gnrm = solver_dot(g_old,g_old);
				solver_axpy(REAL(-1),g,g_old);
				solver_scal(REAL(-1),g_old);
				REAL beta = solver_dot(g,g_old)/gnrm;
				solver_scal(beta,d);
				solver_axpy(REAL(-1),g,d);
			}
			//Debug DEMON Line
			solver_clear(d);
			solver_axpy(REAL(-1),g,d);


			//Start back-stabbing line search
			REAL nabla = REAL(1);
			REAL c1 = REAL(1e-4);
			REAL c2 = REAL(0.5);
			REAL rho = solver_dot(d,g);
			solver_clear(g_old);
			for (int n  = 0; n < csOperators_.size(); n++){
				csOperators_[n]->apply(x,g_old,true);
			}
			REAL f0 = solver_asum(g_old);
			REAL  reg_norm  = norm_mult_M(x);
			f0 += reg_norm*reg_norm+solver_dot(&encoding_space,&encoding_space);
			REAL f_rhs = f0+c1*nabla*rho;
			ARRAY_TYPE x2 = *x;
			solver_axpy(nabla,d,&x2);
			solver_clear(g_old);
			for (int n  = 0; n < csOperators_.size(); n++){
							csOperators_[n]->apply(&x2,g_old,true);
			}

			this->encoding_operator_->mult_M(d,&encoding_space2,false);
			encoding_space3 = encoding_space;
			solver_axpy(nabla,&encoding_space2,&encoding_space3);
			REAL f_lhs = solver_asum(g_old);
			reg_norm  = norm_mult_M(&x2);
			f_lhs += reg_norm*reg_norm+solver_dot(&encoding_space3,&encoding_space3);
			while (f_lhs > f_rhs){
				std::cout << f0 << " " << f_lhs << " " << f_rhs << " "<< nabla << " "  << rho<< std::endl;
				nabla = c2*nabla;
				x2 = *x;
				solver_axpy(nabla,d,&x2);
				solver_clear(g_old);
				for (int n  = 0; n < csOperators_.size(); n++){
								csOperators_[n]->apply(&x2,g_old,true);
				}

				encoding_space3 = encoding_space;
				solver_axpy(nabla,&encoding_space2,&encoding_space3);
				f_lhs = solver_asum(g_old);
				reg_norm  = norm_mult_M(&x2);
				f_lhs += reg_norm*reg_norm+solver_dot(&encoding_space3,&encoding_space3);
				f_rhs = f0+c1*nabla*rho;
			}
			solver_axpy(nabla,d,x);

			if (non_negativity_constraint_) solver_threshold(0,x);


			ARRAY_TYPE * tmp;
			tmp=g_old;
			g_old=g;
			g=tmp;
			REAL rel_res = (solver_nrm2(d)*nabla/r0) ;
			if (i==0) rel_res=REAL(1);
			if( this->output_mode_ >= solver<ARRAY_TYPE,ARRAY_TYPE>::OUTPUT_VERBOSE ){
			  std::cout << "Iteration " <<i << ". Relative residual change: " <<  rel_res << std::endl;
			}
			iteration_callback(x,i,reg_res,data_res);
			if (rel_res < tc_tolerance_) break;
		}


		delete g,g_old;

		return boost::shared_ptr<ARRAY_TYPE>(x);



	}
	virtual ~nonLinearCGSolver(){};
	  // Set preconditioner
	  //
	  /*virtual void set_preconditioner( boost::shared_ptr< cgPreconditioner<ARRAY_TYPE> > precond ) {
	    precond_ = precond;
	  }*/

	  // Set/get maximally allowed number of iterations
	  //
	  virtual void set_max_iterations( unsigned int iterations ) { iterations_ = iterations; }
	  virtual unsigned int get_max_iterations() { return iterations_; }

	  // Set/get tolerance threshold for termination criterium
	  //
	  virtual void set_tc_tolerance( REAL tolerance ) { tc_tolerance_ = tolerance; }
	  virtual REAL get_tc_tolerance() { return tc_tolerance_; }

	  virtual void set_non_negativity_constraint(bool non_negativity_constraint){
		  non_negativity_constraint_=non_negativity_constraint;
	  }

	  virtual void set_dump_residual(bool dump_res){
		  dump_residual = dump_res;
	  }
	  // Pure virtual functions defining core solver functionality
	  // Implemented on the host/device respectively in a derived class


	  virtual ELEMENT_TYPE solver_dot( ARRAY_TYPE*, ARRAY_TYPE* ) = 0;
	  virtual bool solver_clear( ARRAY_TYPE* in_out,ELEMENT_TYPE v = ELEMENT_TYPE(0)  ) = 0;
	  virtual REAL solver_nrm2(ARRAY_TYPE*)=0;
	  virtual REAL solver_asum(ARRAY_TYPE*)=0;
	  virtual bool solver_axpy( ELEMENT_TYPE, ARRAY_TYPE*, ARRAY_TYPE* ) = 0;
	  virtual bool solver_non_negativity_filter(ARRAY_TYPE*,ARRAY_TYPE*)=0;
	  virtual bool solver_threshold(ELEMENT_TYPE, ARRAY_TYPE*)=0;
	  virtual bool solver_scal( ELEMENT_TYPE, ARRAY_TYPE* ) = 0;
	  virtual void iteration_callback(ARRAY_TYPE*,int i,REAL,REAL){};



protected:



  REAL norm_mult_M(ARRAY_TYPE *in){
	  REAL result=REAL(0);

	  for (int i = 0; i < this->regularization_operators_.size(); i++){
		  std::vector<unsigned int> dims = *this->regularization_operators_[i]->get_codomain_dimensions();
		  if (dims.size() < 1){ // If codomain size not set, assume same dimesnions as *in
			  dims = *in->get_dimensions().get();
		  }
		  ARRAY_TYPE out;
		  out.create(&dims);

		  this->regularization_operators_[i]->mult_M(in,&out);
		  REAL nrm = solver_nrm2(&out);
		  result += nrm*nrm*this->regularization_operators_[i]->get_weight();
	  }
	  return sqrt(result);


  }

  // Perform mult_MH_M of the encoding and regularization matrices
  //
  REAL mult_MH_M( ARRAY_TYPE *in, ARRAY_TYPE *out,bool accumulate=false )
  {
    // Basic validity checks
    //

    if( !in || !out ){
      this->solver_error( "Error: nonLinearCGSolver::mult_MH_M : invalid input pointer(s)" );
      return false;
    }

    if( in->get_number_of_elements() != out->get_number_of_elements() ){
      this->solver_error( "Error: nonLinearCGSolver::mult_MH_M : array dimensionality mismatch" );
      return false;
    }

    // Intermediate storage
    //

    ARRAY_TYPE q;
    if( !q.create( in->get_dimensions().get() )) {
      this->solver_error( "Error: nonLinearCGSolver::mult_MH_M : memory allocation failed for temporary image" );
      return false;
    }

    // Start by clearing the output
    //
    if (!accumulate){
    if( !solver_clear( out )){
      this->solver_error( "Error: nonLinearCGSolver::mult_MH_M : failed to clear output" );
      return false;
    }}

    // Apply encoding operator
    //
/*
    if( this->encoding_operator_->mult_MH_M( in, &q, false ) < 0 ) {
      this->solver_error( "Error: nonLinearCGSolver::mult_MH_M : failed to apply encoding operator" );
      return false;
    }

    if( !solver_axpy( this->encoding_operator_->get_weight(), &q, out )) {
      this->solver_error( "Error: nonLinearCGSolver::mult_MH_M : failed to add result from  encoding operator" );
      return false;
    }
*/
    // Iterate over regularization operators
    //
    REAL res = 0;
    for( unsigned int i=0; i<this->regularization_operators_.size(); i++ ){
    	if (dump_residual){
    		this->regularization_operators_[i]->mult_M( in, &q, false );
    		res += pow(solver_nrm2(&q),2);
    		this->regularization_operators_[i]->mult_MH( &q, out, true );
    	}
      if( this->regularization_operators_[i]->mult_MH_M( in, &q, false ) < 0 ) {
	this->solver_error( "Error: nonLinearCGSolver::mult_MH_M : failed to apply linear operator" );
	return false;
      }

      if( !solver_axpy( this->regularization_operators_[i]->get_weight(), &q, out )) {
	this->solver_error( "Error: nonLinearCGSolver::mult_MH_M : failed to add result from linear operator" );
	return false;
      }
    }

    return res;
  }

  // Perform mult_MH_M of the encoding and regularization matrices
  //
  bool mult_M( ARRAY_TYPE *in, ARRAY_TYPE *out,bool accumulate=false )
  {
    // Basic validity checks
    //

    if( !in || !out ){
      this->solver_error( "Error: nonLinearCGSolver::mult_M : invalid input pointer(s)" );
      return false;
    }

    if( in->get_number_of_elements() != out->get_number_of_elements() ){
      this->solver_error( "Error: nonLinearCGSolver::mult_M : array dimensionality mismatch" );
      return false;
    }

    // Intermediate storage
    //

    ARRAY_TYPE q;
    if( !q.create( in->get_dimensions().get() )) {
      this->solver_error( "Error: nonLinearCGSolver::mult_MH_M : memory allocation failed for temporary image" );
      return false;
    }

    // Start by clearing the output
    //
    if (!accumulate){
    if( !solver_clear( out )){
      this->solver_error( "Error: nonLinearCGSolver::mult_MH_M : failed to clear output" );
      return false;
    }}

    // Apply encoding operator
    //
/*
    if( this->encoding_operator_->mult_MH_M( in, &q, false ) < 0 ) {
      this->solver_error( "Error: nonLinearCGSolver::mult_MH_M : failed to apply encoding operator" );
      return false;
    }

    if( !solver_axpy( this->encoding_operator_->get_weight(), &q, out )) {
      this->solver_error( "Error: nonLinearCGSolver::mult_MH_M : failed to add result from  encoding operator" );
      return false;
    }
*/
    // Iterate over regularization operators
    //

    for( unsigned int i=0; i<this->regularization_operators_.size(); i++ ){

      if( this->regularization_operators_[i]->mult_M( in, &q, false ) < 0 ) {
	this->solver_error( "Error: nonLinearCGSolver::mult_M : failed to apply linear operator" );
	return false;
      }

      if( !solver_axpy( this->regularization_operators_[i]->get_weight(), &q, out )) {
	this->solver_error( "Error: nonLinearCGSolver::mult_M : failed to add result from linear operator" );
	return false;
      }
    }

  }

protected:

  // Preconditioner
  //boost::shared_ptr< cgPreconditioner<ARRAY_TYPE> > precond_;
  // Maximum number of iterations
  unsigned int iterations_;
  bool non_negativity_constraint_;
  REAL tc_tolerance_;
  bool dump_residual;
  int restart_interval;
  // Internal variables.
  // Exposed by member functions to allow access for the termination criterium callbacks.

  std::vector<boost::shared_ptr<sparsifyingOperator<REAL,ARRAY_TYPE> > > csOperators_;
};
