#pragma once


#include "linearOperator.h"
#include "linearOperatorSolver.h"
#include "cgPreconditioner.h"
#include "real_utilities.h"

#include <vector>
#include <iostream>
#include "encodingOperatorContainer.h"

namespace Gadgetron {
template <class ARRAY_TYPE> class lsqrSolver: public linearOperatorSolver<ARRAY_TYPE>
{
protected:
	typedef typename ARRAY_TYPE::element_type ELEMENT_TYPE;
	typedef typename realType<ELEMENT_TYPE>::Type REAL;
public:

	lsqrSolver()  {
		iterations_ = 10;
		tc_tolerance_ = (REAL)1e-3;

	}

	virtual ~lsqrSolver() {}
/*
	virtual int set_preconditioner( boost::shared_ptr< cgPreconditioner<ARRAY_TYPE> > precond ) {
		precond_ = precond;
		return 0;
	}
*/
	virtual void set_tc_tolerance( REAL tolerance ) { tc_tolerance_ = tolerance; }
	virtual REAL get_tc_tolerance() { return tc_tolerance_; }

	virtual void set_max_iterations( unsigned int iterations ) { iterations_ = iterations; }
	virtual unsigned int get_max_iterations() { return iterations_; }





	virtual boost::shared_ptr<ARRAY_TYPE> solve( ARRAY_TYPE *b )
  {

		boost::shared_ptr< std::vector<size_t> > image_dims = this->encoding_operator_->get_domain_dimensions();
		if( image_dims->size() == 0 ){
			throw std::runtime_error( "Error: cgSolver::compute_rhs : encoding operator has not set domain dimension" );
		}




		ARRAY_TYPE * x = new ARRAY_TYPE(image_dims);
		clear(x);


		encodingOperatorContainer<ARRAY_TYPE> enc_op;
		boost::shared_ptr<ARRAY_TYPE> u;

		{
			enc_op.add_operator(this->encoding_operator_);
			for (unsigned int i =0; i < this->regularization_operators_.size(); i++)
				enc_op.add_operator(this->regularization_operators_[i]);
			std::vector<ARRAY_TYPE*> encspace(this->regularization_operators_.size()+1,NULL);
			encspace[0] = b;
			u = enc_op.create_codomain(encspace);
		}



		//Initialise u vector
		REAL beta = 0;

		beta = nrm2(u.get());
		*u *= REAL(1)/beta;

		//Initialise v vector
		REAL alpha = 0;
		ARRAY_TYPE v(*x); //v vector is in image space

		clear(&v);


		enc_op.mult_MH(u.get(),&v);


		alpha = nrm2(&v);

		v *= REAL(1)/alpha;

		//Initialise w vector
		ARRAY_TYPE w(v);

		//Initialise phibar
		REAL phibar = beta;
		REAL phibar0 = phibar;

		//Initialise rhobar
		REAL rhobar = alpha;
		REAL rhobar0 = alpha;

		REAL cg_res = alpha;

		REAL rnorm = beta;

		REAL xnorm = 0;
		REAL anorm = 0;
		REAL arnorm = alpha*beta;


		for (int it = 0; it < iterations_; it ++){
			beta = REAL(0);

			*u *= -alpha;

			enc_op.mult_M(&v,u.get(),true);

			beta =nrm2(u.get());
			*u *= REAL(1)/beta;

			v *= -beta;

			enc_op.mult_MH(u.get(),&v,true);
			alpha = nrm2(&v);

			v *= REAL(1)/alpha;


			//Construct and apply next orthogonal transformation
			REAL rho = std::sqrt(norm(rhobar)+norm(beta));
			REAL c = rhobar/rho;
			REAL s = beta/rho;
			REAL theta = s*alpha;
			rhobar = -c*alpha;
			REAL phi = c*phibar;
			phibar *= s;


			//Update x, w
			axpy(phi/rho,&w,x);  //x = x + phi/rho * w

			w *= -theta/rho;
			w += v;

			//Check for convergence

			//rhobar is a good approximation of the euclidian norm of the residual, so we check for that

			if( this->output_mode_ >= solver<ARRAY_TYPE,ARRAY_TYPE>::OUTPUT_VERBOSE ){
				std::cout << "Iteration " <<it << ". Relative residual: " <<  rhobar/rhobar0 << std::endl;
			}

		}



		return boost::shared_ptr<ARRAY_TYPE>(x);

}



protected:

	//boost::shared_ptr< cgPreconditioner<ARRAY_TYPE> > precond_;
	unsigned int iterations_;
	REAL tc_tolerance_;

};

}
