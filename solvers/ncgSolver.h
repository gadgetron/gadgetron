#pragma once

#include "linearOperatorSolver.h"
#include "real_utilities.h"
#include "complext.h"
#include "cgPreconditioner.h"
#include <vector>
#include <iostream>
#include <numeric>

namespace Gadgetron{
/** Nonlinear conjugate gradient solver.
 * Adapted from Y.H. Dai & Y. Yuan 2001 "An Efficient Hybrid Conjugate Gradient Method for Unconstrained Optimization"
 * Annals of Operations Research, March 2001, Volume 103, Issue 1-4, pp 33-47
 *
 */

template <class ARRAY_TYPE> class ncgSolver : public gpSolver<ARRAY_TYPE>
{
protected:
	typedef typename ARRAY_TYPE::element_type ELEMENT_TYPE;
	typedef typename realType<ELEMENT_TYPE>::Type REAL;
	typedef ARRAY_TYPE ARRAY_CLASS;
	typedef gpSolver<ARRAY_TYPE> GP;
	typedef typename gpSolver<ARRAY_TYPE>::l1GPRegularizationOperator l1GPRegularizationOperator;

public:

	ncgSolver(): gpSolver<ARRAY_TYPE>() {
		iterations_ = 10;
		tc_tolerance_ = (REAL)1e-8;
		non_negativity_constraint_=false;
		dump_residual = false;
		threshold= REAL(1e-8);
		barrier_threshold=1e4;
	}

	virtual ~ncgSolver(){}

	virtual boost::shared_ptr<ARRAY_TYPE> solve(ARRAY_TYPE* in)
									{
		if( this->encoding_operator_.get() == 0 ){
			throw std::runtime_error("Error: ncgSolver::compute_rhs : no encoding operator is set" );
		}

		// Get image space dimensions from the encoding operator
		//

		boost::shared_ptr< std::vector<unsigned int> > image_dims = this->encoding_operator_->get_domain_dimensions();
		if( image_dims->size() == 0 ){
			throw std::runtime_error("Error: ncgSolver::compute_rhs : encoding operator has not set domain dimension" );
		}

		ARRAY_TYPE * x = new ARRAY_TYPE;
		x->create(image_dims.get());



		ARRAY_TYPE * g = new ARRAY_TYPE;
		g->create(image_dims.get());
		ARRAY_TYPE *  g_old = new ARRAY_TYPE;
		g_old->create(image_dims.get());

		if (this->x0_.get()){
			*x = *(this->x0_.get());
		} else  {
			clear(x);
		}

		REAL rho = REAL(0.5);
		REAL delta = REAL(0.001);
		REAL sigma = REAL(0.4);

		std::vector<ARRAY_TYPE> regEnc;

		for (int i = 0; i < this->regularization_operators_.size(); i++){
			regEnc.push_back(ARRAY_TYPE(this->regularization_operators_[i]->get_codomain_dimensions()));
			if (reg_priors[i].get()){
				regEnc.back() = *reg_priors[i];
				regEnc.back() *= -std::sqrt(this->regularization_operators_[i]->get_weight());
			}

		}
		std::vector<ARRAY_TYPE> regEnc2 = regEnc;

		ARRAY_TYPE d(image_dims.get());
		clear(&d);
		ARRAY_TYPE encoding_space(in->get_dimensions().get());

		ARRAY_TYPE gtmp(image_dims.get());
		ARRAY_TYPE encoding_space2(in->get_dimensions().get());
		REAL reg_res,data_res;

		if( this->output_mode_ >= solver<ARRAY_TYPE,ARRAY_TYPE>::OUTPUT_VERBOSE ){
			std::cout << "Iterating..." << std::endl;
		}
		REAL grad0;
		for (int i = 0; i < iterations_; i++){
			if (i==0){
				if (this->x0_.get()){
					this->encoding_operator_->mult_M(x,&encoding_space);

				} else clear(&encoding_space);
				encoding_space -= *in;
				this->encoding_operator_->mult_MH(&encoding_space,g);

				*g *=  this->encoding_operator_->get_weight();
				data_res = std::sqrt(this->encoding_operator_->get_weight())*real(dot(&encoding_space,&encoding_space));

				calc_regMultM(x,regEnc);
				for (int n = 0; n < regEnc.size(); n++)
					if (reg_priors[n].get())
						axpy(-std::sqrt(this->regularization_operators_[n]->get_weight()),reg_priors[n].get(),&regEnc[n]);
				this->add_gradient(x,g);
				add_linear_gradient(regEnc,g);
				reg_res=REAL(0);

			}else {
				data_res = real(dot(&encoding_space,&encoding_space));
			}


			//REAL grad_norm = nrm2(g);
			REAL grad_norm=REAL(1);
			if (non_negativity_constraint_) solver_non_negativity_filter(x,g);
			
			/*if (i==0) grad0=::abs((*g)[amax(g)]);
			  else grad_norm = ::abs((*g)[amax(g)])/(grad0); //L-inf norm*/
			if (i==0) grad0 = dot(g,g);
			else grad_norm = dot(g,g)/grad0;
			if( this->output_mode_ >= solver<ARRAY_TYPE,ARRAY_TYPE>::OUTPUT_VERBOSE ){

				std::cout << "Iteration " <<i << ". Gradient norm: " <<  grad_norm << std::endl;
				std::cout << "Data residual: " << data_res << std::endl;
			}
			if (i == 0){
				d -= *g;
			} else {
				REAL g_old_norm = dot(g_old,g_old);
				ELEMENT_TYPE ggold = dot(g,g_old);
				*g_old -= *g;
				ELEMENT_TYPE gg = dot(g,g);
				//ELEMENT_TYPE beta = std::max(REAL(0),-dot(g,g_old)/g_old_norm); //PRP ste[
				ELEMENT_TYPE betaDy = -gg/dot(&d,g_old);
				ELEMENT_TYPE betaHS = dot(g,g_old)/dot(&d,g_old);

				ELEMENT_TYPE beta = std::max(REAL(0),std::min(betaDy,betaHS)); //Hybrid step size from Dai and Yuan 2001

				std::cout << "Beta " << beta << std::endl;
				//ELEMENT_TYPE beta(0);

				d *= beta;
				d -= *g;

			}

			this->encoding_operator_->mult_M(&d,&encoding_space2);
			//this->encoding_operator_->mult_MH(&encoding_space2,&gtmp);
			calc_regMultM(&d,regEnc2);

			REAL alpha0 = REAL(1);
			if (this->operators.size() == 0) alpha0 = -(dot(&encoding_space,&encoding_space2)+calc_dot(regEnc,regEnc2))/(dot(&encoding_space2,&encoding_space2)+calc_dot(regEnc2,regEnc2));
			//REAL alpha0 = REAL(1);
			REAL alpha;
			REAL alpha_old;
			REAL old_norm = functionValue(&encoding_space,regEnc,x);

			REAL gd = dot(g,&d);

			*g_old = *g;

			bool wolfe=false;

			int k = 0;

			alpha_old = 0;
			ARRAY_TYPE x2(*x);

			while (not wolfe){

				alpha=alpha0*std::pow(rho,k);
				axpy(alpha-alpha_old,&encoding_space2,&encoding_space);
				reg_axpy(alpha-alpha_old,regEnc2,regEnc);
				axpy(alpha-alpha_old,&d,&x2);

				//axpy(alpha-alpha_old,&gtmp,g);

				if (functionValue(&encoding_space,regEnc,&x2) <= old_norm+alpha*delta*gd) wolfe = true;//Strong Wolfe condition..
				//if ((dot(g,&d)) >= sigma*gd) wolfe =false; //So... officially this is part of the strong Wolfe condition. For semi-linear problems our initial step size should be sufficient.
				k++;
				//std::cout << "Res: " << dot(&encoding_space,&encoding_space)+calc_dot(regEnc,regEnc) << " Target: " << old_norm+alpha*delta*gd << std::endl;
				//				std::cout << "Step2: " << dot(&gdiff,&d) << " Target " << sigma*gd  << std::endl;
				if (alpha == 0) throw std::runtime_error("Wolfe line search failed");
				alpha_old = alpha;
			}

			axpy(alpha,&d,x);
			if (non_negativity_constraint_){
				clamp_min(x,ELEMENT_TYPE(0));
				this->encoding_operator_->mult_M(x,&encoding_space);
				encoding_space -= *in;
				calc_regMultM(x,regEnc);
				for (int n = 0; n < regEnc.size(); n++)
					if (reg_priors[n].get())
						axpy(-std::sqrt(this->regularization_operators_[n]->get_weight()),reg_priors[n].get(),&regEnc[n]);
			}
			this->encoding_operator_->mult_MH(&encoding_space,g);
			this->add_gradient(x,g);
			add_linear_gradient(regEnc,g);


			iteration_callback(x,i,data_res,reg_res);

			
			if (grad_norm < tc_tolerance_)  break;
		}
		delete g,g_old;

		return boost::shared_ptr<ARRAY_TYPE>(x);
	}



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
	// Set preconditioner
	//

	virtual void set_preconditioner( boost::shared_ptr< cgPreconditioner<ARRAY_TYPE> > precond ) {
		precond_ = precond;
	}

	virtual void add_regularization_operator( boost::shared_ptr< linearOperator< ARRAY_TYPE> > op)
	{
		if( !op.get() ){
			throw std::runtime_error( "Error: linearOperatorSolver::add_regularization_operator : NULL operator provided" );
		}
		this->regularization_operators_.push_back(op);
		reg_priors.push_back(boost::shared_ptr<ARRAY_TYPE>((ARRAY_TYPE*)0));
	}

	virtual void add_regularization_operator( boost::shared_ptr< linearOperator< ARRAY_TYPE> > op,boost::shared_ptr<ARRAY_TYPE> prior)
	{
		if( !op.get() ){
			throw std::runtime_error( "Error: linearOperatorSolver::add_regularization_operator : NULL operator provided" );
		}

		this->regularization_operators_.push_back(op);
		reg_priors.push_back(prior);
	}

	virtual void add_regularization_operator(boost::shared_ptr< linearOperator<ARRAY_TYPE> > op, int L_norm ){
		if (L_norm==1){

			this->operators.push_back(boost::shared_ptr< l1GPRegularizationOperator>(new l1GPRegularizationOperator(op)));
		}else{
			add_regularization_operator(op);
		}
	}


	virtual void add_regularization_operator(boost::shared_ptr< linearOperator<ARRAY_TYPE> > op, boost::shared_ptr<ARRAY_TYPE> prior, int L_norm ){
		if (L_norm==1){
			this->operators.push_back(boost::shared_ptr<l1GPRegularizationOperator>(new l1GPRegularizationOperator(op,prior)));
		}else{
			add_regularization_operator(op,prior);
		}
	}


protected:
	typedef typename std::vector<boost::shared_ptr<linearOperator<ARRAY_TYPE> > >::iterator  csIterator;
	typedef typename std::vector< std::vector<boost::shared_ptr<linearOperator<ARRAY_TYPE> > > >::iterator csGroupIterator;

	virtual void solver_non_negativity_filter(ARRAY_TYPE*,ARRAY_TYPE*)=0;
	virtual void iteration_callback(ARRAY_TYPE*,int i,REAL,REAL){};




	ELEMENT_TYPE calc_dot(std::vector<ARRAY_TYPE>& x,std::vector<ARRAY_TYPE>& y){
		ELEMENT_TYPE res(0);
		for (int  i = 0; i < x.size(); i++)
			res += dot(&x[i],&y[i]);
		return res;
	}

	void add_linear_gradient(std::vector<ARRAY_TYPE>& elems, ARRAY_TYPE* g){
		ARRAY_TYPE tmp(g->get_dimensions());
		for (int i = 0; i <elems.size(); i++){
			this->regularization_operators_[i]->mult_MH(&elems[i],&tmp);
			axpy(std::sqrt(this->regularization_operators_[i]->get_weight()),&tmp,g);
		}
	}

	void calc_regMultM(ARRAY_TYPE* x,std::vector<ARRAY_TYPE>& elems){
		for (int i = 0; i <elems.size(); i++){
			this->regularization_operators_[i]->mult_M(x,&elems[i]);
			elems[i] *= std::sqrt(this->regularization_operators_[i]->get_weight());
		}
	}

	void reg_axpy(REAL alpha, std::vector<ARRAY_TYPE>& x, std::vector<ARRAY_TYPE>& y){
		for (int i = 0; i <x.size(); i++){
			axpy(alpha,&x[i],&y[i]);

		}
	}

	REAL functionValue(ARRAY_TYPE* encoding_space,std::vector<ARRAY_TYPE>& regEnc, ARRAY_TYPE * x){
		REAL res= std::sqrt(this->encoding_operator_->get_weight())*dot(encoding_space,encoding_space);

		for (int i = 0; i  < this->operators.size(); i++){
					res += this->operators[i]->magnitude(x);
		}

		res += calc_dot(regEnc,regEnc);
		return res;

	}






protected:

	// Preconditioner
	//boost::shared_ptr< cgPreconditioner<ARRAY_TYPE> > precond_;
	// Maximum number of iterations
	unsigned int iterations_;
	bool non_negativity_constraint_;
	REAL tc_tolerance_;
	REAL threshold;
	bool dump_residual;

	REAL barrier_threshold;
	// Preconditioner

	std::vector<boost::shared_ptr<ARRAY_TYPE> > reg_priors;
	boost::shared_ptr< cgPreconditioner<ARRAY_TYPE> > precond_;

};
}
