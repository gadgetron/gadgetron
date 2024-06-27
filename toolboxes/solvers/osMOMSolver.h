/*
 * osMOMSolver.h
 *
 *  Created on: Mar 23, 2015
 *      Author: u051747
 */
//Based on Donghwan Kim; Ramani, S.; Fessler, J.A., "Combining Ordered Subsets and Momentum for Accelerated X-Ray CT Image Reconstruction,"
#pragma once
#include "subsetOperator.h"
#include "solver.h"
#include <numeric>
#include <vector>
#include <functional>
#include <boost/iterator/counting_iterator.hpp>
#include <boost/make_shared.hpp>

namespace Gadgetron{
template <class ARRAY_TYPE> class osMOMSolver : public solver< ARRAY_TYPE,ARRAY_TYPE> {
	typedef typename ARRAY_TYPE::element_type ELEMENT_TYPE;
	typedef typename realType<ELEMENT_TYPE>::Type REAL;
public:
	osMOMSolver() :solver< ARRAY_TYPE,ARRAY_TYPE>() {
		_iterations=10;
		_beta = REAL(1);
		_alpha = 0.2;
		_gamma = 0;
		non_negativity_=false;
		reg_steps_=2;
		_kappa = REAL(1);
	}
	virtual ~osMOMSolver(){};

	void set_max_iterations(int i){_iterations=i;}
	int get_max_iterations(){return _iterations;}
	void set_non_negativity_constraint(bool neg=true){non_negativity_=neg;}
	/**
	 * @brief Sets the weight of each step in the SART iteration
	 * @param beta
	 */
	void set_beta(REAL beta){_beta = beta;}
	void set_gamma(REAL gamma){_gamma = gamma;}
	void set_kappa(REAL kappa){_kappa = kappa;}

	/**
	 * Sets the preconditioning image. In most cases this is not needed, and the preconditioning is calculated based on the system transform
	 * @param precon_image
	 */
	void set_preconditioning_image(boost::shared_ptr<ARRAY_TYPE> precon_image){
		this->preconditioning_image_ = precon_image;
	}


	void set_reg_steps(unsigned int reg_steps){ reg_steps_ = reg_steps;}

	boost::shared_ptr<ARRAY_TYPE> solve(ARRAY_TYPE* in){
		//boost::shared_ptr<ARRAY_TYPE> rhs = compute_rhs(in);
		if( this->encoding_operator_.get() == 0 ){
			throw std::runtime_error( "Error: cgSolver::compute_rhs : no encoding operator is set" );
			return boost::shared_ptr<ARRAY_TYPE>();
		}

		// Get image space dimensions from the encoding operator
		//

		boost::shared_ptr< std::vector<size_t> > image_dims = this->encoding_operator_->get_domain_dimensions();
		if( image_dims->size() == 0 ){
			throw std::runtime_error( "Error: cgSolver::compute_rhs : encoding operator has not set domain dimension" );
			return boost::shared_ptr<ARRAY_TYPE>();
		}

		ARRAY_TYPE * z = new ARRAY_TYPE(*image_dims);
		if (this->x0_.get()){
			*z = *(this->x0_.get());
		} else  {
			clear(z);
		}

		ARRAY_TYPE * x = new ARRAY_TYPE(*z);
		ARRAY_TYPE * xold = new ARRAY_TYPE(*z);
		GINFO("DEBUG DINGO 1\n");
		std::vector<boost::shared_ptr<ARRAY_TYPE> > subsets = this->encoding_operator_->projection_subsets(in);

		GINFO("DEBUG DINGO 2\n");
		ARRAY_TYPE tmp_projection(*in->get_dimensions());
		std::vector<boost::shared_ptr<ARRAY_TYPE> > tmp_projections = this->encoding_operator_->projection_subsets(&tmp_projection);

		GINFO("DEBUG DINGO 3\n");
		boost::shared_ptr<ARRAY_TYPE> precon_image;
		if (preconditioning_image_)
			precon_image = preconditioning_image_;
		else {
			precon_image = boost::make_shared<ARRAY_TYPE>(*image_dims);
			fill(precon_image.get(),ELEMENT_TYPE(1));
			this->encoding_operator_->mult_M(precon_image.get(),&tmp_projection,false);
			this->encoding_operator_->mult_MH(&tmp_projection,precon_image.get(),false);
			clamp_min(precon_image.get(),REAL(1e-6));
			reciprocal_inplace(precon_image.get());
			//ones_image *= (ELEMENT_TYPE) this->encoding_operator_->get_number_of_subsets();
		}
		ARRAY_TYPE tmp_image(*image_dims);


		GINFO("DEBUG DINGO 4\n");
		REAL t = 1;
		REAL told = 1;
		if( this->output_mode_ >= solver<ARRAY_TYPE,ARRAY_TYPE>::OUTPUT_VERBOSE ){
			GINFO( "osMOM setup done, starting iterations:\n");
		}

		std::vector<int> isubsets(boost::counting_iterator<int>(0), boost::counting_iterator<int>(this->encoding_operator_->get_number_of_subsets()));
		REAL kappa_int = _kappa;
		REAL step_size;
		for (int i =0; i < _iterations; i++){
			for (int isubset = 0; isubset < this->encoding_operator_->get_number_of_subsets(); isubset++){

				t = 0.5*(1+std::sqrt(1+4*t*t));
				int subset = isubsets[isubset];
				this->encoding_operator_->mult_M(x,tmp_projections[subset].get(),subset,false);
				*tmp_projections[subset] -= *subsets[subset];
				*tmp_projections[subset] *= ELEMENT_TYPE(-1);
				if( this->output_mode_ >= solver<ARRAY_TYPE,ARRAY_TYPE>::OUTPUT_VERBOSE ){
					GINFO_STREAM("Iteration %d Subset %d " << i << " Subset " << subset << " Update norm: " << nrm2(tmp_projections[subset].get()) << std::endl);
				}

				this->encoding_operator_->mult_MH(tmp_projections[subset].get(),&tmp_image,subset,false);
				tmp_image *= *precon_image;
				axpy(REAL(_beta/(1+_gamma*i))*this->encoding_operator_->get_number_of_subsets(),&tmp_image,x);
				if (i ==0 && isubset == 0)
					step_size = _alpha*nrm2(x)/this->encoding_operator_->get_number_of_subsets();


				for (auto op : regularization_operators){
					/*
									for (auto i = 0u; i < reg_steps_; i++){
										op->gradient(x,&tmp_image);
										tmp_image *= REAL(1)/nrm2(&tmp_image);
										axpy(-step_size*op->get_weight(),&tmp_image,x);
									}
					 */


					GINFO("Reg: %d\n", op->magnitude(x));

									tmp_image /= nrm2(&tmp_image);
									auto reg_val = op->magnitude(x);
									GINFO("Reg val: %f\n", reg_val);
									ARRAY_TYPE y = *x;
									axpy(-kappa_int,&tmp_image,&y);


									while(op->magnitude(&y) > reg_val){

										kappa_int /= 2;
										axpy(kappa_int,&tmp_image,&y);
										GINFO("Kappa: %d\n", kappa_int);
									}

									reg_val = op->magnitude(&y);
					 *x = y;

				}

				//axpy(REAL(_beta),&tmp_image,x);
				if (non_negativity_){
					clamp_min(x,REAL(0));
				}


				*z = *x;
				*z *= 1+(told-1)/t;
				axpy(-(told-1)/t,xold,z);
				std::swap(x,xold);
				*x = *z;

				told = t;


				//step_size *= 0.99;

			}
			//std::reverse(isubsets.begin(),isubsets.end());
			//std::random_shuffle(isubsets.begin(),isubsets.end());
			/*
			ARRAY_TYPE tmp_proj(*in);
			clear(&tmp_proj);
			this->encoding_operator_->mult_M(x,&tmp_proj,false);
			tmp_proj -= *in;


			std::stringstream ss;
			ss << "osMOM-" << i << ".real";

			write_nd_array<ELEMENT_TYPE>(x,ss.str().c_str());

			//calc_regMultM(x,regEnc);
			//REAL f = functionValue(&tmp_proj,regEnc,x);
			GINFO_STREAM("Function value: " << dot(&tmp_proj,&tmp_proj) << std::endl);
			 */
		}
		delete x,xold;


		return boost::shared_ptr<ARRAY_TYPE>(z);
	}

	void set_encoding_operator(boost::shared_ptr<subsetOperator<ARRAY_TYPE> > encoding_operator){ encoding_operator_ = encoding_operator; }
	virtual void add_nonlinear_operator(boost::shared_ptr< generalOperator<ARRAY_TYPE> > op ){
		regularization_operators.push_back(op);
	}


protected:
	int _iterations;
	REAL _beta, _gamma, _alpha, _kappa;
	bool non_negativity_;
	unsigned int reg_steps_;
	boost::shared_ptr<subsetOperator<ARRAY_TYPE> > encoding_operator_;
	std::vector<boost::shared_ptr<generalOperator<ARRAY_TYPE>>> regularization_operators;
	boost::shared_ptr<ARRAY_TYPE> preconditioning_image_;

};
}
