#pragma once
#include "subsetOperator.h"
#include "GadgetronException.h"
#include "solver.h"
#include <numeric>
#include <vector>
#include <functional>
#include <boost/iterator/counting_iterator.hpp>

namespace Gadgetron{
template <class ARRAY_TYPE> class osGPBBSolver : public solver< ARRAY_TYPE,ARRAY_TYPE> {
	typedef typename ARRAY_TYPE::element_type ELEMENT_TYPE;
	typedef typename realType<ELEMENT_TYPE>::Type REAL;
	public:
		osGPBBSolver() :solver< ARRAY_TYPE,ARRAY_TYPE>() {
			_iterations=10;
			_beta = REAL(0.8);
			non_negativity_=false;
		}
		virtual ~osGPBBSolver(){};

		void set_max_iterations(int i){_iterations=i;}
		int get_max_iterations(){return _iterations;}
		void set_non_negativity_constraint(bool neg=true){non_negativity_=neg;}
		/**
		 * @brief Sets the weight of each step in the SART iteration
		 * @param beta
		 */
		void set_beta(REAL beta){_beta = beta;}

		boost::shared_ptr<ARRAY_TYPE> solve(ARRAY_TYPE* in){
			//boost::shared_ptr<ARRAY_TYPE> rhs = compute_rhs(in);
			if( this->encoding_operator_.get() == 0 ){
			 BOOST_THROW_EXCEPTION(runtime_error( "Error: cgSolver::compute_rhs : no encoding operator is set" ));
			  return boost::shared_ptr<ARRAY_TYPE>();
			}

			// Get image space dimensions from the encoding operator
			//

			boost::shared_ptr< std::vector<unsigned int> > image_dims = this->encoding_operator_->get_domain_dimensions();
			if( image_dims->size() == 0 ){
			 BOOST_THROW_EXCEPTION(runtime_error( "Error: cgSolver::compute_rhs : encoding operator has not set domain dimension" ));
			  return boost::shared_ptr<ARRAY_TYPE>();
			}

			ARRAY_TYPE * x = new ARRAY_TYPE(image_dims);
			ARRAY_TYPE x_old(image_dims);

			ARRAY_TYPE * g = new ARRAY_TYPE(image_dims);
			ARRAY_TYPE *  g_old = new ARRAY_TYPE(image_dims);
			if (this->x0_.get()) *x = *(this->x0_.get());
			else clear(x);

			std::vector<boost::shared_ptr<ARRAY_TYPE> > subsets = this->encoding_operator_->projection_subsets(in);

			ARRAY_TYPE encoding_space(in->get_dimensions());
			REAL reg_res,data_res;
			std::vector<boost::shared_ptr<ARRAY_TYPE> > encoding_spaces = this->encoding_operator_->projection_subsets(&encoding_space);
			if( this->output_mode_ >= solver<ARRAY_TYPE,ARRAY_TYPE>::OUTPUT_VERBOSE ){
				std::cout << "osGPBB setup done, starting iterations:" << std::endl;
			}

			for (int i =0; i < _iterations; i++){
				for (int subset = 0; subset < this->encoding_operator_->get_number_of_subsets(); subset++){
					this->encoding_operator_->mult_M(x,encoding_spaces[subset].get(),subset,false);
					*encoding_spaces[subset] -= *subsets[subset];

					if( this->output_mode_ >= solver<ARRAY_TYPE,ARRAY_TYPE>::OUTPUT_VERBOSE ){
									std::cout << "Iteration " <<i << " Subset " << subset << " Update norm: " << nrm2(encoding_spaces[subset].get()) << std::endl;
					}

					this->encoding_operator_->mult_MH(encoding_spaces[subset].get(),g,subset,false);
					*g *= this->encoding_operator_->get_weight();


					//this->add_gradient(x,g);
					if (non_negativity_) solver_non_negativity_filter(x,g);

					REAL nabla;
					if (i<1 && subset == 0){
						ARRAY_TYPE tmp_encoding(subsets[subset]->get_dimensions());
						this->encoding_operator_->mult_M(g,&tmp_encoding,i,false);
						if (this->x0_.get()){
							nabla = dot(encoding_spaces[i].get(),&tmp_encoding)/dot(&tmp_encoding,&tmp_encoding);
						} else {
							nabla = -dot(subsets[subset].get(),&tmp_encoding)/dot(&tmp_encoding,&tmp_encoding);
						}
					} else {
						x_old -= *x;
						*g_old -= *g;
						REAL xx = dot(&x_old,&x_old);
						REAL gx = dot(g_old,&x_old);

						REAL nabla1 = xx/gx;

						/* This is the code that enables the adaptive step size.
							 REAL gg = dot(g_old,&x_old);
							 REAL nabla2 = gx/gg;
							 if ((nabla2/nabla1) < 0.5) nabla = nabla2;
							 else nabla = nabla1;*/
						nabla = nabla1;
					}
					ARRAY_TYPE * tmp;
					tmp=g_old;
					g_old=g;
					g=tmp;

					x_old = *x;
					axpy(-nabla,g_old,x);
					if (non_negativity_) clamp_min(x,ELEMENT_TYPE(0));
				}
			}
			delete g,g_old;

			return boost::shared_ptr<ARRAY_TYPE>(x);
		}

		void set_encoding_operator(boost::shared_ptr<subsetOperator<ARRAY_TYPE> > encoding_operator){ encoding_operator_ = encoding_operator; }


protected:
	virtual void solver_non_negativity_filter(ARRAY_TYPE*,ARRAY_TYPE*)=0;
	int _iterations;
	REAL _beta;
	bool non_negativity_;
	boost::shared_ptr<subsetOperator<ARRAY_TYPE> > encoding_operator_;

};
}
