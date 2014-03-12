#pragma once
#include "subsetOperator.h"
#include "solver.h"
#include <numeric>
#include <vector>
#include <functional>
#include <boost/iterator/counting_iterator.hpp>
#include "gpSolver.h"

namespace Gadgetron{
template <class ARRAY_TYPE> class oscGPBBSolver : public gpSolver< ARRAY_TYPE> {
	typedef typename ARRAY_TYPE::element_type ELEMENT_TYPE;
	typedef typename realType<ELEMENT_TYPE>::Type REAL;
	public:
		oscGPBBSolver() : gpSolver< ARRAY_TYPE>() {
			_iterations=10;
			_beta = REAL(0.8);
			non_negativity_=false;
		}
		virtual ~oscGPBBSolver(){};


		void set_max_iterations(int i){_iterations=i;}
		int get_max_iterations(){return _iterations;}
		void set_non_negativity_constraint(bool neg=true){non_negativity_=neg;}


		boost::shared_ptr<ARRAY_TYPE> solve(ARRAY_TYPE* in){
			//boost::shared_ptr<ARRAY_TYPE> rhs = compute_rhs(in);
			if( this->encoding_operator_.get() == 0 ){
			  throw std::runtime_error( "Error: cgSolver::compute_rhs : no encoding operator is set" );
			  return boost::shared_ptr<ARRAY_TYPE>();
			}

			// Get image space dimensions from the encoding operator
			//

			boost::shared_ptr< std::vector<unsigned int> > image_dims = this->encoding_operator_->get_domain_dimensions();
			if( image_dims->size() == 0 ){
			  throw std::runtime_error( "Error: cgSolver::compute_rhs : encoding operator has not set domain dimension" );
			  return boost::shared_ptr<ARRAY_TYPE>();
			}
			std::ofstream textFile("residual.txt");
			ARRAY_TYPE * x = new ARRAY_TYPE(image_dims);


			ARRAY_TYPE  g(image_dims);
			ARRAY_TYPE g_old(image_dims);
			if (this->x0_.get()) *x = *(this->x0_.get());
			else clear(x);

			ARRAY_TYPE x_old(*x);

			std::vector<boost::shared_ptr<ARRAY_TYPE> > subsets = this->encoding_operator_->projection_subsets(in);

			ARRAY_TYPE encoding_space(in->get_dimensions());
			REAL reg_res,data_res;
			std::vector<boost::shared_ptr<ARRAY_TYPE> > encoding_spaces = this->encoding_operator_->projection_subsets(&encoding_space);
			if( this->output_mode_ >= solver<ARRAY_TYPE,ARRAY_TYPE>::OUTPUT_VERBOSE ){
				std::cout << "oscGPBB setup done, starting iterations:" << std::endl;
			}

			std::vector<int> isubsets(boost::counting_iterator<int>(0), boost::counting_iterator<int>(this->encoding_operator_->get_number_of_subsets()));
			const int m = this->encoding_operator_->get_number_of_subsets();
			for (int i =0; i < _iterations; i++){
				REAL nabla;
				REAL nabla0;

				for (int isubset = 0; isubset < m; isubset++){
					int subset = isubsets[isubset];
					this->encoding_operator_->mult_M(x,encoding_spaces[subset].get(),subset,false);
					*encoding_spaces[subset] -= *subsets[subset];

					if( this->output_mode_ >= solver<ARRAY_TYPE,ARRAY_TYPE>::OUTPUT_VERBOSE ){
									std::cout << "Iteration " <<i << " Subset " << subset << " Update norm: " << nrm2(encoding_spaces[subset].get()) << std::endl;
					}

					this->encoding_operator_->mult_MH(encoding_spaces[subset].get(),&g,subset,false);
					g *= this->encoding_operator_->get_weight();



					this->add_gradient(x,&g);
					if (non_negativity_) solver_non_negativity_filter(x,&g);

					if (isubset == 0){
					if (i<1 ){
						ARRAY_TYPE tmp_encoding(subsets[subset]->get_dimensions());
						this->encoding_operator_->mult_M(&g,&tmp_encoding,i,false);
						if (this->x0_.get()){
							nabla0 = dot(encoding_spaces[i].get(),&tmp_encoding)/dot(&tmp_encoding,&tmp_encoding);
						} else {
							nabla0 = -dot(subsets[subset].get(),&tmp_encoding)/dot(&tmp_encoding,&tmp_encoding);
						}
					}

					 nabla = nabla0/(REAL(i)/15+1);
					}

					axpy(-nabla,&g,x);
					if (non_negativity_) clamp_min(x,ELEMENT_TYPE(0));


				}
				std::reverse(isubsets.begin(),isubsets.end());
				//DEBUG DEMON CODE BEGINS HERE
				ARRAY_TYPE tmp_proj(*in);
				clear(&tmp_proj);
				this->encoding_operator_->mult_M(x,&tmp_proj,false);
				tmp_proj -= *in;
				REAL residual = dot(&tmp_proj,&tmp_proj);
				textFile << residual << std::endl;

				iteration_callback(x,i);

			}


			return boost::shared_ptr<ARRAY_TYPE>(x);
		}

		void set_encoding_operator(boost::shared_ptr<subsetOperator<ARRAY_TYPE> > encoding_operator){ encoding_operator_ = encoding_operator; }


protected:
	virtual void solver_non_negativity_filter(ARRAY_TYPE*,ARRAY_TYPE*)=0;
	virtual void iteration_callback(ARRAY_TYPE*, int i){};
	int _iterations;
	REAL _beta;
	bool non_negativity_;
	boost::shared_ptr<subsetOperator<ARRAY_TYPE> > encoding_operator_;

};
}
