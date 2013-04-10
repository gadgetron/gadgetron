#pragma once
#include "subsetOperator.h"
#include "GadgetronException.h"
#include "solver.h"
#include <numeric>
#include <vector>
#include <functional>
#include <boost/iterator/counting_iterator.hpp>

namespace Gadgetron{
template <class ARRAY_TYPE> class osSARTSolver : public solver< ARRAY_TYPE,ARRAY_TYPE> {
	typedef typename ARRAY_TYPE::element_type ELEMENT_TYPE;
	typedef typename realType<ELEMENT_TYPE>::Type REAL;
	public:
		osSARTSolver() :solver< ARRAY_TYPE,ARRAY_TYPE>() {
			_iterations=10;
			_beta = REAL(0.8);
			non_negativity_=false;
		}
		virtual ~osSARTSolver(){};

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

			ARRAY_TYPE * x = new ARRAY_TYPE;
					x->create(image_dims.get());
			if (this->x0_.get()){
						*x = *(this->x0_.get());
			} else  {
				clear(x);
			}

			std::vector<boost::shared_ptr<ARRAY_TYPE> > subsets = this->encoding_operator_->projection_subsets(in);

			ARRAY_TYPE ones_projection(in->get_dimensions().get());
			ARRAY_TYPE tmp_image(image_dims.get());
			tmp_image.fill(ELEMENT_TYPE(1));
			this->encoding_operator_->mult_M(&tmp_image,&ones_projection,false);
			clamp_min(&ones_projection,ELEMENT_TYPE(1e-6));
			reciprocal_inplace(&ones_projection);

			std::vector<boost::shared_ptr<ARRAY_TYPE> > ones_projections = this->encoding_operator_->projection_subsets(&ones_projection);

			ARRAY_TYPE tmp_projection(in->get_dimensions());
			tmp_projection.fill(ELEMENT_TYPE(1));
			std::vector<boost::shared_ptr<ARRAY_TYPE> > tmp_projections = this->encoding_operator_->projection_subsets(&tmp_projection);

			std::vector<ARRAY_TYPE> ones_images;

			for (int i = 0; i < this->encoding_operator_->get_number_of_subsets(); i++){
				ones_images.push_back(ARRAY_TYPE(image_dims.get()));
				this->encoding_operator_->mult_MH(tmp_projections[i].get(),&ones_images.back(),i,false);
				clamp_min(&ones_images.back(),ELEMENT_TYPE(1e-6));
				reciprocal_inplace(&ones_images.back());
			}
			if( this->output_mode_ >= solver<ARRAY_TYPE,ARRAY_TYPE>::OUTPUT_VERBOSE ){
				std::cout << "osSART setup done, starting iterations:" << std::endl;
			}

			for (int i =0; i < _iterations; i++){
				for (int subset = 0; subset < this->encoding_operator_->get_number_of_subsets(); subset++){
					this->encoding_operator_->mult_M(x,tmp_projections[subset].get(),subset,false);
					*tmp_projections[subset] -= *subsets[subset];
					*tmp_projections[subset] *= ELEMENT_TYPE(-1);
					if( this->output_mode_ >= solver<ARRAY_TYPE,ARRAY_TYPE>::OUTPUT_VERBOSE ){
									std::cout << "Iteration " <<i << " Subset " << subset << " Update norm: " << nrm2(tmp_projections[subset].get()) << std::endl;
					}

					*tmp_projections[subset] *= *ones_projections[subset];
					this->encoding_operator_->mult_MH(tmp_projections[subset].get(),&tmp_image,subset,false);
					tmp_image *= ones_images[subset];
					*x += tmp_image;
					if (non_negativity_){
						clamp_min(x,ELEMENT_TYPE(0));
					}

				}
			}


			return boost::shared_ptr<ARRAY_TYPE>(x);
		}

		void set_encoding_operator(boost::shared_ptr<subsetOperator<ARRAY_TYPE> > encoding_operator){ encoding_operator_ = encoding_operator; }


protected:
	int _iterations;
	REAL _beta;
	bool non_negativity_;
	boost::shared_ptr<subsetOperator<ARRAY_TYPE> > encoding_operator_;

};
}
