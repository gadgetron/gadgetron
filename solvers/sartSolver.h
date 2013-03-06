#pragma once
#include "linearSolver.h"

namespace Gadgetron{
template <class ARRAY_TYPE> class sartSolver
: public linearSolver<ARRAY_TYPE> {
	typedef typename ARRAY_TYPE::element_type ELEMENT_TYPE;
	typedef typename realType<ELEMENT_TYPE>::type REAL;
	public:
		sartSolver(): linearSolver<ARRAY_TYPE>() {
			_iterations=10;
			_beta = REAL(0.8);
			non_negativity_=false;
		}
		virtual ~sartSolver(){};

		void set_max_iterations(int i){_iterations=i;}
		void set_non_negativity_constraint(bool neg){non_negativity_=neg;}

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
				x->clear();
			}



			ARRAY_TYPE ones_projection(in->get_dimensions().get());

			ARRAY_TYPE tmp_image(image_dims.get());

			tmp_image.fill(ELEMENT_TYPE(1));
			this->encoding_operator_->mult_M(&tmp_image,&ones_projection,false);

			clamp_min(&ones_projection,ELEMENT_TYPE(1e-6));
			ones_projection.reciprocal();



			ARRAY_TYPE ones_image(image_dims.get());


			ARRAY_TYPE tmp_projection(in->get_dimensions().get());
			tmp_projection.fill(ELEMENT_TYPE(1));
			this->encoding_operator_->mult_MH(&tmp_projection,&ones_image,false);
			clamp_min(&ones_image,ELEMENT_TYPE(1e-6));
			ones_image.reciprocal();

			for (int i =0; i < _iterations; i++){
				this->encoding_operator_->mult_M(x,&tmp_projection,false);
				axpy(ELEMENT_TYPE(-1),in,&tmp_projection);
				tmp_projection *= ELEMENT_TYPE(-1);
				tmp_projection *= ones_projection;
				this->encoding_operator_->mult_MH(&tmp_projection,&tmp_image,false);
				tmp_image *= ones_image;
				axpy(ELEMENT_TYPE(1),&tmp_image,x);
				if( this->output_mode_ >= solver<ARRAY_TYPE,ARRAY_TYPE>::OUTPUT_VERBOSE ){
							  std::cout << "Iteration " <<i << " Update norm: " << nrm2(&tmp_image) << std::endl;
				}
				if (non_negativity_){
					clamp_min(x,ELEMENT_TYPE(0));
				}
			}


			return boost::shared_ptr<ARRAY_TYPE>(x);
		}
		virtual bool add_regularization_operator( boost::shared_ptr< linearOperator<ARRAY_TYPE> > op){
			std::cout << "Error sartSolver::add_regularization_operator - SART Solver does NOT support regularization operators." << std::endl;
			return false;
		}


protected:
	int _iterations;
	REAL _beta;
	bool non_negativity_;
};
}
