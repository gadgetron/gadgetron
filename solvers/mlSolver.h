#pragma once
#include "linearSolver.h"

namespace Gadgetron{
template <class ARRAY_TYPE> class mlSolver
: public linearSolver<ARRAY_TYPE> {
	typedef typename ARRAY_TYPE::element_type ELEMENT_TYPE;
		typedef typename realType<ELEMENT_TYPE>::type REAL;
	public:
		mlSolver(): linearSolver<ARRAY_TYPE>() {
			_iterations=10;
		}
		virtual ~mlSolver(){};

		void set_max_iterations(int i){_iterations=i;}

		boost::shared_ptr<ARRAY_TYPE> solve(ARRAY_TYPE* in){
			//boost::shared_ptr<ARRAY_TYPE> rhs = compute_rhs(in);
			if( this->encoding_operator_.get() == 0 ){
			 BOOST_THROW_EXCEPTION(runtime_error( "Error: mlSolver::solve : no encoding operator is set" ));
			  return boost::shared_ptr<ARRAY_TYPE>();
			}

			// Get image space dimensions from the encoding operator
			//

			boost::shared_ptr< std::vector<unsigned int> > image_dims = this->encoding_operator_->get_domain_dimensions();
			if( image_dims->size() == 0 ){
			 BOOST_THROW_EXCEPTION(runtime_error( "Error: mlSolver::solve : encoding operator has not set domain dimension" ));
			  return boost::shared_ptr<ARRAY_TYPE>();
			}

			ARRAY_TYPE * x = new ARRAY_TYPE;
					x->create(image_dims.get());
			if (this->x0_.get()){
						*x = *(this->x0_.get());
			} else  {
				x->fill(ELEMENT_TYPE(1));
			}

			ARRAY_TYPE tmp_image;
			tmp_image.create(image_dims.get());
			
			ARRAY_TYPE ones_image;
			ones_image.create(image_dims.get());



			ARRAY_TYPE tmp_projection;
			tmp_projection.create(in->get_dimensions().get());
			tmp_projection.fill(ELEMENT_TYPE(1));
			this->encoding_operator_->mult_MH(&tmp_projection,&ones_image,false);

			clamp_min(&ones_image,ELEMENT_TYPE(1e-6));
			ones_image.reciprocal();

			clamp_max(&ones_image,ELEMENT_TYPE(1e6));

			for (int i =0; i < _iterations; i++){
				this->encoding_operator_->mult_M(x,&tmp_projection,false);
				tmp_projection += ELEMENT_TYPE(1);
				tmp_projection.reciprocal();
				tmp_projection *= *in;
				this->encoding_operator_->mult_MH(&tmp_projection,&tmp_image,false);
				tmp_image *= ones_image;
				*x *= tmp_image;
				std::cout << "Iteration: " << i << std::endl;

				
			}

			return boost::shared_ptr<ARRAY_TYPE>(x);
		}
		virtual void add_regularization_operator( boost::shared_ptr< linearOperator<ARRAY_TYPE> > op){
			std::cout << "Error mlSolver::add_regularization_operator - SART Solver does NOT support regularization operators." << std::endl;
		}


protected:
	int _iterations;


};
}
