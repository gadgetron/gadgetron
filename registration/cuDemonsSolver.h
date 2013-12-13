#pragma once
#include "cuNDArray_math.h"
#include "cuNDArray_utils.h"
#include "multiresRegistrationSolver.h"
#include "cuGaussianFilterOperator.h"
#include "vector_td.h"
#include "cuNDArray.h"

namespace Gadgetron{
template<class T, unsigned int D> class cuDemonsSolver : public multiresRegistrationSolver<cuNDArray<T>, D>{


public:
	cuDemonsSolver() : alpha(1.0),beta(1e-6),sigma(5.0){};
	virtual ~cuDemonsSolver(){};

	virtual void compute( cuNDArray<T> *fixed_image, cuNDArray<T> *moving_image, cuNDArray<T> *stencil_image, boost::shared_ptr<cuNDArray<T> > &result ){


		cuNDArray<T> def_moving(*moving_image);

		cuGaussianFilterOperator<T,D> gauss;
		gauss.set_sigma(sigma);
		std::vector<unsigned int> dims = *moving_image->get_dimensions();

		dims.push_back(D);

		if (!result.get()){
			result = boost::shared_ptr<cuNDArray<T> >(new cuNDArray<T>(&dims));
			clear(result.get());
		}

		for (int i = 0; i < this->max_num_iterations_per_level_; i++){
			//Calculate the gradients

			boost::shared_ptr<cuNDArray<T> > update = demonicStep(fixed_image,&def_moving);
			if (sigma > 0){
				cuNDArray<T> blurred_update(update->get_dimensions());
				gauss.mult_M(update.get(),&blurred_update);
				std::cout << "Update step: " << nrm2(&blurred_update) << std::endl;

				*result -= blurred_update;

			} else *result -= *update;




			def_moving = *this->deform(moving_image,result);

		}






	}

	void set_sigma(T _sigma){
		sigma = _sigma;

			}
	void set_alpha(T _alpha){
		alpha = _alpha;
	}



protected:
	boost::shared_ptr<cuNDArray<T> > demonicStep(cuNDArray<T>*,cuNDArray<T>*);


	T sigma,alpha,beta;
};
}
