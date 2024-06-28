#pragma once

#include "opticalFlowOperator.h"
#include "cuPartialDerivativeOperator.h"
#include "cuCgSolver.h"
namespace Gadgetron{


template<class T, unsigned int D> class cuCGHSOFSolver : public multiresRegistrationSolver<cuNDArray<T>, D>{
public:
	cuCGHSOFSolver(){

		OF = boost::shared_ptr<OFOp>(new OFOp);
		solver = boost::shared_ptr< cuCgSolver<T> >(new cuCgSolver<T>);
		solver->set_encoding_operator(OF);
		for (unsigned int i = 0; i < D; i++){
			boost::shared_ptr<cuPartialDerivativeOperator<T,D> > dx(new cuPartialDerivativeOperator<T,D>(i));
			solver->add_regularization_operator(dx);
			ops.push_back(dx);
		}
	}

	virtual ~cuCGHSOFSolver(){};
	typedef opticalFlowOperator<cuNDArray<T>,cuPartialDerivativeOperator<T,D>,D> OFOp;

	virtual void compute( cuNDArray<T> *fixed_image, cuNDArray<T> *moving_image, cuNDArray<T> *stencil_image, boost::shared_ptr<cuNDArray<T> > &result )
  {
		std::vector<size_t> dims = fixed_image->get_dimensions();
		OF->set_codomain_dimensions(dims);    
		OF->set_images(fixed_image,moving_image);

		for (int i = 0; i < ops.size(); i++){
				ops[i]->set_domain_dimensions(dims);
				ops[i]->set_codomain_dimensions(dims);
				ops[i]->set_weight(_alpha);
		}

		dims.push_back(D);
		OF->set_domain_dimensions(dims);
		cuNDArray<T> It(*fixed_image);
		It -= *moving_image;
		boost::shared_ptr<cuNDArray<T> > resOp = solver->solve(&It);

		if (result.get()) *result += *resOp;
		else result = resOp;
	}

	void set_alpha(T alpha){
		_alpha = alpha;
	}

	boost::shared_ptr< cuCgSolver<T> > get_solver(){
		return solver;
	}

protected:

	T _alpha;
	boost::shared_ptr< cuCgSolver<T> > solver;
	boost::shared_ptr<OFOp> OF;
	std::vector<boost::shared_ptr<cuPartialDerivativeOperator<T,D> >  >ops;
};




}
