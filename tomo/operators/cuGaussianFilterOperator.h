#pragma once
#include "linearOperator.h"
#include "cuNDArray.h"
#include "cuNDArray_math.h"
namespace Gadgetron{


template<class T, unsigned int D> class cuGaussianFilterOperator : public linearOperator<cuNDArray<T> > {
	typedef typename realType<T>::Type REAL;
public:
	void mult_M( cuNDArray<T>* in, cuNDArray<T>* out, bool accumulate = false);
	void mult_MH( cuNDArray<T>* in, cuNDArray<T>* out, bool accumulate = false);
	virtual boost::shared_ptr< linearOperator<cuNDArray<T> > > clone()
		{
			return linearOperator< cuNDArray<T> >::clone(this);
		}
	void set_sigma(REAL sigma){ _sigma=sigma;}
protected:

	REAL _sigma;


};

}
