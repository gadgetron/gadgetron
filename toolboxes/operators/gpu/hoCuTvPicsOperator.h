#pragma once


#include "hoCuNDArray.h"
#include "hoCuTVOperator.h"

namespace Gadgetron{


template<class T, unsigned int D> class hoCuTVPicsOperator : public generalOperator<hoCuNDArray<T> >{
private:
	typedef typename realType<T>::Type REAL;
public:
	hoCuTVPicsOperator():generalOperator<hoCuNDArray<T> >(){

	}
	void set_prior(boost::shared_ptr<hoCuNDArray<T> > _prior){
		prior = _prior;
	}

	virtual void gradient(hoCuNDArray<T>* in, hoCuNDArray<T>* out, bool accumulate){
		hoCuNDArray<T> tmp = *in;
		tmp -= *prior;
		op.gradient(&tmp,out, accumulate);
	}

	void set_limit(REAL limit){
		op.set_limit(limit);
	}

	virtual void set_weight(REAL weight){
		op.set_weight(weight);
	}
	virtual REAL get_weight(){
		return op.get_weight();
	}
	virtual ~hoCuTVPicsOperator(){}

protected:
	hoCuTVOperator<T,D> op;
	boost::shared_ptr<hoCuNDArray<T> > prior;
};

}
