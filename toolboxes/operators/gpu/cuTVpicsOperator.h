#pragma once


#include "cuNDArray.h"
#include "cuTVOperator.h"

namespace Gadgetron{


template<class T, unsigned int D> class cuTVpicsOperator : public generalOperator<cuNDArray<T> >{
private:
	typedef typename realType<T>::Type REAL;
public:
	cuTVpicsOperator():generalOperator<cuNDArray<T> >(){

	}
	void set_prior(boost::shared_ptr<cuNDArray<T> > _prior){
		prior = _prior;
	}

	virtual void gradient(cuNDArray<T>* in, cuNDArray<T>* out, bool accumulate){
		cuNDArray<T> tmp = *in;
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

	boost::shared_ptr<cuNDArray<T> > get_prior(){
		return prior;
	}
	virtual ~cuTVpicsOperator(){}

protected:
	cuTVOperator<T,D> op;
	boost::shared_ptr<cuNDArray<T> > prior;
};

}
