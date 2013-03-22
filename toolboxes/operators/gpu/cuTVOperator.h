#pragma once
#include "cuNDArray.h"
#include "generalOperator.h"
#include "complext.h"
namespace Gadgetron{

template<class T, unsigned int D> class cuTVOperator : public generalOperator<cuNDArray<T> > {
private:
	typedef typename realType<T>::type REAL;
public:
	cuTVOperator() : generalOperator<cuNDArray<T> >(){
		limit_ = REAL(1e-8);
	}
	virtual ~cuTVOperator(){};
	void set_limit(REAL limit){
		limit_ = limit;
	}
	virtual void gradient(cuNDArray<T>*,cuNDArray<T>*, bool accumulate=false);
protected:
	REAL limit_;
};
}
