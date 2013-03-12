#pragma once
#include "cuNDArray.h"
#include "generalOperator.h"
#include "complext.h"
namespace Gadgetron{
template<class T, unsigned int D> class cuTV1DOperator : public generalOperator< cuNDArray<T> >{
private:
	typedef typename realType<T>::type REAL;
public:
	cuTV1DOperator() : generalOperator< cuNDArray<T> >(){
		limit_ = REAL(1e-8);

	}
	virtual ~cuTV1DOperator(){};
	void set_limit(REAL limit){
		limit_ = limit;
	}
	virtual void gradient(cuNDArray<T>*,cuNDArray<T>*, bool accumulate=false);
protected:
	REAL limit_;

};
}
