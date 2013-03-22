#pragma once

#include "picsOperator.h"
#include "hoNDArray.h"
#include "hoNDArray_vector_td_utilities.h"
template<class REAL, class T, unsigned int D> class hoTVpicsOperator : public picsOperator<REAL,T,hoNDArray<T> >{
public:
	hoTVpicsOperator():picsOperator<REAL,T,hoNDArray<T> >(){
		set_CSOperator(boost::shared_ptr<sparsifyingOperator<REAL,hoNDArray<T> > >(new hoTVOperator<REAL,T,D>));
	}
	virtual ~hoTVpicsOperator(){}

	void set_limit(REAL limit){
		hoTVOperator<REAL,T,D>* op = dynamic_cast<hoTVOperator<REAL,T,D>* >( this->csOperator.get());
		op->set_limit(limit);
	}
protected:
	virtual bool pics_axpy(T a,hoNDArray<T>* x,hoNDArray<T>* y){
		return hoNDA_axpy<T>(a,x,y);

	}
	virtual bool pics_scal(T a,hoNDArray<T>* x){
			return hoNDA_scal(a,x);

	}
};
