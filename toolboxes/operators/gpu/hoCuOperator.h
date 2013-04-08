#pragma once
#include "hoNDArray_operators.h"
#include "linearOperator.h"

#include <boost/shared_ptr.hpp>
namespace Gadgetron{


template<class T > class hoCuOperator : public linearOperator<hoCuNDArray<T> > {

	public:
		hoCuOperator(){};
		hoCuOperator(boost::shared_ptr<linearOperator<hoNDArray<T> > > _op): op(_op) {};
		virtual ~hoCuOperator(){};

		virtual void mult_M(hoCuNDArray<T>* in, hoCuNDArray<T>* out, bool accumulate=false){
			op->mult_M(in,out,accumulate);
		}
		virtual void mult_MH(hoCuNDArray<T>* in, hoCuNDArray<T>* out, bool accumulate=false){
			op->mult_MH(in,out,accumulate);
		}
		virtual void mult_MH_M(hoCuNDArray<T>* in, hoCuNDArray<T>* out, bool accumulate=false){
			op->mult_MH_M(in,out,accumulate);
		}
		boost::shared_ptr<linearOperator<hoNDArray<T> > > op;
	 virtual boost::shared_ptr< linearOperator< hoCuNDArray<T> > > clone() {
				return linearOperator< hoCuNDArray<T> >::clone(this);
			}
};

template<class T> boost::shared_ptr<linearOperator<hoCuNDArray<T> > > to_hoCu(boost::shared_ptr<linearOperator<hoNDArray<T> > > _op){
	return boost::shared_ptr<hoCuOperator<T> > (new hoCuOperator<T>(_op));
}

}
