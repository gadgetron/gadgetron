/*
 * CSIOperator.h
 *
 *  Created on: Nov 10, 2014
 *      Author: dch
 */

#ifndef CSIOPERATOR_H_
#define CSIOPERATOR_H_

#include "linearOperator.h"
#include "cuNDArray.h"

namespace Gadgetron {

template <class T> class CSIOperator: public Gadgetron::linearOperator<cuNDArray<complext<T>>> {
public:
	CSIOperator();
	virtual ~CSIOperator();
	virtual void mult_M(cuNDArray<complext<T>>* in, cuNDArray<complext<T>>* out,bool accumulate );
	virtual void mult_MH(cuNDArray<complext<T>>* in, cuNDArray<complext<T>>* out,bool accumulate );

	void set_senseOp(boost::shared_ptr<linearOperator<cuNDArray<complext<T>>>> op){ senseOp = op;}

	virtual boost::shared_ptr<linearOperator<cuNDArray<complext<T>>>> clone(){
		return linearOperator<cuNDArray<complext<T>>>::clone(this);
	}
protected:
	boost::shared_ptr<linearOperator<cuNDArray<complext<T>>>> senseOp;

};

} /* namespace Gadgetron */
#endif /* CSIOPERATOR_H_ */
