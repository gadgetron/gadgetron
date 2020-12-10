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
	CSIOperator(T dtt, T dte);
	virtual ~CSIOperator();
	virtual void mult_M(cuNDArray<complext<T>>* in, cuNDArray<complext<T>>* out,bool accumulate );
	virtual void mult_MH(cuNDArray<complext<T>>* in, cuNDArray<complext<T>>* out,bool accumulate );

	void set_senseOp(boost::shared_ptr<linearOperator<cuNDArray<complext<T>>>> op){ senseOp = op;}
	void set_frequencies(std::vector<T>& freq);
	
	cuNDArray<T> get_frequencies(){
		return frequencies;
	}


	T get_echotime(){ return dte_;}
	T get_pointtime(){return dtt_;}

protected:
	boost::shared_ptr<linearOperator<cuNDArray<complext<T>>>> senseOp;
	T dte_; //Time between echoes
	T dtt_; //Time between k-space points
	cuNDArray<T> frequencies;
};

} /* namespace Gadgetron */

#endif /* CSIOPERATOR_H_ */
