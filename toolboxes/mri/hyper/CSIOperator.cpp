/*
 * CSIOperator.cpp
 *
 *  Created on: Nov 10, 2014
 *      Author: dch
 */

#include "CSIOperator.h"
#include "cuNDFFT.h"
#include "cuNDArray_math.h"
namespace Gadgetron {

template<class T> CSIOperator<T>::CSIOperator() {
	// TODO Auto-generated constructor stub

}

template<class T> CSIOperator<T>::~CSIOperator() {
	// TODO Auto-generated destructor stub
}


template<class T> void CSIOperator<T>::mult_MH(cuNDArray<complext<T>> *in , cuNDArray<complext<T>> * out, bool accumulate){
	cuNDArray<complext<T>> tmp(*in);
	std::vector<size_t> flat_dim = {tmp.get_number_of_elements()};
	cuNDArray<complext<T>> flattened(flat_dim,tmp.get_data_ptr());
	cuNDFFT<T>::instance()->fft(&flattened,0u);
	senseOp->mult_MH(&tmp,out,accumulate);

}

template<class T> void CSIOperator<T>::mult_M(cuNDArray<complext<T>> *in , cuNDArray<complext<T>> * out, bool accumulate){
	cuNDArray<complext<T>>  * tmp = out;
	if (accumulate)
		tmp = new cuNDArray<complext<T>>(out->get_dimensions());
	size_t frequencies = *in->get_dimensions()->end();

	senseOp->mult_M(in, tmp,false);
	std::vector<size_t> flat_dim = {tmp->get_number_of_elements()};
	cuNDArray<complext<T>> flattened(flat_dim,tmp->get_data_ptr());
	cuNDFFT<T>::instance()->ifft(&flattened,0u);
	if (accumulate){
		*out += *tmp;
		delete tmp;
	}
}


template class CSIOperator<float>;

} /* namespace Gadgetron */
