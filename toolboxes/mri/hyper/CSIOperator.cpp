/*
 * CSIOperator.cpp
 *
 *  Created on: Nov 10, 2014
 *      Author: dch
 */

#include "CSIOperator.h"
#include "cuNDFFT.h"
#include "cuNDArray_math.h"
#include "CSI_utils.h"

namespace Gadgetron {


template<class T> CSIOperator<T>::CSIOperator() {
	// TODO Auto-generated constructor stub
}

template<class T> CSIOperator<T>::CSIOperator(T dtt, T dte ) : dtt_(dtt), dte_(dte) {
	// TODO Auto-generated constructor stub

}

template<class T> CSIOperator<T>::~CSIOperator() {
	// TODO Auto-generated destructor stub
}


template<class T> void CSIOperator<T>::mult_MH(cuNDArray<complext<T>> *in , cuNDArray<complext<T>> * out, bool accumulate){
	std::vector<size_t> kdim = *in->get_dimensions();
	kdim[2] = frequencies.get_number_of_elements();
	cuNDArray<complext<T>> tmp(kdim);
	//cuNDFFT<float>::instance()->fft(&tmp,0u); //FFT along the TE dimension
	CSI_dft(&tmp,in,&frequencies,dte_,dtt_);
	senseOp->mult_MH(&tmp,out,accumulate);

}

template<class T> void CSIOperator<T>::mult_M(cuNDArray<complext<T>> *in , cuNDArray<complext<T>> * out, bool accumulate){
	auto out_tmp = out;
	if (accumulate) out_tmp = new cuNDArray<complext<T>>(out);
	std::vector<size_t> kdim = *out->get_dimensions();
	kdim[2] = frequencies.get_number_of_elements();
	cuNDArray<complext<T>> tmp(kdim);
	senseOp->mult_M(in, &tmp,false);
	CSI_dftH(&tmp,out_tmp,&frequencies,dte_,dtt_);
	if (accumulate){
		*out = *out_tmp;
		delete out_tmp;
	}
}


template class CSIOperator<float>;

} /* namespace Gadgetron */
