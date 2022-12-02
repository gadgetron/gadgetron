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
	kdim[1] =frequencies.size(); //Set to number of time samples rather than number of frequencies
	cuNDArray<complext<T>> tmp(kdim);
	//senseOp->mult_MH(in,out,accumulate);

	CSI_dft(&tmp,in,&frequencies,dtt_,dte_);
	senseOp->mult_MH(&tmp,out,accumulate);
	//cuNDFFT<float>::instance()->fft(&tmp,0u); //FFT along the TE dimension
}

template <class T> void CSIOperator<T>::set_frequencies(std::vector<T>& freq) {
    frequencies = cuNDArray<T>(freq.size());
    cudaMemcpy(frequencies.data(), freq.data(), frequencies.get_number_of_bytes(), cudaMemcpyKind::cudaMemcpyHostToDevice);
    
}

template<class T> void CSIOperator<T>::mult_M(cuNDArray<complext<T>> *in , cuNDArray<complext<T>> * out, bool accumulate){
	cuNDArray<complext<T>>* out_tmp = out;
	if (accumulate) out_tmp = new cuNDArray<complext<T>>(out->get_dimensions());
	std::vector<size_t> kdim = *out->get_dimensions();
	kdim[1] =frequencies.size();
	cuNDArray<complext<T>> tmp(kdim);
	senseOp->mult_M(in,&tmp,accumulate);

	CSI_dftH(&tmp,out_tmp,&frequencies,dtt_,dte_);
	if (accumulate){
		*out += *out_tmp;
		delete out_tmp;
	}
}


template class EXPORTHYPER CSIOperator<float>;

} /* namespace Gadgetron */
