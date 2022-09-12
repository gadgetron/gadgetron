/*
 * CSfreqOperator.h
 *
 *  Created on: Dec 2, 2014
 *      Author: dch
 */

#ifndef CSFREQOPERATOR_H_
#define CSFREQOPERATOR_H_

#include "linearOperator.h"
#include "cuNDArray.h"
#include "CSI_utils.h"
namespace Gadgetron{
class CSfreqOperator : public linearOperator<cuNDArray<float_complext> > {

public:

	CSfreqOperator(){};
	CSfreqOperator(float dtt_, float dte_) : dtt(dtt_), dte(dte_){

}

	virtual  void mult_M(cuNDArray<float_complext> * in, cuNDArray<float_complext>* out,bool accumulate){
		auto tmp_out  = out;
		if (accumulate) tmp_out = new cuNDArray<float_complext>(*out);
		CSI_dftH(in,tmp_out,&freqs,dte,dtt);
		if (accumulate){
			*out += *tmp_out;
			delete tmp_out;
		}
	}
	virtual  void mult_MH(cuNDArray<float_complext> * in, cuNDArray<float_complext>* out,bool accumulate){
		auto tmp_out  = out;
		if (accumulate) tmp_out = new cuNDArray<float_complext>(*out);
		CSI_dft(tmp_out,in,&freqs,dte,dtt);
		if (accumulate){
			*out += *tmp_out;
			delete tmp_out;
		}
	}
        void set_frequencies(std::vector<float>& freq) {
            freqs = cuNDArray<float>(freq.size());
            cudaMemcpy(freqs.data(), freq.data(), freqs.get_number_of_bytes(), cudaMemcpyKind::cudaMemcpyHostToDevice);
		}

	cuNDArray<float> freqs;
	float dtt,dte;
};
}


#endif /* CSFREQOPERATOR_H_ */
