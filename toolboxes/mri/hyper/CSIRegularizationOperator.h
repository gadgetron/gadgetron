/*
 * CSIRegularizationOperator.h
 *
 *  Created on: Feb 17, 2015
 *      Author: dch
 */

#pragma once
#include "linearOperator.h"
#include "sense_utilities.h"
namespace Gadgetron {

template<class T> class CSIRegularizationOperator: public Gadgetron::linearOperator {
public:
	CSIRegularizationOperator(T dtt, T dte, boost::shared_ptr<cuNDArray<complext<T>>> freq_cal, std::vector<T> & freqs, boost::shared_ptr<cuNDArray<complext<T>> csm): dtt_(dtt), dte_(dte), csm_(csm){
		frequencies = thrust::device_vector<T>(freqs.begin(),freqs.end());
		std::vector<size_t> dims {1,freqs.size(),freq_cal->get_size(2)};
		freq_calibration = boost::make_shared<cuNDArray<complext<T>>>(dims);

		CSI_dft(freq_cal.get(),freq_calibration,&frequencies,dtt_,dte_);

		ncoils = csm_->get_dimensions()->back();

	};
	virtual ~CSIRegularizationOperator();

	virtual void mult_M(cuNDArray<complext<T>>* in, cuNDArray<complext<T>>* out, bool accumulate) override {
		auto in_dims = *in->get_dimension();
		in_dims->push_back(ncoils);
		cuNDArray<complext<T>> coiled(in_dims);
		csm_mult_M(in,&coiled,csm_.get());
		std::vector<size_t> new_dims{coiled.get_number_of_elements()/(ncoils*frequencies.size()),frequencies.size(),ncoils};
		coiled.reshape(new_dims);
		auto meaned = sum(&coiled,0);
		*meaned /= T(new_dims[0]);

	}

	T dtt_, dte_;
	boost::shared_ptr<cuNDArray<complext<T>>> freq_calibration;
	boost::shared_ptr<cuNDArray<complext<T>>> csm_;
	thrust::device_vector<T> frequencies;
	size_t ncoils;
};

} /* namespace Gadgetron */
