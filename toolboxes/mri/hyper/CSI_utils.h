/*
 * CSI_utils.h
 *
 *  Created on: Nov 20, 2014
 *      Author: dch
 */

#pragma once

#include "cuNDArray.h"
#include <thrust/device_vector.h>
#include "gadgetron_toolbox_hyper_export.h"
namespace Gadgetron {

	/**
	 * Performs a non-cartesian discrete fourier transform (DFT) along the specified frequencies, at the specified time intervals.
	 * Note that this should be no used for general purpose DFTS as it will be mindnumbingly slow.
	 * @param kspace The output kspace
	 * @param tspace The input time space
	 * @param frequencies The frequencies on which to do DFT.
	 * @param dtt Time step between points in the first dimension of the tspace
	 * @param dte Time step between points in the second dimension of the tspace
	 */
	template<class T> EXPORTHYPER void CSI_dft(cuNDArray<complext<T> >* kspace, cuNDArray<complext<T> >* tspace, cuNDArray<T>* frequencies, T dtt, T dte);
	/**
	 * Performs the adjoint of the non-cartesian discrete fourier transform.
	 * @param kspace The input kspace
	 * @param tspace The output time space
	 * @param frequencies Frequencies on which to do DFT.
	 * @param dte Time step between points in the first dimension of the tspace
	 * @param dtt Time step between points in the second dimension of the tspace
	 */
	template<class T> EXPORTHYPER void CSI_dftH(cuNDArray<complext<T> >* kspace, cuNDArray<complext<T> >* tspace, cuNDArray<T>* frequencies, T dte, T dtt);

	template<class T> EXPORTHYPER boost::shared_ptr<cuNDArray<complext<T> > 	> calculate_frequency_calibration(cuNDArray<complext<T> >* time_track, cuNDArray<T>* frequencies,cuNDArray<complext<T> > * csm,T dtt,T dte);

	template<class T> EXPORTHYPER void mult_freq(cuNDArray<complext<T> >* in_out, cuNDArray<complext<T> >* freqs, bool conjugate);
}
