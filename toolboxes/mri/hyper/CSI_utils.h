/*
 * CSI_utils.h
 *
 *  Created on: Nov 20, 2014
 *      Author: dch
 */

#pragma once

#include "cuNDArray.h"
#include <thrust/device_vector.h>
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
	template<class T> void CSI_dft(cuNDArray<complext<T> >* kspace, cuNDArray<complext<T> >* tspace, thrust::device_vector<T>* frequencies, T dtt, T dte);
	/**
	 * Performs the adjoint of the non-cartesian discrete fourier transform.
	 * @param kspace The input kspace
	 * @param tspace The output time space
	 * @param frequencies Frequencies on which to do DFT.
	 * @param dte Time step between points in the first dimension of the tspace
	 * @param dtt Time step between points in the second dimension of the tspace
	 */
	template<class T> void CSI_dftH(cuNDArray<complext<T> >* kspace, cuNDArray<complext<T> >* tspace, thrust::device_vector<T>* frequencies, T dte, T dtt);

}
