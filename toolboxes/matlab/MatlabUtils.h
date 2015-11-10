/*
 * MatlabUtils.h
 *
 *  Created on: Dec 5, 2014
 *      Author: dch
 */

#ifndef MATLABUTILS_H_
#define MATLABUTILS_H_

#include "engine.h"
#include "hoNDArray.h"
#include "mri_core_data.h"
namespace Gadgetron{

/**
 * Creates a Matlab array from an hoNDArray
 * @param
 * @return
 */
template<class T> mxArray* hoNDArrayToMatlab(hoNDArray<T>* );

/**
 * Create hoNDArray from a Matlab array. Will attempt type conversion.
 * @param
 * @return
 */
template<class T> hoNDArray<T>* MatlabToHoNDArray(mxArray*);

/**
 * Creates a matlab struct from an IsmrmrdDataBuffer
 * @param buffer
 * @return
 */
mxArray* BufferToMatlabStruct(IsmrmrdDataBuffered* buffer);


IsmrmrdDataBuffered MatlabStructToBuffer(mxArray* mxstruct);

/**
 * Create a Matlab struct from a SamplingDescription
 * @param samp
 * @return
 */
mxArray* samplingdescriptionToMatlabStruct(SamplingDescription* samp);
}
#endif /* MATLABUTILS_H_ */
