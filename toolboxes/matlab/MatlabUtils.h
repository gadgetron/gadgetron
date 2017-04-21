/*
 * MatlabUtils.h
 *
 *  Created on: Dec 5, 2014
 *      Author: dch
 */

#ifndef MATLABUTILS_H_
#define MATLABUTILS_H_

#include "matlab_export.h"
#include "engine.h"
#include "hoNDArray.h"
#include "hoNDImage.h"
#include "mri_core_data.h"

namespace Gadgetron{

/**
 * Creates a Matlab array from an hoNDArray
 * @param
 * @return
 */
template<class T> EXPORTMATLAB mxArray* hoNDArrayToMatlab(hoNDArray<T>* );

/**
 * Create hoNDArray from a Matlab array. Will attempt type conversion.
 * @param
 * @return
 */
template<class T> EXPORTMATLAB hoNDArray<T> MatlabToHoNDArray(mxArray*);
template<class T> EXPORTMATLAB void MatlabToHoNDArray(mxArray* m, hoNDArray<T>& a);

/**
* Creates a Matlab array from an hoNDImage
* @param
* @return
*/
template<class T, unsigned int D> EXPORTMATLAB mxArray* hoNDImageToMatlab(const hoNDImage<T, D>* a, mxArray*& h );

/**
* Create hoNDImage from a Matlab array.
* @param
* @return
*/
template<class T, unsigned int D> EXPORTMATLAB void MatlabToHoNDImage(const mxArray* m, const mxArray* header, hoNDImage<T, D>& a);

/**
 * Creates a matlab struct from an IsmrmrdDataBuffer
 * @param buffer
 * @return
 */
EXPORTMATLAB mxArray* BufferToMatlabStruct(IsmrmrdDataBuffered* buffer, bool omitData = false);

/**
 * LA's data splitting
 * @param
 * @return
 */
EXPORTMATLAB mxArray* GetSplitReconData(IsmrmrdDataBuffered* buffer, size_t index_begin, size_t index_end);

EXPORTMATLAB IsmrmrdDataBuffered MatlabStructToBuffer(mxArray* mxstruct);

/**
 * Create a Matlab struct from a SamplingDescription
 * @param samp
 * @return
 */
EXPORTMATLAB mxArray* samplingdescriptionToMatlabStruct(SamplingDescription* samp);

/**
* Create a Matlab array from std::vector
* @param 
* @return
*/
template<class T> EXPORTMATLAB mxArray* StdVecToMatlab(const std::vector<T>*);

/**
* Create a std vector from a Matlab array
* @param 
* @return
*/
template<class T> EXPORTMATLAB void MatlabToStdVec(const mxArray*, std::vector<T>& a);

/**
* Create a Matlab array from std::string
* @param 
* @return
*/
EXPORTMATLAB mxArray* StdStringToMatlab(const std::string*);

/**
* Create a string from a Matlab array
* @param 
* @return
*/
EXPORTMATLAB std::string MatlabToStdString(const mxArray*);

}
#endif /* MATLABUTILS_H_ */
