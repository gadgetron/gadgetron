/*
 * mri_hdf5_io.h
 *
 *  Created on: Jan 26, 2012
 *      Author: Michael S. Hansen
 */

#ifndef MRI_HDF5_IO_H_
#define MRI_HDF5_IO_H_

#include "hdf5mri_export.h"
#include <H5Cpp.h>
#include <hoNDArray.h>

#ifndef H5_NO_NAMESPACE
	using namespace H5;
#endif

template <class T> EXPORTHDF5MRI int hdf5_append_struct(T* s,
		const char* filename, const char* varname);

template <class STRUCT, class DATATYPE> EXPORTHDF5MRI int hdf5_append_struct_with_data(STRUCT* s,
		hoNDArray<DATATYPE>* a, const char* filename, const char* varname);

#endif /* MRI_HDF5_IO_H_ */
