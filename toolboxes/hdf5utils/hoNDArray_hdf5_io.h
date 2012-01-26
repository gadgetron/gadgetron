/*
 * hoNDArray_hdf5_io.h
 *
 *  Created on: Jan 20, 2012
 *      Author: Michael S. Hansen
 */

#ifndef HONDARRAY_HDF5_IO_H_
#define HONDARRAY_HDF5_IO_H_

#include "hdf5utils_export.h"

#include <hoNDArray.h>

template <class T> EXPORTHDF5UTILS int hoNDArray_hdf5_append(hoNDArray<T>* a,
		const char* filename, const char* varname);


#endif /* HONDARRAY_HDF5_IO_H_ */
