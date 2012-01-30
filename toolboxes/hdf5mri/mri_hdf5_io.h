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


template <class T> EXPORTHDF5MRI boost::shared_ptr<T> hdf5_read_struct(const char* filename, const char* varname,
		unsigned int index = 0);

template <class STRUCT, class DATATYPE> struct header_data_struct
{
	boost::shared_ptr<STRUCT> h;
	boost::shared_ptr< hoNDArray<DATATYPE> > d;
};

template <class STRUCT, class DATATYPE> EXPORTHDF5MRI header_data_struct<STRUCT, DATATYPE>
	hdf5_read_struct_with_data(const char* filename, const char* varname, unsigned index = 0);

#endif /* MRI_HDF5_IO_H_ */
