/*
 * mri_hdf5_io.cpp
 *
 *  Created on: Jan 26, 2012
 *      Author: Michael S. Hansen
 */

#include "mri_hdf5_io.h"

#include <H5Cpp.h>
#include "hdf5_core.h"
#include "hoNDArray_hdf5_io.h"

#include "GadgetMRIHeaders.h"

#include <vector>
#include <complex>

#ifndef H5_NO_NAMESPACE
	using namespace H5;
#endif

template <class T> boost::shared_ptr<CompType> getHDF5CompositeType();
template <class T> boost::shared_ptr<DataType> getHDF5ArrayType(int LENGTH);

template <> boost::shared_ptr<DataType> getHDF5ArrayType<float>(int LENGTH)
{
	std::vector<hsize_t> dims(1,LENGTH);
	hid_t array_tid = H5Tarray_create(H5T_NATIVE_FLOAT, 1, &dims[0]);
	boost::shared_ptr<DataType> ret(new DataType(array_tid));
	return ret;
}

template <> boost::shared_ptr<DataType> getHDF5ArrayType<ACE_UINT16>(int LENGTH)
{
	std::vector<hsize_t> dims(1,LENGTH);
	hid_t array_tid = H5Tarray_create(H5T_NATIVE_USHORT, 1, &dims[0]);
	boost::shared_ptr<DataType> ret(new DataType(array_tid));
	return ret;
}

template <> boost::shared_ptr<CompType> getHDF5CompositeType<LoopCounters>()
{
	boost::shared_ptr<CompType> ret;

	try {
		ret = boost::shared_ptr<CompType>(new CompType(sizeof(LoopCounters)));

		ret->insertMember( "line",        HOFFSET(LoopCounters,line),        PredType::NATIVE_USHORT);
		ret->insertMember( "acquisition", HOFFSET(LoopCounters,acquisition), PredType::NATIVE_USHORT);
		ret->insertMember( "slice",       HOFFSET(LoopCounters,slice),       PredType::NATIVE_USHORT);
		ret->insertMember( "partition",   HOFFSET(LoopCounters,partition),   PredType::NATIVE_USHORT);
		ret->insertMember( "echo",        HOFFSET(LoopCounters,echo),        PredType::NATIVE_USHORT);
		ret->insertMember( "phase",       HOFFSET(LoopCounters,phase),       PredType::NATIVE_USHORT);
		ret->insertMember( "repetition",  HOFFSET(LoopCounters,repetition),  PredType::NATIVE_USHORT);
		ret->insertMember( "set",         HOFFSET(LoopCounters,set),         PredType::NATIVE_USHORT);
		ret->insertMember( "segment",     HOFFSET(LoopCounters,segment),     PredType::NATIVE_USHORT);
		ret->insertMember( "channel",     HOFFSET(LoopCounters,channel),     PredType::NATIVE_USHORT);
	} catch ( ... ) {
		std::cout << "Exception caught while creating HDF5 compound datatype for LoopCounter" << std::endl;
	}
	return ret;
}

template <> boost::shared_ptr<CompType> getHDF5CompositeType<GadgetMessageAcquisition>()
{
	boost::shared_ptr<CompType> ret(new CompType(sizeof(GadgetMessageAcquisition)));
	ret->insertMember( "flags",              HOFFSET(GadgetMessageAcquisition,flags),               PredType::NATIVE_UINT);
	ret->insertMember( "meas_uid",           HOFFSET(GadgetMessageAcquisition,meas_uid),            PredType::NATIVE_UINT);
	ret->insertMember( "scan_counter",       HOFFSET(GadgetMessageAcquisition, scan_counter),       PredType::NATIVE_UINT);
	ret->insertMember( "time_stamp",         HOFFSET(GadgetMessageAcquisition, time_stamp),         PredType::NATIVE_UINT);
	ret->insertMember( "samples",            HOFFSET(GadgetMessageAcquisition, samples),            PredType::NATIVE_USHORT);
	ret->insertMember( "channels",           HOFFSET(GadgetMessageAcquisition, channels),           PredType::NATIVE_USHORT);

	boost::shared_ptr<DataType> position_type = getHDF5ArrayType<float>(3);
	boost::shared_ptr<DataType> quarterion_type = getHDF5ArrayType<float>(4);
	boost::shared_ptr<CompType> loopcounters_type = getHDF5CompositeType<LoopCounters>();

	ret->insertMember( "position",           HOFFSET(GadgetMessageAcquisition, position),       *position_type);
	ret->insertMember( "quarternion",        HOFFSET(GadgetMessageAcquisition, quarternion),        *quarterion_type);
	ret->insertMember( "table_position",     HOFFSET(GadgetMessageAcquisition, table_position),     PredType::NATIVE_FLOAT);
	ret->insertMember( "idx",                HOFFSET(GadgetMessageAcquisition, idx),                *loopcounters_type);
	ret->insertMember( "min_idx",            HOFFSET(GadgetMessageAcquisition, min_idx),            *loopcounters_type);
	ret->insertMember( "max_idx",            HOFFSET(GadgetMessageAcquisition, max_idx),            *loopcounters_type);

	return ret;
}


template <> boost::shared_ptr<CompType> getHDF5CompositeType<GadgetMessageImage>()
{

	boost::shared_ptr<CompType> ret;

	try {
		ret = boost::shared_ptr<CompType>(new CompType(sizeof(GadgetMessageImage)));

		boost::shared_ptr<DataType> matrix_size_type = getHDF5ArrayType<ACE_UINT16>(3);
		boost::shared_ptr<CompType> loopcounters_type = getHDF5CompositeType<LoopCounters>();
		boost::shared_ptr<DataType> position_type = getHDF5ArrayType<float>(3);
		boost::shared_ptr<DataType> quarterion_type = getHDF5ArrayType<float>(4);

		ret->insertMember( "flags",              HOFFSET(GadgetMessageImage,flags),                PredType::NATIVE_UINT);
		ret->insertMember( "matrix_size",        HOFFSET(GadgetMessageImage, matrix_size),         *matrix_size_type);
		ret->insertMember( "channels",           HOFFSET(GadgetMessageImage, channels),            PredType::NATIVE_USHORT);
		ret->insertMember( "position",           HOFFSET(GadgetMessageImage, position),            *position_type);
		ret->insertMember( "quarternion",        HOFFSET(GadgetMessageImage, quarternion),         *quarterion_type);
		ret->insertMember( "table_position",     HOFFSET(GadgetMessageImage, table_position),      PredType::NATIVE_FLOAT);
		ret->insertMember( "data_idx_min",       HOFFSET(GadgetMessageImage, data_idx_min),        *loopcounters_type);
		ret->insertMember( "data_idx_max",       HOFFSET(GadgetMessageImage, data_idx_max),        *loopcounters_type);
		ret->insertMember( "data_idx_current",   HOFFSET(GadgetMessageImage, data_idx_current),    *loopcounters_type);
		ret->insertMember( "time_stamp",         HOFFSET(GadgetMessageImage, time_stamp),          PredType::NATIVE_UINT);
		ret->insertMember( "image_format",       HOFFSET(GadgetMessageImage, image_format),        PredType::NATIVE_USHORT);
		ret->insertMember( "image_type",         HOFFSET(GadgetMessageImage, image_type),          PredType::NATIVE_USHORT);
		ret->insertMember( "image_index",        HOFFSET(GadgetMessageImage, image_index),         PredType::NATIVE_USHORT);
		ret->insertMember( "image_series_index", HOFFSET(GadgetMessageImage, image_series_index),  PredType::NATIVE_USHORT);
	} catch ( ... ) {
		std::cout << "Exception caught while creating HDF5 compound datatype for GadgetMessageImage" << std::endl;
	}


	return ret;
}


template <class T> int hdf5_append_struct(T* s,
		boost::shared_ptr<DataType> datatype,
		const char* filename, const char* varname)

{
	hoNDArray<T> tmp;
	std::vector<unsigned int> dims(1,1);

	tmp.create(&dims, s, false); //This is just a dummy container for the data

	return hdf5_append_array(&tmp, datatype, filename, varname);

	return 0;
}

template <class T> int hdf5_append_struct(T* s,
		const char* filename, const char* varname)

{
	boost::shared_ptr<DataType> datatype = getHDF5CompositeType<T>();
	return hdf5_append_struct(s, datatype, filename, varname);
}

template int hdf5_append_struct(GadgetMessageAcquisition* s, const char* filename, const char* varname);
template int hdf5_append_struct(GadgetMessageImage* s, const char* filename, const char* varname);

template <class T> struct local_hdf5_append_struct
{
	T h;
	hvl_t d;
};

template <class STRUCT, class DATATYPE> int hdf5_append_struct_with_data(STRUCT* s,
		hoNDArray<DATATYPE>* a, const char* filename, const char* varname)
{

	local_hdf5_append_struct<STRUCT> tmp;
    tmp.h = *s;

    tmp.d.len = a->get_number_of_elements();
    tmp.d.p = (void*) a->get_data_ptr();

	boost::shared_ptr<DataType> structdatatype = getHDF5CompositeType<STRUCT>();
	boost::shared_ptr<DataType> vdatatype = getHDF5Type<DATATYPE>();
	vdatatype = boost::shared_ptr<DataType>(new DataType(H5Tvlen_create (vdatatype->getId())));

	CompType* ct = new CompType(sizeof(local_hdf5_append_struct<STRUCT>));
	ct->insertMember( "h", HOFFSET(local_hdf5_append_struct<STRUCT>,h),  *structdatatype);
	ct->insertMember( "d", HOFFSET(local_hdf5_append_struct<STRUCT>,d),  *vdatatype);

	boost::shared_ptr<DataType> datatype(ct);


	return hdf5_append_struct(&tmp, datatype, filename, varname);
}

template int hdf5_append_struct_with_data(GadgetMessageImage* s, hoNDArray<unsigned short>* a,
		                                  const char* filename, const char* varname);

template int hdf5_append_struct_with_data(GadgetMessageImage* s, hoNDArray<float>* a,
		                                  const char* filename, const char* varname);

template int hdf5_append_struct_with_data(GadgetMessageImage* s, hoNDArray< std::complex<float> >* a,
		                                  const char* filename, const char* varname);

template int hdf5_append_struct_with_data(GadgetMessageAcquisition* s, hoNDArray< std::complex<float> >* a,
		                                  const char* filename, const char* varname);
