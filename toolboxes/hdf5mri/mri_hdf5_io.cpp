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

#include "ismrmrd.h"

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
	//hid_t array_tid = H5Tarray_create(H5T_NATIVE_FLOAT, 1, &dims[0]);
	boost::shared_ptr<DataType> ret(new ArrayType(PredType::NATIVE_FLOAT, 1, &dims[0]));
	return ret;
}

template <> boost::shared_ptr<DataType> getHDF5ArrayType<ACE_UINT16>(int LENGTH)
{
	std::vector<hsize_t> dims(1,LENGTH);
	//hid_t array_tid = H5Tarray_create(H5T_NATIVE_USHORT, 1, &dims[0]);
	//boost::shared_ptr<DataType> ret(new DataType(array_tid));
	boost::shared_ptr<DataType> ret(new ArrayType(PredType::NATIVE_USHORT, 1, &dims[0]));
	return ret;
}



template <> boost::shared_ptr<CompType> getHDF5CompositeType<ISMRMRD::ImageHeader>()
{

	boost::shared_ptr<CompType> ret;

	try {
		ret = boost::shared_ptr<CompType>(new CompType(sizeof(ISMRMRD::ImageHeader)));

		boost::shared_ptr<DataType> matrix_size_type = getHDF5ArrayType<ACE_UINT16>(3);
		boost::shared_ptr<DataType> position_type = getHDF5ArrayType<float>(3);
		boost::shared_ptr<DataType> quaterion_type = getHDF5ArrayType<float>(4);

		ret->insertMember( "flags",              HOFFSET(ISMRMRD::ImageHeader,flags),                PredType::NATIVE_UINT);
		ret->insertMember( "matrix_size",        HOFFSET(ISMRMRD::ImageHeader, matrix_size),         *matrix_size_type);
		ret->insertMember( "channels",           HOFFSET(ISMRMRD::ImageHeader, channels),            PredType::NATIVE_USHORT);
		ret->insertMember( "position",           HOFFSET(ISMRMRD::ImageHeader, position),            *position_type);
		ret->insertMember( "quaternion",         HOFFSET(ISMRMRD::ImageHeader, quaternion),          *quaterion_type);
		ret->insertMember( "table_position",     HOFFSET(ISMRMRD::ImageHeader, patient_table_position),      PredType::NATIVE_FLOAT);
		ret->insertMember( "slice",              HOFFSET(ISMRMRD::ImageHeader, slice),               PredType::NATIVE_USHORT);
		ret->insertMember( "contrast",           HOFFSET(ISMRMRD::ImageHeader, contrast),            PredType::NATIVE_USHORT);
		ret->insertMember( "set",                HOFFSET(ISMRMRD::ImageHeader, set),                 PredType::NATIVE_USHORT);
		ret->insertMember( "phase",              HOFFSET(ISMRMRD::ImageHeader, phase),               PredType::NATIVE_USHORT);
		ret->insertMember( "average",            HOFFSET(ISMRMRD::ImageHeader, average),             PredType::NATIVE_USHORT);
		ret->insertMember( "repetition",         HOFFSET(ISMRMRD::ImageHeader, repetition),            PredType::NATIVE_USHORT);
		ret->insertMember( "time_stamp",         HOFFSET(ISMRMRD::ImageHeader, acquisition_time_stamp),          PredType::NATIVE_UINT);
		ret->insertMember( "pmu_time_stamp",     HOFFSET(ISMRMRD::ImageHeader, physiology_time_stamp),      PredType::NATIVE_UINT);
		ret->insertMember( "image_format",       HOFFSET(ISMRMRD::ImageHeader, image_data_type),        PredType::NATIVE_USHORT);
		ret->insertMember( "image_type",         HOFFSET(ISMRMRD::ImageHeader, image_type),          PredType::NATIVE_USHORT);
		ret->insertMember( "image_index",        HOFFSET(ISMRMRD::ImageHeader, image_index),         PredType::NATIVE_USHORT);
		ret->insertMember( "image_series_index", HOFFSET(ISMRMRD::ImageHeader, image_series_index),  PredType::NATIVE_USHORT);
	} catch ( ... ) {
		std::cout << "Exception caught while creating HDF5 compound datatype for ISMRMRD::ImageHeader" << std::endl;
	}


	return ret;
}


template <class T> EXPORTHDF5MRI int hdf5_append_struct(T* s,
		const char* filename, const char* varname)

{
	boost::shared_ptr<DataType> datatype = getHDF5CompositeType<T>();
	return hdf5_append_struct(s, datatype, filename, varname);
}

template EXPORTHDF5MRI int hdf5_append_struct(ISMRMRD::ImageHeader* s, const char* filename, const char* varname);

template <class T> struct local_hdf5_append_struct
{
	T h;
	hvl_t d;
};

template <class STRUCT, class DATATYPE> EXPORTHDF5MRI int hdf5_append_struct_with_data(STRUCT* s,
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

template EXPORTHDF5MRI int hdf5_append_struct_with_data(ISMRMRD::ImageHeader* s, hoNDArray<unsigned short>* a,
		                                  const char* filename, const char* varname);

template EXPORTHDF5MRI int hdf5_append_struct_with_data(ISMRMRD::ImageHeader* s, hoNDArray<float>* a,
		                                  const char* filename, const char* varname);

template EXPORTHDF5MRI int hdf5_append_struct_with_data(ISMRMRD::ImageHeader* s, hoNDArray< std::complex<float> >* a,
		                                  const char* filename, const char* varname);

template <class T> EXPORTHDF5MRI boost::shared_ptr<T> hdf5_read_struct(const char* filename, const char* varname,
		unsigned int index )
{
	boost::shared_ptr<DataType> structdatatype = getHDF5CompositeType<T>();
	boost::shared_ptr< T > ret = hdf5_read_struct<T>(structdatatype, filename,varname, index);

	return ret;
}


template <class STRUCT, class DATATYPE> EXPORTHDF5MRI header_data_struct<STRUCT, DATATYPE>
	hdf5_read_struct_with_data(const char* filename, const char* varname, unsigned index )
{

	boost::shared_ptr<DataType> structdatatype = getHDF5CompositeType<STRUCT>();
	boost::shared_ptr<DataType> vdatatype = getHDF5Type<DATATYPE>();
	vdatatype = boost::shared_ptr<DataType>(new DataType(H5Tvlen_create (vdatatype->getId())));

	CompType* ct = new CompType(sizeof(local_hdf5_append_struct<STRUCT>));
	ct->insertMember( "h", HOFFSET(local_hdf5_append_struct<STRUCT>,h),  *structdatatype);
	ct->insertMember( "d", HOFFSET(local_hdf5_append_struct<STRUCT>,d),  *vdatatype);

	boost::shared_ptr<DataType> datatype(ct);

	boost::shared_ptr< hoNDArray<  local_hdf5_append_struct<STRUCT> > > tmp = hdf5_read_array_slice<  local_hdf5_append_struct<STRUCT> >(datatype, filename, varname, index);

	header_data_struct<STRUCT, DATATYPE > ret;

	if (tmp->get_number_of_elements() != 1) {
		std::cout << "Error reading struct from HDF5 file. Expexting 1 and only 1 return value." << std::endl;
		return ret;
	}

	ret.h = boost::shared_ptr<STRUCT>(new STRUCT);
	memcpy(ret.h.get(), &tmp->get_data_ptr()[0].h, sizeof(STRUCT));

	ret.d = boost::shared_ptr< hoNDArray<DATATYPE> >(new hoNDArray<DATATYPE>());
	std::vector<unsigned int> ndarray_dims(1, tmp->get_data_ptr()[0].d.len);

	//We are taking control of this pointer from the variable length array here and now the NDArray will become responsible for deleting it
	if (!ret.d->create(&ndarray_dims, reinterpret_cast<DATATYPE*>(tmp->get_data_ptr()[0].d.p), true)) {
		std::cout << "Error allocating array for HDF5 read" << std::endl;
		return ret;
	}

	/*
	 * This code copies the data from the hvl_t (it was allocated by the HDF5 library) and then reclaims the memory
	 *
	 * With the statement above, this is not needed, since we hand control of the data to the hoNDArray thus avoiding a copy operation
	 *
	memcpy(ret.d->get_data_ptr(), tmp->get_data_ptr()[0].d.p, tmp->get_data_ptr()[0].d.len*sizeof(DATATYPE));

	std::vector<hsize_t> dims(1,1);
	DataSpace space(1, &dims[0]);
	H5Dvlen_reclaim (datatype->getId(), space.getId(), H5P_DEFAULT, &tmp->get_data_ptr()[0]);
	*/

	return ret;
}

template EXPORTHDF5MRI boost::shared_ptr<ISMRMRD::ImageHeader> hdf5_read_struct<ISMRMRD::ImageHeader>(const char* , const char* , unsigned int);
