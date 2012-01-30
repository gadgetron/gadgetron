/*
 * hoNDArray_hdf5_io.h
 *
 *  Created on: Jan 20, 2012
 *      Author: Michael S. Hansen
 */

#ifndef HONDARRAY_HDF5_IO_H_
#define HONDARRAY_HDF5_IO_H_


#include "hdf5utils_export.h"

#include <complex>
#include <FileInfo.h>
#include <H5Cpp.h>
#include <boost/shared_ptr.hpp>

#include "hdf5_core.h"

#ifndef H5_NO_NAMESPACE
	using namespace H5;
#endif


/**
 *   Base class template for returning HDF5 data type
 */
template <class T> boost::shared_ptr<DataType> getHDF5Type();

/**
 *   Returns HDF5 data type for single precision float
 */
template <> boost::shared_ptr<DataType> getHDF5Type<float>()
{
	boost::shared_ptr<DataType> ret(new DataType(H5Tcopy(H5T_NATIVE_FLOAT)));
	return ret;
}

/**
 *   Returns HDF5 data type for char
 */
template <> boost::shared_ptr<DataType> getHDF5Type<char>()
{
	boost::shared_ptr<DataType> ret(new DataType(H5Tcopy(H5T_NATIVE_CHAR)));
	return ret;
}

/**
 *   Returns HDF5 compound data type for single precision complex float
 */
template <> boost::shared_ptr<DataType> getHDF5Type< std::complex<float> >()
{
	CompType* ct = new CompType(sizeof( std::complex<float> ));
	ct->insertMember( "real",  0,              PredType::NATIVE_FLOAT);
	ct->insertMember( "imag",  sizeof(float),  PredType::NATIVE_FLOAT);
	boost::shared_ptr<DataType> ret(ct);
	return ret;
}

/**
 *   Returns HDF5 data type for unsigned short int
 */
template <> boost::shared_ptr<DataType> getHDF5Type< unsigned short >()
{
	boost::shared_ptr<DataType> ret(new DataType(H5Tcopy(H5T_NATIVE_USHORT)));
	return ret;
}


/**
 *  This is the core function for writing data to HDF 5 files.
 *  All functions that write data to HDF5 file ultimately call this function
 *
 *  Given a @param filename and a @param varname, it will append the array to the variable in @param varname
 *  in the file in @param filename.
 *
 *  If the Array has the dimensions, [128, 128,4], the resulting HDF5 variable would have the dimensions
 *  [128,128,4,1] after appending the first array and [128,128,4,2] after appending the second array and so on.
 *
 *  Remember HDF5 uses "C-style" array dimensions, so in the HDF5 file, the variable would have dimensions [2,4,128,128].
 *
 */
template <class T> int hdf5_append_array(hoNDArray<T>* a,
		boost::shared_ptr<DataType> datatype,
		const char* filename, const char* varname)
{

	boost::shared_ptr<H5File> f = OpenHF5File(filename);
	boost::shared_ptr<DataSet> dataset;

	std::vector<hsize_t> dims;
	std::vector<hsize_t> max_dims;
	if (!HDF5LinkExists(f.get(), varname)) {
		std::vector<hsize_t> dims;
		std::vector<hsize_t> max_dims;
		dims.push_back(1);
		max_dims.push_back(H5S_UNLIMITED);
		boost::shared_ptr< std::vector<unsigned int> > adims = a->get_dimensions();
		for (int i = a->get_number_of_dimensions()-1; i >= 0; i--) {
			dims.push_back(static_cast<unsigned long long int>(a->get_size(i)));
			max_dims.push_back(static_cast<unsigned long long int>(a->get_size(i)));
		}
		try {

			if (HDF5CreateGroupForDataset(f.get(), varname) < 0) {
				std::cout << "Failed to create group in HDF 5 file." << std::endl;
				return -1;
			}

			DataSpace mspace1( dims.size(), &dims[0], &max_dims[0]);

			DSetCreatPropList cparms;
			cparms.setChunk( dims.size(), &dims[0] );

			dataset = boost::shared_ptr<DataSet>(new DataSet(f->createDataSet( varname, *datatype, mspace1, cparms)));
			mspace1 = dataset->getSpace();

			DataSpace mspace2( dims.size(), &dims[0] );

			std::vector<hsize_t> offset(dims.size());
			mspace1.selectHyperslab(H5S_SELECT_SET, &dims[0], &offset[0]);
			dataset->write( a->get_data_ptr(), *datatype, mspace2, mspace1 );

		} catch( Exception& e ) {
			std::cout << "Exception caught while creating HDF5 dataset" << std::endl;
			std::cout << e.getDetailMsg() << std::endl;
			return -1;
		}
	} else {
		try {  // to determine if the dataset exists in the group
			dataset = boost::shared_ptr<DataSet>(new DataSet(f->openDataSet(varname)));

			DataType mtype = dataset->getDataType();
			if (!(mtype == (*datatype))) {
				std::cout << "Attempting to append data to HDF5 dataset with the wrong type" << std::endl;
				return -1;
			}

			DataSpace mspace1 = dataset->getSpace();
			int rank = mspace1.getSimpleExtentNdims();
			std::vector<hsize_t> ddims(rank,0);
			mspace1.getSimpleExtentDims(&ddims[0], NULL);

			boost::shared_ptr< std::vector<unsigned int> > adims = a->get_dimensions();


			if ((ddims.size()-1) != adims->size()) {
				std::cout << "Dimensions in dataset does not match with existing HDF5 dataset" << std::endl;
			}

			dims.push_back(1);
			for (unsigned int i = 1; i < ddims.size(); i++) {
				if ((*adims)[ddims.size()-1-i] != ddims[i]) {
					std::cout << "Error trying to write array to existing HDF5 file. Variable has wrong size." << std::endl;
					std::cout << (*adims)[ddims.size()-1-i] << ", " << ddims[i] << std::endl;
					return -1;
				}
				dims.push_back(ddims[i]);
			}

			std::vector<hsize_t> offset(rank, 0);
			offset[0] = ddims[0];

			ddims[0]++;

			dataset->extend(&ddims[0]);

			DataSpace fspace2 = dataset->getSpace();
			fspace2.selectHyperslab( H5S_SELECT_SET, &dims[0], &offset[0] );

			DataSpace mspace2( rank, &dims[0] );

			dataset->write( a->get_data_ptr(), *datatype, mspace2, fspace2 );

		 }
		 catch( FileIException& not_found_error)
		 {
			 std::cout << "Dataset is not found. At this point, it should have been created!" << std::endl;
			 return -1;
		 }
	}

	return 0;
}


template <class T> boost::shared_ptr< hoNDArray<T> > hdf5_read_array_slice(
		boost::shared_ptr<DataType> datatype,
		const char* filename, const char* varname, unsigned int index = 0)
{
	boost::shared_ptr< hoNDArray<T> > ret(new hoNDArray<T>());

	if (!FileInfo(std::string(filename)).exists()) {
		std::cout << "Trying to open non-existing HDF5 file" << std::endl;
		return ret;
	}

	try {
		boost::shared_ptr<H5File> f = OpenHF5File(filename);
		if (!HDF5LinkExists(f.get(), varname)) {
			std::cout << "Trying to access non-existing variable in HDF5 file." << std::endl;
			return ret;
		}

		DataSet d = f->openDataSet(varname);

		DataSpace dataspace = d.getSpace();

		DataType dtype = d.getDataType();

		if (!(dtype == *datatype)) {
			std::cout << "HDF5 datatype for selected variable does not match signature of reading function" << std::endl;
			return ret;
		}

		int rank = dataspace.getSimpleExtentNdims();
		std::vector<hsize_t> dims(rank,0);

		dataspace.getSimpleExtentDims(&dims[0]);

		if (dims[0] <= index) {
			std::cout << "Attempting to access non-existing hyperslice" << std::endl;
			return ret;
		}

		std::vector<hsize_t> slice_dims(rank,0);
		std::vector<hsize_t> offset(rank,0);
		std::vector<unsigned int> ndarray_dims(rank-1,0);
		slice_dims[0] = 1;
		offset[0] = index;

		for (unsigned int i = 1; i < rank; i++) {
			slice_dims[i] = dims[i];
			ndarray_dims[ndarray_dims.size()-i] = dims[i]; //Flip dimensions for NDArray.
		}

		/*
		for (unsigned int i = 0; i < slice_dims.size(); i++) {
			std::cout << "slice_dims[" << i << "] = " << slice_dims[i] << std::endl;
			std::cout << "offset[" << i << "] = " << offset[i] << std::endl;

		}
		*/

		dataspace.selectHyperslab( H5S_SELECT_SET, &slice_dims[0], &offset[0] );

		DataSpace memspace(rank,&slice_dims[0]);

		//ret = boost::shared_ptr<  hoNDArray<T>() >( new hoNDArray<T>() );

		if (!ret->create(&ndarray_dims)) {
			std::cout << "Failed to create NDArray for HDF5 read" << std::endl;
			return ret;
		}

		//OK finally ready, now read the data.
		d.read(reinterpret_cast<void*>(ret->get_data_ptr()), *datatype, memspace, dataspace);

	} catch (...) {
		std::cout << "Error caught while attempting to read HDF5 file" << std::endl;
		return ret;
	}

	return ret;

}

/**
 * Wrapper function for calling hdf5_read_array_slice
 */
template <class T> boost::shared_ptr< hoNDArray<T> > hdf5_read_array_slice(
		const char* filename, const char* varname, unsigned int index = 0)
{
	boost::shared_ptr<DataType> datatype = getHDF5Type<T>();
	return hdf5_read_array_slice<T>(datatype, filename, varname, index);
}

/**
 * Wrapper function for calling hdf5_append_array
 */
template <class T> EXPORTHDF5UTILS int hdf5_append_array(hoNDArray<T>* a,
		const char* filename, const char* varname)
{
	boost::shared_ptr<DataType> datatype = getHDF5Type<T>();
	return hdf5_append_array(a, datatype, filename, varname);
}

#endif /* HONDARRAY_HDF5_IO_H_ */
