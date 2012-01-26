/*
 * hoNDArray_hdf5_io.cpp
 *
 *  Created on: Jan 20, 2012
 *      Author: Michael S. Hansen
 */

#include "hoNDArray_hdf5_io.h"
#include <complex>
#include <FileInfo.h>
#include <H5Cpp.h>
#include <boost/shared_ptr.hpp>
#include <boost/algorithm/string.hpp>

#ifndef H5_NO_NAMESPACE
	using namespace H5;
#endif

boost::shared_ptr<H5File> OpenHF5File(const char* filename)
{
	boost::shared_ptr<H5File> ret;

	if (FileInfo(std::string(filename)).exists()) {
		ret = boost::shared_ptr<H5::H5File>(new H5File(filename, H5F_ACC_RDWR));
	} else {
		ret = boost::shared_ptr<H5::H5File>(new H5File(filename, H5F_ACC_TRUNC));
	}
	return ret;
}

bool HDF5LinkExists(H5File* f, const char* name)
{
	std::vector<std::string> name_elements;
	std::string splitstr("/");
	std::string namestr(name);
	boost::split(name_elements, namestr, boost::is_any_of(splitstr));
	std::string current_path("");
	for (unsigned int i = 0; i < name_elements.size(); i++) {
		if (name_elements[i].size() > 0) {
			current_path = current_path + std::string("/") + name_elements[i];
			if (!H5Lexists(f->getId(), current_path.c_str(), H5P_DEFAULT )) {
				return false;
			}
		}
	}
	return true;
}

int HDF5CreateGroupForDataset(H5File* f, const char* name)
{
	std::vector<std::string> name_elements;
	std::string splitstr("/");
	std::string namestr(name);
	boost::split(name_elements, namestr, boost::is_any_of(splitstr));
	std::string current_path("");
	for (unsigned int i = 0; i < name_elements.size()-1; i++) { //Skip the last one, we are just creating the group
		if (name_elements[i].size() > 0) {
			current_path = current_path + std::string("/") + name_elements[i];
			if (!H5Lexists(f->getId(), current_path.c_str(), H5P_DEFAULT )) {
				f->createGroup( current_path.c_str());
			}
		}
	}
	return 0;
}

template <class T> boost::shared_ptr<DataType> getHDF5Type();

template <> boost::shared_ptr<DataType> getHDF5Type<float>()
{
	boost::shared_ptr<DataType> ret(new DataType(H5Tcopy(H5T_NATIVE_FLOAT)));
	return ret;
}

template <> boost::shared_ptr<DataType> getHDF5Type< std::complex<float> >()
{
	boost::shared_ptr<DataType> ret(new DataType(H5Tcopy(H5T_NATIVE_FLOAT)));
	return ret;
}

template <> boost::shared_ptr<DataType> getHDF5Type< unsigned short >()
{
	boost::shared_ptr<DataType> ret(new DataType(H5Tcopy(H5T_NATIVE_USHORT)));
	return ret;
}


template <class T> unsigned short getHDF5Components();

template <> unsigned short getHDF5Components<float>() {
	return 1;
}

template <> unsigned short getHDF5Components<unsigned short>() {
	return 1;
}

template <> unsigned short getHDF5Components< std::complex<float> >() {
	return 2;
}

template <class T> int hoNDArray_hdf5_append(hoNDArray<T>* a,
		const char* filename, const char* varname)
{

	boost::shared_ptr<H5File> f = OpenHF5File(filename);
	boost::shared_ptr<DataSet> dataset;
	boost::shared_ptr<DataType> datatype = getHDF5Type<T>();

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
		if (getHDF5Components<T>() > 1) {
			dims.push_back(getHDF5Components<T>());
			max_dims.push_back(getHDF5Components<T>());
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

			if (getHDF5Components<T>() > 1) {
				adims->push_back(getHDF5Components<T>());
			}

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

template int hoNDArray_hdf5_append(hoNDArray<float>* a, const char* filename, const char* varname);
template int hoNDArray_hdf5_append(hoNDArray< std::complex<float> >* a, const char* filename, const char* varname);
template int hoNDArray_hdf5_append(hoNDArray< unsigned short >* a, const char* filename, const char* varname);
