/*
 * hdf5_test.cpp
 *
 *  Created on: Jan 20, 2012
 *      Author: Michael S. Hansen
 */

#include <hoNDArray.h>
#include <hoNDArray_fileio.h>
#include "hoNDArray_hdf5_io.h"

#include <iostream>
#include <complex>

int main(int argc, char** argv) {

	if (argc < 4) {
		std::cout << "Usage: " << std::endl;
		std::cout << " - hdf_test <arrayfilename> <hdf5file> <varname> <type>" << std::endl;
		std::cout << " - types: float (0), cplx float (1), unsigned short (2)" << std::endl;
		return -1;
	}

	int datatype = 0;
	if (argc >= 5) {
		datatype = atoi(argv[4]);
	}


	std::string arrayfilename = std::string(argv[1]);
	std::string hdf5filename = std::string(argv[2]);
	std::string varname = std::string(argv[3]);

	std::cout << "Appending " << arrayfilename
			<< " to HDF5 file " << hdf5filename
			<< " (" << varname << ")" << std::endl;

	if (datatype == 0) {
		boost::shared_ptr< hoNDArray< float > > A = read_nd_array< float >(arrayfilename.c_str());
		if (hoNDArray_hdf5_append(A.get(), hdf5filename.c_str(), varname.c_str()) < 0) {
			std::cout << "Error adding array to HDF5 file" << std::endl;
			return -1;
		} else {
			std::cout << "Array appended to HDF5 file" << std::endl;
		}
	} else if (datatype == 1) {

		boost::shared_ptr< hoNDArray< std::complex<float> > > A = read_nd_array< std::complex<float> >(arrayfilename.c_str());

		if (hoNDArray_hdf5_append(A.get(), hdf5filename.c_str(), varname.c_str()) < 0) {
			std::cout << "Error adding array to HDF5 file" << std::endl;
			return -1;
		} else {
			std::cout << "Array appended to HDF5 file" << std::endl;
		}
	} else if (datatype == 2) {

		boost::shared_ptr< hoNDArray< unsigned short > > A = read_nd_array< unsigned short >(arrayfilename.c_str());

		if (hoNDArray_hdf5_append(A.get(), hdf5filename.c_str(), varname.c_str()) < 0) {
			std::cout << "Error adding array to HDF5 file" << std::endl;
			return -1;
		} else {
			std::cout << "Array appended to HDF5 file" << std::endl;
		}
	} else {
		std::cout << "Unknown type: " << datatype << std::endl;
		return -1;
	}

	return 0;
}
