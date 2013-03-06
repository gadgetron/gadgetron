/*
 * hdf5_utils.h
 *
 *  Created on: Dec 3, 2012
 *      Author: Dae
 */

#pragma once
#include "hoNDArray.h"
#include "hdf5.h"
#include "hdf5_hl.h"
#include <stdlib.h>

namespace Gadgetron{


template<unsigned int D> void saveNDArray2HDF5(hoNDArray<float>* input,std::string filename,vector_td<float,D> dimensions, vector_td<float,D> origin, std::string arguements, unsigned int iterations){
	/* Create a new file using default properties. */
	   hid_t file_id = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

	   hsize_t dims[D];

	   for (int i = 0; i < D; i++) dims[i] = input->get_size(D-i-1);
	   H5LTmake_dataset (file_id, "/image", D, dims, H5T_NATIVE_FLOAT, input->get_data_ptr());
	   H5LTset_attribute_string(file_id,"/image","input",arguements.c_str());

	   H5LTset_attribute_float(file_id,"/image","dimensions",dimensions.vec,D);
	   H5LTset_attribute_float(file_id,"/image","origin",origin.vec,D);
	   H5LTset_attribute_uint(file_id,"/image","iterations",&iterations,1);


	   H5Fclose(file_id);
}
}

