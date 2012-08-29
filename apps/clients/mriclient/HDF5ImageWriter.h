/*
 * HDF5ImageWriter.h
 *
 *  Created on: Jan 25, 2012
 *      Author: Michael S. Hansen
 */

#ifndef HDF5IMAGEWRITER_H_
#define HDF5IMAGEWRITER_H_

#include "ImageWriter.h"
#include "hoNDArray_hdf5_io.h"
#include "mri_hdf5_io.h"
#include <sstream>

template <typename T> class HDF5ImageWriter : public ImageWriter<T>
{

public:
	HDF5ImageWriter(std::string filename, std::string groupname)
	: ImageWriter<T>()
	, file_name_(filename)
	, group_name_(groupname)
	{

	}

	virtual int process_image(ISMRMRD::ImageHeader* img_head,
			hoNDArray< T >* data)
	{
		try {
		HDF5Exclusive lock; //This will ensure threadsafe access to HDF5
	    std::stringstream st;
	    st << img_head->image_series_index;
		std::string varname = group_name_ + std::string("/") + std::string("data_") + st.str();

		if (!(hdf5_append_array(data, file_name_.c_str(), varname.c_str()) == 0)) {
			GADGET_DEBUG1("File is not good for writing\n");
			return GADGET_FAIL;
		}

		varname = group_name_ + std::string("/") + std::string("header_") + st.str();

		if (!(hdf5_append_struct(img_head, file_name_.c_str(), varname.c_str()) == 0)) {
			GADGET_DEBUG1("File is not good for writing\n");
			return GADGET_FAIL;
		}

		/*
		varname = group_name_ + std::string("/") + std::string("comb_") + st.str();

		if (!(hdf5_append_struct_with_data(img_head, data, file_name_.c_str(), varname.c_str()) == 0)) {
			GADGET_DEBUG1("File is not good for writing\n");
			return GADGET_FAIL;
		}
		*/
		} catch (...) {
			GADGET_DEBUG1("Error attempting to append images to HDF5 file\n");
			return GADGET_FAIL;
		}

		return GADGET_OK;
	}

protected:
	std::string group_name_;
	std::string file_name_;
};


#endif /* HDF5IMAGEWRITER_H_ */
