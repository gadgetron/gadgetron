/*
 * HDF5ImageWriter.h
 *
 *  Created on: Jan 25, 2012
 *      Author: Michael S. Hansen
 */

#ifndef HDF5IMAGEWRITER_H_
#define HDF5IMAGEWRITER_H_

#include "ImageWriter.h"

#include <ismrmrd_hdf5.h>
#include <sstream>

namespace Gadgetron{
template <typename T> class HDF5ImageWriter : public ImageWriter<T>
{

public:
	HDF5ImageWriter(std::string filename, std::string groupname)
	: ImageWriter<T>()
	, file_name_(filename)
	, group_name_(groupname)
	, dataset_(filename.c_str(), groupname.c_str())
	{

	}

	virtual int process_image(ISMRMRD::ImageHeader* img_head,
			hoNDArray< T >* data)
	{
		try {
			ISMRMRD::HDF5Exclusive lock; //This will ensure threadsafe access to HDF5
			std::stringstream st1;
			st1 << "image_" << img_head->image_series_index << ".head";
			std::string head_varname = st1.str();

			std::stringstream st2;
			st2 << "image_" << img_head->image_series_index << ".img";
			std::string img_varname = st2.str();

			if (dataset_.appendImageHeader(*img_head, head_varname.c_str()) < 0) {
				GADGET_DEBUG1("Failed to write image header\n");
				return GADGET_FAIL;
			}

            std::vector<unsigned long long> dim = *data->get_dimensions();
            std::vector<unsigned int> dim2(dim.size());

            unsigned long long ii;
            for ( ii=0; ii<dim.size(); ii++ )
            {
                dim2[ii] = dim[ii];
            }

			if (dataset_.appendArray(dim2,data->get_data_ptr(), img_varname.c_str())  < 0) {
				GADGET_DEBUG1("Failed to write image data\n");
				return GADGET_FAIL;
			};
		} catch (...) {
			GADGET_DEBUG1("Error attempting to append images to HDF5 file\n");
			return GADGET_FAIL;
		}

		return GADGET_OK;
	}

protected:
	std::string group_name_;
	std::string file_name_;
	ISMRMRD::IsmrmrdDataset dataset_;
};

}
#endif /* HDF5IMAGEWRITER_H_ */
