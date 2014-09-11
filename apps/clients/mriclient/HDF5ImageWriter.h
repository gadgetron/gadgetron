/*
* HDF5ImageWriter.h
*
*  Created on: Jan 25, 2012
*      Author: Michael S. Hansen
*/

#ifndef HDF5IMAGEWRITER_H_
#define HDF5IMAGEWRITER_H_

#include "ImageWriter.h"

#include <ismrmrd_dataset.h>
#include <sstream>

namespace Gadgetron
{
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
                std::stringstream st1;
                st1 << "image_" << img_head->image_series_index;
                std::string image_varname = st1.str();

		// TODO this makes a copy of the data
		// what's the best way to do it without copies?
		ISMRMRD::Image img;
                img.setHead(*img_head);
                memcpy(img.getData(), data->get_data_ptr(), img.getDataSize());

                if (dataset_.appendImage(image_varname, ISMRMRD::ISMRMRD_BLOCKMODE_ARRAY, img) < 0) {
                    GADGET_DEBUG1("Failed to write image.\n");
                    return GADGET_FAIL;
                }

            } catch (...) {
                GADGET_DEBUG1("Error attempting to append images to HDF5 file\n");
                return GADGET_FAIL;
            }

            return GADGET_OK;
        }

    protected:
        std::string group_name_;
        std::string file_name_;
        ISMRMRD::Dataset dataset_;
    };
}

#endif /* HDF5IMAGEWRITER_H_ */
