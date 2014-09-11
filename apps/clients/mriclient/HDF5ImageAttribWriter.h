/*
*  HDF5ImageAttribWriter.h
*
*  Created on: April 27, 2014
*  Author: Hui Xue
*/

#ifndef HDF5IMAGEATTRIBWRITER_H_
#define HDF5IMAGEATTRIBWRITER_H_

#include "ImageAttribWriter.h"

#include <ismrmrd_dataset.h>
#include <sstream>

namespace Gadgetron
{
    template <typename T> class HDF5ImageAttribWriter : public ImageAttribWriter<T>
    {
    public:

        typedef ImageAttribWriter<T> BaseClass;
        typedef typename BaseClass::size_t_type size_t_type;

        HDF5ImageAttribWriter(std::string filename, std::string groupname, std::string prefix=std::string())
            : ImageAttribWriter<T>()
            , file_name_(filename)
            , group_name_(groupname)
            , prefix_(prefix)
            , dataset_(filename.c_str(), groupname.c_str())
        {

        }

        virtual int process_image(ISMRMRD::ImageHeader* img_head, hoNDArray< T >* data, GtImageAttribType* img_attrib)
        {
            try
            {
                // image data role
                std::vector<std::string> dataRole;
                if ( !img_attrib->attributeString_.get(GTPLUS_DATA_ROLE, dataRole) )
                {
                    dataRole.push_back("Image");
                }

                long long imageNumber;
                img_attrib->attributeInteger_.get(GTPLUS_IMAGENUMBER, 0, imageNumber);

                long long cha, slc, e2, con, phs, rep, set, ave;
                img_attrib->attributeInteger_.get(GTPLUS_CHA,        0, cha);
                img_attrib->attributeInteger_.get(GTPLUS_SLC,        0, slc);
                img_attrib->attributeInteger_.get(GTPLUS_E2,         0, e2);
                img_attrib->attributeInteger_.get(GTPLUS_CONTRAST,   0, con);
                img_attrib->attributeInteger_.get(GTPLUS_PHASE,      0, phs);
                img_attrib->attributeInteger_.get(GTPLUS_REP,        0, rep);
                img_attrib->attributeInteger_.get(GTPLUS_SET,        0, set);
                img_attrib->attributeInteger_.get(GTPLUS_AVERAGE,    0, ave);

                std::ostringstream ostr;

                if ( !prefix_.empty() )
                {
                    ostr << prefix_ << "_";
                }

                size_t n;
                for ( n=0; n<dataRole.size(); n++ )
                {
                    ostr << dataRole[n] << "_";
                }

                ostr << "SLC" << slc << "_"
                     << "E2" << e2 << "_"
                     << "CON" << con << "_"
                     << "PHS" << phs << "_"
                     << "REP" << rep << "_"
                     << "SET" << set << "_" 
                     << "AVE" << ave << "_" 
                     << "num";

                std::string filename = ostr.str();

		// TODO: maybe give the user some debug info.
		// Otherwise the above should be removed

                char* buf = NULL;
                size_t_type len(0);
                if ( !img_attrib->serialize(buf, len) )
                {
                    GADGET_DEBUG1("Failed to serialize image attributes\n");
                    return GADGET_FAIL;
                }
                std::string attrib = std::string(buf+sizeof(size_t_type));
                delete [] buf;

                std::stringstream st1;
                st1 << "image_" << img_head->image_series_index;
                std::string image_varname = st1.str();

		// TODO this makes a copy of the data
		// what's the best way to do it without copies?
		ISMRMRD::Image img;
		img.setHead(*img_head);
		img.setAttributeString(attrib);
                memcpy(img.getData(), data->get_data_ptr(), img.getDataSize());

                if (dataset_.appendImage(image_varname, ISMRMRD::ISMRMRD_BLOCKMODE_ARRAY, img) < 0) {
                    GADGET_DEBUG1("Failed to write image.\n");
                    return GADGET_FAIL;
                }

            }
            catch (...)
            {
                GADGET_DEBUG1("Error attempting to append images to HDF5 file\n");
                return GADGET_FAIL;
            }

            return GADGET_OK;
        }

    protected:
        std::string group_name_;
        std::string file_name_;
        std::string prefix_;
        ISMRMRD::Dataset dataset_;
    };
}

#endif /* HDF5IMAGEATTRIBWRITER_H_ */
