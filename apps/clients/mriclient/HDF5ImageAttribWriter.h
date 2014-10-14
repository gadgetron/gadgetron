/*
*  HDF5ImageAttribWriter.h
*
*  Created on: April 27, 2014
*  Author: Hui Xue
*/

#ifndef HDF5IMAGEATTRIBWRITER_H_
#define HDF5IMAGEATTRIBWRITER_H_

#include "ImageAttribWriter.h"

#include <ismrmrd/dataset.h>
#include <sstream>

namespace Gadgetron
{
    template <typename T> class HDF5ImageAttribWriter : public ImageAttribWriter<T>
    {
    public:

        typedef ImageAttribWriter<T> BaseClass;
        typedef typename BaseClass::size_t_type size_t_type;

        HDF5ImageAttribWriter(std::string filename, std::string groupname, ACE_Thread_Mutex& mtx, std::string prefix=std::string())
            : ImageAttribWriter<T>()
            , file_name_(filename)
            , group_name_(groupname)
            , prefix_(prefix)
            , mtx_(&mtx)
            , dataset_(filename.c_str(), groupname.c_str())
        {

        }

        virtual int process_image(ISMRMRD::ImageHeader* img_head, hoNDArray< T >* data, ISMRMRD::MetaContainer* img_attrib)
        {
            try
            {
                size_t n;

                // image data role
                std::vector<std::string> dataRole;

                size_t num = img_attrib->length(GTPLUS_DATA_ROLE);

                if ( num == 0 )
                {
                    dataRole.push_back("Image");
                }
                else
                {
                    dataRole.resize(num);
                    for ( n=0; n<num; n++ )
                    {
                        dataRole[n] = std::string( img_attrib->as_str(GTPLUS_DATA_ROLE, n) );
                    }
                }

                long imageNumber;
                imageNumber = img_attrib->as_long(GTPLUS_IMAGENUMBER, 0);

                long cha, slc, e2, con, phs, rep, set, ave;
                cha = img_attrib->as_long(GTPLUS_CHA,        0);
                slc = img_attrib->as_long(GTPLUS_SLC,        0);
                e2  = img_attrib->as_long(GTPLUS_E2,         0);
                con = img_attrib->as_long(GTPLUS_CONTRAST,   0);
                phs = img_attrib->as_long(GTPLUS_PHASE,      0);
                rep = img_attrib->as_long(GTPLUS_REP,        0);
                set = img_attrib->as_long(GTPLUS_SET,        0);
                ave = img_attrib->as_long(GTPLUS_AVERAGE,    0);

                std::ostringstream ostr;

                if ( !prefix_.empty() )
                {
                    ostr << prefix_ << "_";
                }

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

                ACE_DEBUG( (LM_DEBUG, ACE_TEXT("Writing image %s\n"), filename.c_str()) );
		// Otherwise the above should be removed

                char* buf = NULL;
                size_t_type len(0);

                try
                {
                    std::stringstream str;
                    ISMRMRD::serialize( *img_attrib, str);
                    std::string attribContent = str.str();
                    len = attribContent.length()+1;

                    buf = new char[len];
                    GADGET_CHECK_THROW(buf != NULL);

                    memset(buf, '\0', sizeof(char)*len);
                    memcpy(buf, attribContent.c_str(), len-1);
                }
                catch(...)
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
		ISMRMRD::Image<T> img;
		img.setHead(*img_head);
		img.setAttributeString(attrib);
                memcpy(img.getData(), data->get_data_ptr(), img.getDataSize());

                {
                    ACE_GUARD_RETURN(ACE_Thread_Mutex, guard, *mtx_, -1);
                    try {
                        dataset_.appendImage(image_varname, ISMRMRD::ISMRMRD_BLOCKMODE_ARRAY, img);
                    }
                    catch (...)
                    {
                        GADGET_DEBUG1("Failed to write image.\n");
                        return GADGET_FAIL;
                    }
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

        ACE_Thread_Mutex* mtx_;
    };
}

#endif /* HDF5IMAGEATTRIBWRITER_H_ */
