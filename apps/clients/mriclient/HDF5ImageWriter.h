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

                std::vector<size_t> dim = *data->get_dimensions();
                std::vector<unsigned int> dim2(dim.size());

                size_t ii;
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
                ISMRMRD::HDF5Exclusive lock; //This will ensure threadsafe access to HDF5

                // image data role
                std::vector<std::string> dataRole;
                if ( !img_attrib->attribute4_.get(GTPLUS_DATA_ROLE, dataRole) )
                {
                    dataRole.push_back("Image");
                }

                long long imageNumber;
                img_attrib->attribute1_.get(GTPLUS_IMAGENUMBER, 0, imageNumber);

                long long cha, slc, e2, con, phs, rep, set;
                img_attrib->attribute1_.get(GTPLUS_CHA,        0, cha);
                img_attrib->attribute1_.get(GTPLUS_SLC,        0, slc);
                img_attrib->attribute1_.get(GTPLUS_E2,         0, e2);
                img_attrib->attribute1_.get(GTPLUS_CONTRAST,   0, con);
                img_attrib->attribute1_.get(GTPLUS_PHASE,      0, phs);
                img_attrib->attribute1_.get(GTPLUS_REP,        0, rep);
                img_attrib->attribute1_.get(GTPLUS_SET,        0, set);

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
                     << "num";

                std::string filename = ostr.str();

                ACE_DEBUG( (LM_DEBUG, ACE_TEXT("Writing image %s\n"), filename.c_str()) );

                std::stringstream st1;
                st1 << filename << img_head->image_index << ".head";
                std::string head_varname = st1.str();

                std::stringstream st2;
                st2 << filename << img_head->image_index << ".img";
                std::string img_varname = st2.str();

                std::stringstream st3;
                st3 << filename << img_head->image_index << ".attrib";
                std::string meta_varname = st3.str();

                if (dataset_.appendImageHeader( *img_head, head_varname.c_str()) < 0)
                {
                    GADGET_DEBUG1("Failed to write image header\n");
                    return GADGET_FAIL;
                }

                char* buf = NULL;
                size_t_type len(0);

                if ( !img_attrib->serialize(buf, len) )
                {
                    GADGET_DEBUG1("Failed to serialize image attributes\n");
                    return GADGET_FAIL;
                }

                std::string attrib = std::string(buf+sizeof(size_t_type));

                /*if (dataset_.appendImageAttrib(attrib, meta_varname.c_str()) < 0)
                {
                    GADGET_DEBUG1("Failed to write image attributes\n");
                    return GADGET_FAIL;
		}*/

                delete [] buf;

                std::vector<size_t> dim = *data->get_dimensions();
                std::vector<unsigned int> dim2(dim.size());

                size_t ii;
                for ( ii=0; ii<dim.size(); ii++ )
                {
                    dim2[ii] = dim[ii];
                }

                if (dataset_.appendArray(dim2, data->get_data_ptr(), img_varname.c_str())  < 0)
                {
                    GADGET_DEBUG1("Failed to write image data\n");
                    return GADGET_FAIL;
                }

                // still store in the conventional way
                std::stringstream st5;
                st5 << "image_" << img_head->image_series_index << ".img";
                img_varname = st5.str();

                if (dataset_.appendArray(dim2, data->get_data_ptr(), img_varname.c_str())  < 0)
                {
                    GADGET_DEBUG1("Failed to write image data\n");
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
        ISMRMRD::IsmrmrdDataset dataset_;
    };
}

#endif /* HDF5IMAGEWRITER_H_ */
