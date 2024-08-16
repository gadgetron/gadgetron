

/** \file   mri_core_stream.h
    \brief  Implementation a class to help stream out the reconstruction intermediate data and results

    The hoNDArray is streamed out as the IO class utility function does. 

    The images and headers are streamed out as  the ISMRMRD image and header.

    ToBeImplemented: The waveforms are streamed out as the ISMRMRD waveform. 

    All ismrmrd structures follow the convention to use the ISMRMRD::ProtocolSerializer and ProtocolDeserializer.

    \author Hui Xue
*/

#pragma once

#include <fstream>
#include "ismrmrd/ismrmrd.h"
#include "ismrmrd/xml.h"
#include "ismrmrd/meta.h"
#include "ismrmrd/serialization.h"
#include "ismrmrd/serialization_iostream.h"
#include "hoNDArray.h"
#include "mri_core_def.h"
#include "mri_core_data.h"
#include "io/primitives.h"
#include "io/ismrmrd_types.h"

namespace Gadgetron 
{
    template <typename T> 
    void convert_hoNDArray_to_ismrmrd_ndarray(const hoNDArray<T>& ho_arr, ISMRMRD::NDArray<T>& arr);

    template <typename T> 
    void convert_ismrmrd_ndarray_to_hoNDArray(const ISMRMRD::NDArray<T>& arr, hoNDArray<T>& ho_arr);

    class GenericReconIsmrmrdStreamer
    {
    public:

        GenericReconIsmrmrdStreamer();
        GenericReconIsmrmrdStreamer(const std::map<std::string, std::string>& parameters);

        ~GenericReconIsmrmrdStreamer();

        bool verbose_;

        void initialize_stream_name_buffer(const std::map<std::string, std::string>& parameters);
        void initialize_stream_name_buffer(const std::map<std::string, std::string>& parameters, const std::string& name);

        void close_stream_buffer();

        void stream_ismrmrd_header(const ISMRMRD::IsmrmrdHeader& hdr);

        void stream_ismrmrd_waveform(const std::vector< ISMRMRD::Waveform>& wav);

        // stream of ND array buffer
        template <typename DataType> 
        void stream_to_array_buffer(const std::string& name, const hoNDArray<DataType>& data)
        {
            std::string buf_name;
            std::shared_ptr<std::ofstream> os = this->find_and_open_stream(name, buf_name);
            if (os && os->is_open())
            {
                GDEBUG_STREAM("Generic recon, continue streaming the array data to the buffer " << buf_name);
                ISMRMRD::NDArray<DataType> arr;
                convert_hoNDArray_to_ismrmrd_ndarray(data, arr);

                ISMRMRD::OStreamView ws(*os);
                ISMRMRD::ProtocolSerializer serializer(ws);
                serializer.serialize(arr);
            }
        }

        // stream to ismrmrd image stream, e.g. for coil maps, g-maps and reconed images

        template <typename DataType> 
        void stream_images(std::ofstream& os, const std::vector< ISMRMRD::Image<DataType> >& ismrmrd_images)
        {
            ISMRMRD::OStreamView ws(os);
            ISMRMRD::ProtocolSerializer serializer(ws);

            for (auto im : ismrmrd_images)
            {
                serializer.serialize(im);
            }
            os.flush();
        }

        template <typename DataType> 
        void stream_to_ismrmrd_image_buffer(const std::string& name, const hoNDArray<DataType>& img, const hoNDArray< ISMRMRD::ImageHeader >& headers, const std::vector< ISMRMRD::MetaContainer >& meta)
        {
            std::string buf_name;
            std::shared_ptr<std::ofstream> os = this->find_and_open_stream(name, buf_name);
            if (os && os->is_open())
            {
                // convert images to one or more ismrmrd images
                std::vector< ISMRMRD::Image<DataType> > ismrmrd_images;

                size_t RO = img.get_size(0);
                size_t E1 = img.get_size(1);
                size_t E2 = img.get_size(2);
                size_t CHA = img.get_size(3);

                size_t N = img.get_size(4);
                size_t S = img.get_size(5);
                size_t SLC = img.get_size(6);

                GDEBUG_STREAM("Generic recon, convert recon images to ismrmd images for " << name << " [RO E1 E2 CHA N S SLC] = [" << RO << " " << E1 << " " << E2 << " " << CHA << " " << N << " " << S << " " << SLC << "] - " << buf_name);

                ismrmrd_images.resize(N*S*SLC);

                for (auto slc=0; slc<SLC; slc++)
                {
                    for (auto s=0; s<S; s++)
                    {
                        for (auto n=0; n<N; n++)
                        {
                            size_t ind = n+s*N+slc*N*S;

                            ISMRMRD::Image<DataType>& a_img = ismrmrd_images[ind];

                            a_img.resize(RO, E1, E2, CHA);
                            memcpy(a_img.getDataPtr(), &img(0, 0, 0, 0, n, s, slc), sizeof(DataType)*RO*E1*E2*CHA);

                            ISMRMRD::ImageHeader hd = headers(n, s, slc);
                            hd.data_type = Gadgetron::Core::IO::ismrmrd_data_type<DataType>();
                            a_img.setHead(hd);

                            std::ostringstream str;
                            ISMRMRD::serialize(meta[ind], str);
                            a_img.setAttributeString(str.str());
                        }
                    }
                }

                GDEBUG_STREAM("Generic recon, continue streaming the data to the ismrmrd image buffer " << buf_name << " for " << ismrmrd_images.size() << " images ...");
                stream_images(*os, ismrmrd_images);
            }
        }

    protected:
        std::map<std::string, std::pair<std::string, std::shared_ptr<std::ofstream> > > buffer_names_;

        std::shared_ptr<std::ofstream> find_and_open_stream(const std::string& name, std::string& buf_name);
    };
}