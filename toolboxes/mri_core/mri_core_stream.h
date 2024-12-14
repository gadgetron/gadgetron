/** \file   mri_core_stream.h
    \brief  Implementation a class to help stream out the reconstruction intermediate data and results

    The images, headers, and hoNDArrays are streamed out using the MRD protocol.

    ToBeImplemented: The waveforms are streamed out as the MRD Waveform.

    \author Hui Xue
*/

#pragma once

#include <fstream>

#include "hoNDArray.h"
#include "mri_core_def.h"
#include "io/primitives.h"

#include "mrd/binary/protocols.h"


namespace Gadgetron
{
    class GenericReconMrdStreamer
    {
    public:

        GenericReconMrdStreamer();
        GenericReconMrdStreamer(const std::map<std::string, std::string>& parameters);

        ~GenericReconMrdStreamer();

        bool verbose_;

        void initialize_stream_name_buffer(const std::map<std::string, std::string>& parameters);
        void initialize_stream_name_buffer(const std::map<std::string, std::string>& parameters, const std::string& name);

        void close_stream_buffer();

        void stream_mrd_header(const mrd::Header& hdr);

        void stream_mrd_waveforms(const std::vector<mrd::WaveformUint32>& wavs);

        // stream of ND array buffer
        template <typename DataType>
        void stream_to_array_buffer(const std::string& name, const hoNDArray<DataType>& data)
        {
            std::string buf_name;
            auto writer = this->find_and_open_stream(name, buf_name);
            if (writer)
            {
                GDEBUG_STREAM("Generic recon, continue streaming the array data to the buffer " << buf_name);
                writer->WriteData(mrd::ArrayComplexFloat(data));
            }
            else
            {
                GWARN_CONDITION_STREAM(this->verbose_, "Generic reconstruction chain, the pre-set buffer names do not include " << name << "; the data will not be saved into the buffer ...");
            }
        }

        // stream to mrd image stream, e.g. for coil maps, g-maps and reconed images
        template <typename DataType>
        void stream_to_mrd_image_buffer(const std::string& name, const hoNDArray<DataType>& img, const hoNDArray< mrd::ImageHeader >& headers, const hoNDArray< mrd::ImageMeta >& meta)
        {
            std::string buf_name;
            auto writer = find_and_open_stream(name, buf_name);

            if (writer) {
                // convert images to one or more mrd Images
                std::vector< mrd::Image<DataType> > mrd_images;

                size_t RO = img.get_size(0);
                size_t E1 = img.get_size(1);
                size_t E2 = img.get_size(2);
                size_t CHA = img.get_size(3);

                size_t N = img.get_size(4);
                size_t S = img.get_size(5);
                size_t SLC = img.get_size(6);

                GDEBUG_STREAM("Generic recon, convert recon images to mrd images for " << name << " [RO E1 E2 CHA N S SLC] = [" << RO << " " << E1 << " " << E2 << " " << CHA << " " << N << " " << S << " " << SLC << "] - " << buf_name);

                mrd_images.resize(N*S*SLC);

                for (auto slc=0; slc<SLC; slc++)
                {
                    for (auto s=0; s<S; s++)
                    {
                        for (auto n=0; n<N; n++)
                        {
                            size_t ind = n+s*N+slc*N*S;

                            mrd::Image<DataType>& a_img = mrd_images[ind];

                            a_img.data.create(RO, E1, E2, CHA);
                            memcpy(a_img.data.data(), &img(0, 0, 0, 0, n, s, slc), sizeof(DataType)*RO*E1*E2*CHA);

                            a_img.head = headers(n, s, slc);
                            a_img.meta = meta(n, s, slc);
                        }
                    }
                }

                GDEBUG_STREAM("Generic recon, continue streaming the data to the mrd image buffer " << buf_name << " for " << mrd_images.size() << " images ...");
                for (auto im : mrd_images)
                {
                    this->buffer_names_[name].second->WriteData(im);
                }
            }
            else
            {
                GWARN_CONDITION_STREAM(this->verbose_, "Generic reconstruction chain, the pre-set buffer names do not include " << name << "; the data will not be saved into the buffer ...");
            }
        }

    protected:
        std::map<std::string, std::pair<std::string, std::shared_ptr<mrd::binary::MrdWriter> > > buffer_names_;
        std::shared_ptr<mrd::binary::MrdWriter> find_and_open_stream(const std::string& name, std::string& buf_name, std::optional<mrd::Header> header = std::nullopt);
    };
}