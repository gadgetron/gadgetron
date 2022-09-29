#include <string>
#include <memory>
#include <mutex>
#include <variant>
#include <iostream>

#include <boost/asio.hpp>
#include <boost/filesystem.hpp>

#include <ismrmrd/ismrmrd.h>
#include <ismrmrd/dataset.h>
#include <ismrmrd/meta.h>
#include <ismrmrd/xml.h>
#include <ismrmrd/waveform.h>

#include "MessageID.h"

namespace
{

class GadgetronClientImageMessageReader
{

public:
    GadgetronClientImageMessageReader(std::string filename, std::string groupname)
        : file_name_(filename)
        , group_name_(groupname)
    {

    }

    ~GadgetronClientImageMessageReader() {
    }

    template <typename T>
    void read_data_attrib(std::istream& stream, const ISMRMRD::ImageHeader& h, ISMRMRD::Image<T>& im)
    {
        im.setHead(h);

        typedef unsigned long long size_t_type;

        //Read meta attributes
        size_t_type meta_attrib_length;
        stream.read(reinterpret_cast<char*>(&meta_attrib_length), sizeof(size_t_type));

        if (meta_attrib_length>0)
        {
            std::vector<char> data(meta_attrib_length);
            stream.read(data.data(), meta_attrib_length);
            std::string meta_attrib(data.data(), data.size());
            im.setAttributeString(meta_attrib);
        }

        //Read image data
        stream.read(reinterpret_cast<char*>(im.getDataPtr()), im.getDataSize());
        {
            if (!dataset_) {

                {
                    dataset_ = std::shared_ptr<ISMRMRD::Dataset>(new ISMRMRD::Dataset(file_name_.c_str(), group_name_.c_str(), true)); // create if necessary
                }
            }

            std::stringstream st1;
            st1 << "image_" << h.image_series_index;
            std::string image_varname = st1.str();

            {
                dataset_->appendImage(image_varname, im);
            }
        }
    }

    void read(std::istream& stream)
    {
        ISMRMRD::ImageHeader h;
        stream.read(reinterpret_cast<char*>(&h),sizeof(ISMRMRD::ImageHeader));

        if (h.data_type == ISMRMRD::ISMRMRD_USHORT)
        {
            ISMRMRD::Image<unsigned short> im;
            this->read_data_attrib(stream, h, im);
        }
        else if (h.data_type == ISMRMRD::ISMRMRD_SHORT)
        {
            ISMRMRD::Image<short> im;
            this->read_data_attrib(stream, h, im);
        }
        else if (h.data_type == ISMRMRD::ISMRMRD_UINT)
        {
            ISMRMRD::Image<unsigned int> im;
            this->read_data_attrib(stream, h, im);
        }
        else if (h.data_type == ISMRMRD::ISMRMRD_INT)
        {
            ISMRMRD::Image<int> im;
            this->read_data_attrib(stream, h, im);
        }
        else if (h.data_type == ISMRMRD::ISMRMRD_FLOAT)
        {
            ISMRMRD::Image<float> im;
            this->read_data_attrib(stream, h, im);
        }
        else if (h.data_type == ISMRMRD::ISMRMRD_DOUBLE)
        {
            ISMRMRD::Image<double> im;
            this->read_data_attrib(stream, h, im);
        }
        else if (h.data_type == ISMRMRD::ISMRMRD_CXFLOAT)
        {
            ISMRMRD::Image< std::complex<float> > im;
            this->read_data_attrib(stream, h, im);
        }
        else if (h.data_type == ISMRMRD::ISMRMRD_CXDOUBLE)
        {
            ISMRMRD::Image< std::complex<double> > im;
            this->read_data_attrib(stream, h, im);
        }
        else
        {
            throw std::runtime_error("Invalid image data type ... ");
        }
    }

protected:
    std::string group_name_;
    std::string file_name_;
    std::shared_ptr<ISMRMRD::Dataset> dataset_;
};

}


class ImgMsgAdapter
{
public:
    explicit ImgMsgAdapter() = default;

    void convert(
        std::istream& input_stream,
        const std::string& ismrmrd_filepath,
        const std::string& output_group = "/dataset")
    {
        input_stream.exceptions(std::istream::failbit | std::istream::badbit | std::istream::eofbit);

        read_messages(input_stream, ismrmrd_filepath, output_group);
    }

private:
    void read_ismrmrd_header(std::istream& input_stream)
    {
        ISMRMRD::IsmrmrdHeader hdr;

        uint32_t hdr_size = 0;
        input_stream.read(reinterpret_cast<char*>(&hdr_size), sizeof(uint32_t));
        if(hdr_size > 0)
        {
            std::vector<char> data(hdr_size);
            input_stream.read(data.data(), hdr_size);
            ISMRMRD::deserialize(std::string(data.data(), data.size()).c_str(), hdr);
        }
    }

    void read_messages(std::istream& input_stream, const std::string& ismrmrd_filepath, const std::string& output_group)
    {
        auto output = GadgetronClientImageMessageReader(ismrmrd_filepath, output_group);
        bool closed = false;

        while (!input_stream.eof() && !closed)
        {
            ::Gadgetron::Core::MessageID id = ::Gadgetron::Core::MessageID::CLOSE;
            input_stream.read(reinterpret_cast<char*>(&id), sizeof(::Gadgetron::Core::MessageID));

            switch(id)
            {
                case ::Gadgetron::Core::MessageID::HEADER:
                {
                    read_ismrmrd_header(input_stream);
                    break;
                }
                case ::Gadgetron::Core::MessageID::GADGET_MESSAGE_ISMRMRD_IMAGE:
                {
                    output.read(input_stream);
                    break;
                }
                case ::Gadgetron::Core::MessageID::CLOSE:
                {
                    closed = true;
                    break;
                }
                case ::Gadgetron::Core::MessageID::ERROR:
                {
                    throw std::runtime_error("Got error while processing input stream");
                }
                default:
                {
                    throw std::runtime_error("Unsupported message ID: " + std::to_string(id));
                }
            }
        }
    }
};
