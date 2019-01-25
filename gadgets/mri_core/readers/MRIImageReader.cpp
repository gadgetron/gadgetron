#include "MRIImageReader.h"
#include <ismrmrd/ismrmrd.h>
#include <ismrmrd/meta.h>
#include "io/primitives.h"
#include <boost/variant.hpp>
/*


class GadgetronClientImageMessageReader : public GadgetronClientMessageReader
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
    void read_data_attrib(tcp::socket* stream, const ISMRMRD::ImageHeader& h, ISMRMRD::Image<T>& im)
    {
        im.setHead(h);

        typedef unsigned long long size_t_type;

        //Read meta attributes
        size_t_type meta_attrib_length;
        boost::asio::read(*stream, boost::asio::buffer(&meta_attrib_length, sizeof(size_t_type)));

        if (meta_attrib_length>0)
        {
            std::string meta_attrib(meta_attrib_length, 0);
            boost::asio::read(*stream, boost::asio::buffer(const_cast<char*>(meta_attrib.c_str()), meta_attrib_length));
            im.setAttributeString(meta_attrib);
        }

        //Read image data
        boost::asio::read(*stream, boost::asio::buffer(im.getDataPtr(), im.getDataSize()));
        {
            if (!dataset_) {

                {
                    mtx.lock();
                    dataset_ = boost::shared_ptr<ISMRMRD::Dataset>(new ISMRMRD::Dataset(file_name_.c_str(), group_name_.c_str(), true)); // create if necessary 
                    mtx.unlock();
                }
            }

            std::stringstream st1;
            st1 << "image_" << h.image_series_index;
            std::string image_varname = st1.str();

            {
                mtx.lock();
                //TODO should this be wrapped in a try/catch?
                dataset_->appendImage(image_varname, im);
                mtx.unlock();
            }
        }
    }

    virtual void read(tcp::socket* stream) 
    {
        //Read the image headerfrom the socket
        ISMRMRD::ImageHeader h;
        boost::asio::read(*stream, boost::asio::buffer(&h,sizeof(ISMRMRD::ImageHeader)));

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
            throw GadgetronClientException("Invalide image data type ... ");
        }
    }



 */

namespace Gadgetron {

    namespace {
        template<class T>
        Core::Message
        combine_and_read(std::istream &stream, ISMRMRD::ImageHeader header,
                         Core::optional<ISMRMRD::MetaContainer> meta) {

            using namespace Gadgetron::Core;
            auto array = std::make_unique<hoNDArray<T>>(header.matrix_size[0], header.matrix_size[1],
                                                        header.matrix_size[2], header.channels);
            IO::read(stream, *array);

            return Core::Message(std::move(header), std::move(array), std::move(meta));


        }

        using ISMRMRD_TYPES = boost::variant<uint16_t, int16_t, uint32_t, int32_t, float, double, std::complex<float>, std::complex<double>>;
        static const auto ismrmrd_type_map = std::unordered_map<uint16_t, ISMRMRD_TYPES>{
                {ISMRMRD::ISMRMRD_USHORT,   uint16_t()},
                {ISMRMRD::ISMRMRD_SHORT,    int16_t()},
                {ISMRMRD::ISMRMRD_INT,      int32_t()},
                {ISMRMRD::ISMRMRD_UINT,     uint32_t()},
                {ISMRMRD::ISMRMRD_FLOAT,    float()},
                {ISMRMRD::ISMRMRD_DOUBLE,   double()},
                {ISMRMRD::ISMRMRD_CXFLOAT,  std::complex<float>()},
                {ISMRMRD::ISMRMRD_CXDOUBLE, std::complex<double>()}
        };

    }


    Core::Message MRIImageReader::read(std::istream &stream) {
        using namespace Gadgetron::Core;

        auto header = IO::read<ISMRMRD::ImageHeader>(stream);
        typedef unsigned long long size_t_type;

        //Read meta attributes

        auto meta_attrib_length = IO::read<size_t>(stream);

        optional<ISMRMRD::MetaContainer> meta;
        if (meta_attrib_length > 0) {
            auto buffer = std::make_unique<char[]>(meta_attrib_length + 1);

            stream.read(buffer.get(), meta_attrib_length + 1);

            meta = ISMRMRD::MetaContainer();

            ISMRMRD::deserialize(buffer.get(), *meta);
        }

        return boost::apply_visitor([&](auto type_tag) {
            return combine_and_read<decltype(type_tag)>(stream, std::move(header), std::move(meta));
        }, ismrmrd_type_map.at(header.data_type));
    }


    uint16_t MRIImageReader::slot() {
        return 1009;
    }
  GADGETRON_READER_EXPORT(MRIImageReader)
}
