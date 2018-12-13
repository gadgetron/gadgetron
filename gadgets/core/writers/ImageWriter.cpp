
#include <ismrmrd/meta.h>
#include <boost/optional.hpp>

#include "ImageWriter.h"

namespace {

    using namespace Gadgetron;
    using namespace Gadgetron::Core;

    template<class T>
class TypedImageWriter : public TypedWriter<ISMRMRD::ImageHeader, hoNDArray<T>, boost::optional<ISMRMRD::MetaContainer>> {
    public:
        void serialize(
                std::ostream &stream,
                std::unique_ptr<ISMRMRD::ImageHeader> &&header,
                std::unique_ptr<hoNDArray<T>> &&data,
                std::unique_ptr<boost::optional<ISMRMRD::MetaContainer>> &&meta
        ) override {

            uint16_t message_id = 1022;
            stream.write(reinterpret_cast<char *>(&message_id), sizeof(message_id));
            stream.write(reinterpret_cast<char *>(header.get()), sizeof(*header));

            std::string serialized_meta;

            if(*meta) {
                std::stringstream meta_stream;

                ISMRMRD::serialize(meta->value(), meta_stream);

                serialized_meta = meta_stream.str();
            }

            uint64_t meta_size = serialized_meta.size();
            stream.write(reinterpret_cast<char *>(&meta_size), sizeof(meta_size));
            stream.write(serialized_meta.c_str(), meta_size);

            stream.write(reinterpret_cast<char *>(data->get_data_ptr()), data->get_number_of_bytes());
        }
    };

    const std::vector<std::shared_ptr<Writer>> writers {
        std::make_shared<TypedImageWriter<float>>(),
        std::make_shared<TypedImageWriter<double>>(),
        std::make_shared<TypedImageWriter<std::complex<float>>>(),
        std::make_shared<TypedImageWriter<std::complex<double>>>(),
        std::make_shared<TypedImageWriter<unsigned short>>(),
        std::make_shared<TypedImageWriter<short>>(),
        std::make_shared<TypedImageWriter<int>>(),
         std::make_shared<TypedImageWriter<unsigned int>>()
    };
}


namespace Gadgetron::Core::Writers {

    bool ImageWriter::accepts(const Message &message) {
        return std::any_of(writers.begin(), writers.end(),
                [&](auto &writer) { return writer->accepts(message); }
        );
    }

    void ImageWriter::write(std::ostream &stream, std::unique_ptr<Message> &&message) {
        for (auto &writer : writers){
            if (writer->accepts(*message)) {
                writer->write(stream, std::move(message));
                return;
            }
        }
    }

    GADGETRON_WRITER_EXPORT(ImageWriter)
}

