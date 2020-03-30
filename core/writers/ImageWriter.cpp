
#include <ismrmrd/meta.h>
#include <ismrmrd/ismrmrd.h>
#include <boost/optional.hpp>
#include <io/primitives.h>

#include "MessageID.h"
#include "ImageWriter.h"

namespace {

    using namespace Gadgetron;
    using namespace Gadgetron::Core;

    template<class T> inline constexpr uint16_t ismrmrd_data_type(){ return 0;}
    template<> inline constexpr uint16_t ismrmrd_data_type<unsigned short>(){return ISMRMRD::ISMRMRD_USHORT;}
    template<> inline constexpr uint16_t ismrmrd_data_type<unsigned int>(){return ISMRMRD::ISMRMRD_UINT;}
    template<> inline constexpr uint16_t ismrmrd_data_type<int>(){return ISMRMRD::ISMRMRD_INT;}
    template<> inline constexpr uint16_t ismrmrd_data_type<float>(){return ISMRMRD::ISMRMRD_FLOAT;}
    template<> inline constexpr uint16_t ismrmrd_data_type<double>(){return ISMRMRD::ISMRMRD_DOUBLE;}
    template<> inline constexpr uint16_t ismrmrd_data_type<std::complex<float>>(){return ISMRMRD::ISMRMRD_CXFLOAT;}
    template<> inline constexpr uint16_t ismrmrd_data_type<std::complex<double>>(){return ISMRMRD::ISMRMRD_CXDOUBLE;}
    template<> inline constexpr uint16_t ismrmrd_data_type<complext<float>>(){return ISMRMRD::ISMRMRD_CXFLOAT;}
    template<> inline constexpr uint16_t ismrmrd_data_type<complext<double>>(){return ISMRMRD::ISMRMRD_CXDOUBLE;}


    template<class T>
    class TypedImageWriter : public TypedWriter<ISMRMRD::ImageHeader, hoNDArray<T>, Core::optional<ISMRMRD::MetaContainer>> {
    public:
        void serialize(
                std::ostream &stream,
                const ISMRMRD::ImageHeader& header,
                const hoNDArray<T>& data,
                const optional<ISMRMRD::MetaContainer>& meta
        ) override {
            std::string serialized_meta;
            if(meta) {
                std::stringstream meta_stream;
                ISMRMRD::serialize(*meta, meta_stream);
                serialized_meta = meta_stream.str();
            }
            uint64_t meta_size = serialized_meta.size();

            auto corrected_header = header;
            corrected_header.data_type = ismrmrd_data_type<T>();
            corrected_header.attribute_string_len = meta_size;

            auto message_id = GADGET_MESSAGE_ISMRMRD_IMAGE;
            IO::write(stream, message_id);
            IO::write(stream, corrected_header);
            IO::write(stream, meta_size);
            stream.write(serialized_meta.c_str(), meta_size);
            IO::write(stream, data.get_data_ptr(), data.get_number_of_elements());
        }
    };

    const std::vector<std::shared_ptr<Writer>> writers {
        std::make_shared<TypedImageWriter<float>>(),
        std::make_shared<TypedImageWriter<double>>(),
        std::make_shared<TypedImageWriter<std::complex<float>>>(),
        std::make_shared<TypedImageWriter<std::complex<double>>>(),
        std::make_shared<TypedImageWriter<unsigned short>>(),
        std::make_shared<TypedImageWriter<short>>(),
        std::make_shared<TypedImageWriter<unsigned int>>(),
        std::make_shared<TypedImageWriter<int>>()
    };
}


namespace Gadgetron::Core::Writers {

    bool ImageWriter::accepts(const Message &message) {
        return std::any_of(writers.begin(), writers.end(),
                [&](auto &writer) { return writer->accepts(message); }
        );
    }

    void ImageWriter::write(std::ostream &stream, Message message) {
        for (auto &writer : writers){
            if (writer->accepts(message)) {
                writer->write(stream, std::move(message));
                return;
            }
        }
    }

    GADGETRON_WRITER_EXPORT(ImageWriter)
}

