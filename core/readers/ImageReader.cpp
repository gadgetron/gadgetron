
#include "ImageReader.h"
#include "MessageID.h"

#include "io/primitives.h"

namespace {
    using namespace Gadgetron;
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


    using image_datatypes = Core::variant<unsigned short, unsigned int, int, float, double, std::complex<float>,std::complex<double>>;
    std::map<uint16_t,image_datatypes> ismrmrd_to_variant = {
        {ISMRMRD::ISMRMRD_USHORT, (unsigned short)0},
        {ISMRMRD::ISMRMRD_UINT, (unsigned int)0},
        {ISMRMRD::ISMRMRD_INT, (int)0},
        {ISMRMRD::ISMRMRD_FLOAT, (float)0},
        {ISMRMRD::ISMRMRD_DOUBLE, (double)0},
        {ISMRMRD::ISMRMRD_CXFLOAT, std::complex<float>(0)},
        {ISMRMRD::ISMRMRD_CXDOUBLE, std::complex<double>(0)}
    };

    template<class T>
    Core::Message read_image_message(std::istream& stream, ISMRMRD::ImageHeader header, Core::optional<ISMRMRD::MetaContainer> meta, T type_tag){
        auto image_data = hoNDArray<T>(header.matrix_size[0],header.matrix_size[1],header.matrix_size[2],header.channels);
        Gadgetron::Core::IO::read(stream,image_data.data(),image_data.size());
        return Core::Message(header,std::move(image_data),std::move(meta));
    }

    Core::optional<ISMRMRD::MetaContainer> parse_meta(const std::string &serialized_meta) {

        if (serialized_meta.empty()) return Core::none;

        ISMRMRD::MetaContainer meta;
        ISMRMRD::deserialize(serialized_meta.c_str(), meta);

        return meta;
    }
}

Gadgetron::Core::Message Gadgetron::Core::Readers::ImageReader::read(std::istream& stream) {

    auto header = IO::read<ISMRMRD::ImageHeader>(stream);
    auto serialized_meta = IO::read_string_from_stream<uint64_t>(stream);

    auto meta = parse_meta(serialized_meta);

    auto datatype = ismrmrd_to_variant.at(header.data_type);
    return Core::visit([&](auto tag){return read_image_message(stream,header,std::move(meta),tag);}, datatype);
}

uint16_t Gadgetron::Core::Readers::ImageReader::slot() {
    return GADGET_MESSAGE_ISMRMRD_IMAGE;
}


namespace Gadgetron::Core::Readers{
    GADGETRON_READER_EXPORT(ImageReader)
}