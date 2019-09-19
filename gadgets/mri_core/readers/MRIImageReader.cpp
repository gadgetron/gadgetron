#include "MRIImageReader.h"
#include <ismrmrd/ismrmrd.h>
#include <ismrmrd/meta.h>
#include "io/primitives.h"
#include "MessageID.h"

namespace Gadgetron {

    namespace {
        template<class T>
        Core::Message
        combine_and_read(std::istream &stream, ISMRMRD::ImageHeader header,
                         Core::optional<ISMRMRD::MetaContainer> meta) {

            using namespace Gadgetron::Core;
            auto array = hoNDArray<T>(
                    header.matrix_size[0],
                    header.matrix_size[1],
                    header.matrix_size[2],
                    header.channels
            );

            IO::read(stream, array.data(), array.size());
            return Core::Message(std::move(header), std::move(array), std::move(meta));
        }

        using ISMRMRD_TYPES = Core::variant<uint16_t, int16_t, uint32_t, int32_t, float, double, std::complex<float>, std::complex<double>>;
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

            stream.read(buffer.get(), meta_attrib_length);

            meta = ISMRMRD::MetaContainer();

            ISMRMRD::deserialize(buffer.get(), *meta);
        }

        return Core::visit([&](auto type_tag) {
            return combine_and_read<decltype(type_tag)>(stream, std::move(header), std::move(meta));
        }, ismrmrd_type_map.at(header.data_type));
    }

    uint16_t MRIImageReader::slot() {
        return Core::MessageID::GADGET_MESSAGE_ISMRMRD_IMAGE;
    }

    GADGETRON_READER_EXPORT(MRIImageReader)
}
