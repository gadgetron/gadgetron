
#include <ismrmrd/meta.h>
#include <ismrmrd/ismrmrd.h>
#include <boost/optional.hpp>
#include <io/primitives.h>

#include "MessageID.h"
#include "ImageWriter.h"

namespace {

    using namespace Gadgetron;
    using namespace Gadgetron::Core;



    template<class T>
    class TypedImageWriter : public TypedWriter<ISMRMRD::ImageHeader, hoNDArray<T>, Core::optional<ISMRMRD::MetaContainer>> {
    public:
        void serialize(
                std::ostream &stream,
                const ISMRMRD::ImageHeader& header,
                const hoNDArray<T>& data,
                const optional<ISMRMRD::MetaContainer>& meta
        ) override {
            auto message_id = GADGET_MESSAGE_ISMRMRD_IMAGE;
            IO::write(stream, message_id);
            IO::write(stream,Image<T>{header,data,meta});
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

