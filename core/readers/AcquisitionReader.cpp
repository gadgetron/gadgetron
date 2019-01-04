
#include "io/readers.h"
#include "mri_core_data.h"

#include "AcquisitionReader.h"

namespace Gadgetron::Core::Readers {


    std::unique_ptr<Core::Message> AcquisitionReader::read(std::istream &stream) {

        using namespace Core;
        using namespace std::literals;

        std::vector<std::unique_ptr<Message>> messages;

        auto header = ISMRMRD::AcquisitionHeader{};

        IO::read(stream, header);

        messages.emplace_back(std::make_unique<TypedMessage<ISMRMRD::AcquisitionHeader>>(header));

        if (header.trajectory_dimensions) {
            auto trajectory = std::make_unique<hoNDArray<float>>(header.trajectory_dimensions,
                                                                 header.number_of_samples);
            IO::read(stream, *trajectory);
            messages.emplace_back(std::make_unique<TypedMessage<hoNDArray<float>>>(std::move(trajectory)));
        }

        {
            auto data = std::make_unique<hoNDArray<std::complex<float>>>(header.number_of_samples,
                                                                         header.active_channels);
            IO::read(stream, *data);
            messages.emplace_back(std::make_unique<TypedMessage<hoNDArray<std::complex<float>>>>(std::move(data)));
        }

        return std::make_unique<Core::MessageTuple>(std::move(messages));
    }

    uint16_t AcquisitionReader::slot() {
        return 1008;
    }

    GADGETRON_READER_EXPORT(AcquisitionReader)
}


