
#include "io/readers.h"
#include "mri_core_data.h"

#include "WaveformReader.h"

namespace Gadgetron::Core::Readers {

    std::unique_ptr<Core::Message> WaveformReader::read(std::istream& stream) {

        using namespace Core;
        using namespace std::literals;


        auto header = std::make_unique<ISMRMRD::WaveformHeader>();
        IO::read(stream, *header);

        auto data = std::make_unique<hoNDArray<uint32_t>>(header->number_of_samples, header->channels);

        IO::read(stream, *data);

        return std::make_unique<MessageTuple>(std::move(header), std::move(data));
    }

    uint16_t WaveformReader::slot() {
        return 1026;
    }

    GADGETRON_READER_EXPORT(WaveformReader)
}