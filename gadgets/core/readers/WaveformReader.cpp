
#include "io/readers.h"
#include "mri_core_data.h"

#include "WaveformReader.h"

namespace Gadgetron::Core::Readers {

    std::unique_ptr<Core::Message> WaveformReader::read(std::istream& stream) {

        using namespace Core;
        using namespace std::literals;

        auto wave = std::make_unique<Waveform>();

        IO::read(stream,wave->header);

        wave->data = hoNDArray<uint32_t>(wave->header.number_of_samples,wave->header.channels);

        IO::read(stream,wave->data);

        return std::unique_ptr<Message>(new TypedMessage<Waveform>(std::move(wave)));
    }

    uint16_t WaveformReader::port() {
        return 1026;
    }

    GADGETRON_READER_EXPORT(WaveformReader)
}