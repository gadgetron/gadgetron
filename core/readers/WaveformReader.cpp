
#include "io/primitives.h"
#include "mri_core_data.h"

#include "MessageID.h"
#include "WaveformReader.h"

#include <ismrmrd/waveform.h>

namespace Gadgetron::Core::Readers {

    Core::Message WaveformReader::read(std::istream& stream) {

        using namespace Core;
        using namespace std::literals;

        auto header = IO::read<ISMRMRD::WaveformHeader>(stream);
        auto data = hoNDArray<uint32_t>(header.number_of_samples, header.channels);

        IO::read(stream, data.data(), data.size());

        return Message(header, std::move(data));
    }

    uint16_t WaveformReader::slot() {
        return GADGET_MESSAGE_ISMRMRD_WAVEFORM;
    }

    GADGETRON_READER_EXPORT(WaveformReader)
}