#include "WaveformWriter.h"
#include "io/primitives.h"
#include "MessageID.h"

namespace Gadgetron::Core::Writers {

    void WaveformWriter::serialize(
            std::ostream &stream,
            const ISMRMRD::WaveformHeader& header,
            const hoNDArray<uint32_t>& array
    ) {

        IO::write(stream, GADGET_MESSAGE_ISMRMRD_WAVEFORM);
        IO::write(stream, header);
        IO::write(stream, array.get_data_ptr(), array.get_number_of_elements());
    }

    GADGETRON_WRITER_EXPORT(WaveformWriter)
}

