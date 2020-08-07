#include "AcquisitionWriter.h"
#include "io/primitives.h"
#include "MessageID.h"

namespace Gadgetron::Core::Writers {

    void AcquisitionWriter::serialize(
            std::ostream &stream,
            const ISMRMRD::AcquisitionHeader& header,
            const Gadgetron::hoNDArray<std::complex<float>>& data,
            const Core::optional<Gadgetron::hoNDArray<float>>& trajectory
    ) {
        IO::write(stream, GADGET_MESSAGE_ISMRMRD_ACQUISITION);
        IO::write(stream, header);
        if (trajectory)
            IO::write(stream, trajectory->get_data_ptr(), trajectory->get_number_of_elements());
        IO::write(stream, data.get_data_ptr(), data.get_number_of_elements());
    }

    GADGETRON_WRITER_EXPORT(AcquisitionWriter);
}


