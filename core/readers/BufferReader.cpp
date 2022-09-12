#include "BufferReader.h"
#include "mri_core_data.h"
#include "io/primitives.h"
#include "io/adapt_struct.h"

GADGETRON_ADAPT_STRUCT(Gadgetron::IsmrmrdReconBit, GADGETRON_ACCESS_ELEMENT(data_), GADGETRON_ACCESS_ELEMENT(ref_))

GADGETRON_ADAPT_STRUCT(Gadgetron::IsmrmrdDataBuffered, GADGETRON_ACCESS_ELEMENT(data_), GADGETRON_ACCESS_ELEMENT(trajectory_), GADGETRON_ACCESS_ELEMENT(density_), GADGETRON_ACCESS_ELEMENT(headers_), GADGETRON_ACCESS_ELEMENT(sampling_))
GADGETRON_ADAPT_STRUCT(Gadgetron::IsmrmrdReconData, GADGETRON_ACCESS_ELEMENT(rbit_))


Gadgetron::Core::Message Gadgetron::Core::Readers::BufferReader::read(std::istream &stream) {
    return Message(IO::read<IsmrmrdReconData>(stream));
}

uint16_t Gadgetron::Core::Readers::BufferReader::slot() {
    return MessageID::GADGET_MESSAGE_RECONDATA;
}

namespace Gadgetron::Core::Readers{
    GADGETRON_READER_EXPORT(BufferReader)
}

