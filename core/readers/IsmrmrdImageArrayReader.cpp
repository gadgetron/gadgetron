#include "IsmrmrdImageArrayReader.h"

#include "BufferReader.h"

#include "io/ismrmrd_types.h"
#include "mri_core_data.h"
#include "io/adapt_struct.h"

GADGETRON_ADAPT_STRUCT(Gadgetron::IsmrmrdImageArray, GADGETRON_ACCESS_ELEMENT(data_), GADGETRON_ACCESS_ELEMENT(headers_), GADGETRON_ACCESS_ELEMENT(meta_), GADGETRON_ACCESS_ELEMENT(waveform_), GADGETRON_ACCESS_ELEMENT(acq_headers_))

Gadgetron::Core::Message Gadgetron::Core::Readers::IsmrmrdImageArrayReader::read(std::istream& stream) {
    return Message(IO::read<IsmrmrdImageArray>(stream));
}
uint16_t Gadgetron::Core::Readers::IsmrmrdImageArrayReader::slot() {
    return MessageID::GADGET_MESSAGE_ISMRMRD_IMAGE_ARRAY;
}

namespace Gadgetron::Core::Readers{
    GADGETRON_READER_EXPORT(IsmrmrdImageArrayReader)
}
