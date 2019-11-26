#include "IsmrmrdImageArrayWriter.h"
#include "io/ismrmrd_types.h"
#include "io/adapt_struct.h"
#include "MessageID.h"

GADGETRON_ADAPT_STRUCT(Gadgetron::IsmrmrdImageArray,GADGETRON_ACCESS_ELEMENT(data_), GADGETRON_ACCESS_ELEMENT(headers_),GADGETRON_ACCESS_ELEMENT(meta_), GADGETRON_ACCESS_ELEMENT(waveform_), GADGETRON_ACCESS_ELEMENT(acq_headers_))

void Gadgetron::Core::Writers::IsmrmrdImageArrayWriter::serialize(
    std::ostream& stream, const Gadgetron::IsmrmrdImageArray& image_array) {
    IO::write(stream, MessageID::GADGET_MESSAGE_ISMRMRD_IMAGE_ARRAY);
    IO::write(stream,image_array);
}

namespace Gadgetron::Core::Writers{
    GADGETRON_WRITER_EXPORT(IsmrmrdImageArrayWriter)
}
