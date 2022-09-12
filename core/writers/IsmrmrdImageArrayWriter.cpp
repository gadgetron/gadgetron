#include "IsmrmrdImageArrayWriter.h"
#include "io/ismrmrd_types.h"
#include "io/adapt_struct.h"
#include "MessageID.h"

void Gadgetron::Core::Writers::IsmrmrdImageArrayWriter::serialize(
    std::ostream& stream, const Gadgetron::IsmrmrdImageArray& image_array) {
    IO::write(stream, MessageID::GADGET_MESSAGE_ISMRMRD_IMAGE_ARRAY);
    IO::write(stream,image_array);
}

namespace Gadgetron::Core::Writers{
    GADGETRON_WRITER_EXPORT(IsmrmrdImageArrayWriter)
}
