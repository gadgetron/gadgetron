//
// Created by dchansen on 2/25/19.
//

#include "IsmrmrdImageArrayWriter.h"
#include "io/ismrmrd_types.h"
#include "MessageID.h"

BOOST_HANA_ADAPT_STRUCT(Gadgetron::IsmrmrdImageArray, data_, headers_, meta_, waveform_, acq_headers_);



void Gadgetron::Core::Writers::IsmrmrdImageArrayWriter::serialize(
    std::ostream& stream, const Gadgetron::IsmrmrdImageArray& image_array) {
    IO::write(stream, MessageID::GADGET_MESSAGE_ISMRMRD_IMAGE_ARRAY);
    IO::write(stream,image_array);

}

namespace Gadgetron::Core::Writers{
    GADGETRON_WRITER_EXPORT(IsmrmrdImageArrayWriter)
}