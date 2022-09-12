//
// Created by david on 2/11/2019.
//

#include "BufferWriter.h"
#include "Types.h"
#include "io/primitives.h"
#include <boost/hana/adapt_struct.hpp>
#include "MessageID.h"
#include "io/adapt_struct.h"


GADGETRON_ADAPT_STRUCT(Gadgetron::IsmrmrdReconBit, GADGETRON_ACCESS_ELEMENT(data_), GADGETRON_ACCESS_ELEMENT(ref_))

GADGETRON_ADAPT_STRUCT(Gadgetron::IsmrmrdDataBuffered, GADGETRON_ACCESS_ELEMENT(data_), GADGETRON_ACCESS_ELEMENT(trajectory_), GADGETRON_ACCESS_ELEMENT(density_),GADGETRON_ACCESS_ELEMENT(headers_),GADGETRON_ACCESS_ELEMENT(sampling_))
GADGETRON_ADAPT_STRUCT(Gadgetron::IsmrmrdReconData,GADGETRON_ACCESS_ELEMENT(rbit_))


void
Gadgetron::Core::Writers::BufferWriter::serialize(std::ostream &stream, const Gadgetron::IsmrmrdReconData &reconData) {
    static_assert(!Gadgetron::Core::is_trivially_copyable_v<IsmrmrdReconData>);
    GDEBUG("Sending out reconData\n");
    IO::write(stream,MessageID::GADGET_MESSAGE_RECONDATA);
    IO::write(stream, reconData);
}

namespace Gadgetron::Core::Writers {
    GADGETRON_WRITER_EXPORT(BufferWriter)
}