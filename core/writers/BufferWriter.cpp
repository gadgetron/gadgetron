//
// Created by david on 2/11/2019.
//

#include "BufferWriter.h"
#include "Types.h"
#include "io/primitives.h"
#include <boost/hana/adapt_struct.hpp>
#include "MessageID.h"

BOOST_HANA_ADAPT_STRUCT(Gadgetron::IsmrmrdReconBit, data_, ref_);
BOOST_HANA_ADAPT_STRUCT(Gadgetron::IsmrmrdDataBuffered, data_, trajectory_, density_, headers_, sampling_);
BOOST_HANA_ADAPT_STRUCT(Gadgetron::IsmrmrdReconData, rbit_);


void
Gadgetron::Core::Writers::BufferWriter::serialize(std::ostream &stream, const Gadgetron::IsmrmrdReconData &reconData) {

    static_assert(!std::is_trivially_copyable_v<IsmrmrdReconData>);
    IO::write(stream,MessageID::GADGET_MESSAGE_ISMRMRD_BUFFER);
    IO::write(stream, reconData);

}
namespace Gadgetron::Core::Writers {
    GADGETRON_WRITER_EXPORT(BufferWriter)
}