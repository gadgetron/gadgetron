//
// Created by david on 2/12/2019.
//

#include "BufferReader.h"
#include "mri_core_data.h"
#include "io/primitives.h"
BOOST_HANA_ADAPT_STRUCT(Gadgetron::IsmrmrdReconBit, data_, ref_);
BOOST_HANA_ADAPT_STRUCT(Gadgetron::IsmrmrdDataBuffered, data_, trajectory_, density_, headers_, sampling_);
BOOST_HANA_ADAPT_STRUCT(Gadgetron::IsmrmrdReconData, rbit_);

Gadgetron::Core::Message Gadgetron::Core::Readers::BufferReader::read(std::istream &stream) {
    return Message(IO::read<IsmrmrdReconData>(stream));
}

namespace Gadgetron::Core::Readers{
    GADGETRON_READER_EXPORT(BufferReader)
}
