#include "AcquisitionBucketReader.h"

#include "MessageID.h"

using namespace Gadgetron;
using namespace Gadgetron::Core;

namespace {


}

namespace Gadgetron::Core::Readers {

    Message AcquisitionBucketReader::read(std::istream& stream) {
        return Message();
    }

    uint16_t AcquisitionBucketReader::slot() {
        return GADGET_MESSAGE_BUCKET;
    }

    GADGETRON_READER_EXPORT(AcquisitionBucketReader)
}