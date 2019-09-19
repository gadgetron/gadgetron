#include "AcquisitionBucketReader.h"

#include "MessageID.h"
#include "io/primitives.h"

using namespace Gadgetron;
using namespace Gadgetron::Core;


namespace {

    struct bundle_meta {
        uint64_t count;
        struct {
            uint64_t header;
            uint64_t data;
            uint64_t trajectory;
        } nbytes;
    };

    struct stats_meta {
        uint64_t nbytes;
    };

    struct waveform_meta {
        uint64_t count;
        struct {
            uint64_t header;
            uint64_t data;
        } nbytes;
    };

    struct bucket_meta {
        bundle_meta   data, reference;
        stats_meta    data_stats, reference_stats;
        waveform_meta waveforms;
    };

}

namespace Gadgetron::Core::Readers {

    Message AcquisitionBucketReader::read(std::istream& stream) {


        auto meta = Core::IO::read<bucket_meta>(stream);



        return Message();
    }

    uint16_t AcquisitionBucketReader::slot() {
        return GADGET_MESSAGE_BUCKET;
    }

    GADGETRON_READER_EXPORT(AcquisitionBucketReader)
}