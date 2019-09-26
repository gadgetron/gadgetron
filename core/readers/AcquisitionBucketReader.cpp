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
        bundle_meta data, reference;
        stats_meta data_stats, reference_stats;
        waveform_meta waveforms;
    };

    std::vector<Core::Acquisition> read_acquisitions(std::istream& stream, bundle_meta sizes) {
        auto acqs = std::vector<Core::Acquisition>(sizes.count);
        for (auto& acq : acqs) {
            auto& head = std::get<ISMRMRD::AcquisitionHeader>(acq);
            head       = IO::read<ISMRMRD::AcquisitionHeader>(stream);
            auto& data = std::get<hoNDArray<std::complex<float>>>(acq);
            data       = hoNDArray<std::complex<float>>(head.number_of_samples, head.active_channels);
            if (head.trajectory_dimensions > 0) {
                auto& traj = std::get<Core::optional<hoNDArray<float>>>(acq);
                traj       = hoNDArray<float>(head.trajectory_dimensions, head.number_of_samples);
            }
        }


        for (auto& acq : acqs) {
            auto& traj = std::get<Core::optional<hoNDArray<float>>>(acq);
            if (traj)
                stream.read(reinterpret_cast<char*>(traj->data()), traj->get_number_of_bytes());
        }
        for (auto& acq : acqs) {
            auto& data = std::get<hoNDArray<std::complex<float>>>(acq);
            stream.read(reinterpret_cast<char*>(data.data()), data.get_number_of_bytes());
        }

    }

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