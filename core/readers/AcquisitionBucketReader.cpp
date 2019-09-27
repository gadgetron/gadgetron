#include "AcquisitionBucketReader.h"
#include <mri_core_acquisition_bucket.h>

#include "MessageID.h"
#include "io/primitives.h"
#include "io/adapt_struct.h"

using namespace Gadgetron;
using namespace Gadgetron::Core;


GADGETRON_ADAPT_STRUCT(AcquisitionBucketStats,
    GADGETRON_ACCESS_ELEMENT(kspace_encode_step_1),
    GADGETRON_ACCESS_ELEMENT(kspace_encode_step_2),
    GADGETRON_ACCESS_ELEMENT(slice),
    GADGETRON_ACCESS_ELEMENT(phase),
    GADGETRON_ACCESS_ELEMENT(contrast),
    GADGETRON_ACCESS_ELEMENT(repetition),
    GADGETRON_ACCESS_ELEMENT(set),
    GADGETRON_ACCESS_ELEMENT(segment),
    GADGETRON_ACCESS_ELEMENT(average)
)
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
        return acqs;

    }

    std::vector<Core::Waveform> read_waveforms(std::istream& stream, waveform_meta sizes){
        auto wavs = std::vector<Core::Waveform >(sizes.count);
        for(auto& wav : wavs){
            auto& head = std::get<0>(wav);
            head = IO::read<ISMRMRD::WaveformHeader>(stream);
            auto& data = std::get<1>(wav);
            data = hoNDArray<uint32_t>(head.number_of_samples,head.channels);
        }

        for (auto& wav : wavs){
            auto& data = std::get<1>(wav);
            stream.read(reinterpret_cast<char*>(data.data()),data.get_number_of_bytes());
        }
        return wavs;
    }

}

namespace Gadgetron::Core::Readers {

    Message AcquisitionBucketReader::read(std::istream& stream) {

        auto meta = Core::IO::read<bucket_meta>(stream);
        auto bucket = AcquisitionBucket{};
        bucket.data_ = read_acquisitions(stream,meta.data);
        bucket.datastats_ = IO::read<std::vector<AcquisitionBucketStats>>(stream);
        bucket.ref_ = read_acquisitions(stream,meta.reference);
        bucket.refstats_ = IO::read<std::vector<AcquisitionBucketStats>>(stream);
        bucket.waveform_ = read_waveforms(stream,meta.waveforms);


        return Message(std::move(bucket));
    }

    uint16_t AcquisitionBucketReader::slot() {
        return GADGET_MESSAGE_BUCKET;
    }

    GADGETRON_READER_EXPORT(AcquisitionBucketReader)
}