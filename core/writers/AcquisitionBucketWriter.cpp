#include "AcquisitionBucketWriter.h"

#include "io/primitives.h"
#include "io/adapt_struct.h"

#include "MessageID.h"

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
        bundle_meta   data, reference;
        stats_meta    data_stats, reference_stats;
        waveform_meta waveforms;
    };

    struct acquisition_streams {
        std::stringstream header, data, trajectory, stats;
    };

    struct waveform_streams {
        std::stringstream header, data;
    };

    template<class T>
    void write_array_data(std::ostream &stream, const hoNDArray<T> &array) {
        IO::write(stream, array.data(), array.get_number_of_elements());
    }

    void serialize_acquisitions(
            acquisition_streams &streams,
            const std::vector<Core::Acquisition> &acquisitions,
            const std::vector<AcquisitionBucketStats> &stats
    ) {
        for (const auto& acq : acquisitions) {
            IO::write(streams.header, std::get<0>(acq));
            write_array_data(streams.data, std::get<1>(acq));

            if (std::get<2>(acq)) write_array_data(streams.trajectory, *std::get<2>(acq));
        }

        IO::write(streams.stats,stats);
    }

    void serialize_waveforms(
            waveform_streams &streams,
            const std::vector<Core::Waveform> &waveforms
    ) {
        for (const auto& waveform : waveforms) {
            IO::write(streams.header, std::get<0>(waveform));
            IO::write(streams.data, std::get<1>(waveform).data(), std::get<1>(waveform).size());
        }
    }

    bucket_meta initialize_bucket_meta(const Gadgetron::AcquisitionBucket&bucket) {
        bucket_meta meta{};
        meta.data.count = bucket.data_.size();
        meta.reference.count = bucket.ref_.size();
        meta.waveforms.count = bucket.waveform_.size();
        return meta;
    }

    void update_meta(
            bundle_meta &meta,
            stats_meta &stats,
            acquisition_streams &streams
    ) {
        meta.nbytes.header     = streams.header.tellp();
        meta.nbytes.data       = streams.data.tellp();
        meta.nbytes.trajectory = streams.trajectory.tellp();
        stats.nbytes           = streams.stats.tellp();
    }

    void update_meta(
            waveform_meta &meta,
            waveform_streams &streams
    ) {
        meta.nbytes.header = streams.header.tellp();
        meta.nbytes.data   = streams.data.tellp();
    }

    void write_stringstream(std::ostream &output, std::stringstream &strstream) {
        auto str = strstream.str();
        output.write(str.data(), str.size());
    }

    void write_streams(std::ostream &output, acquisition_streams &streams) {
        write_stringstream(output, streams.header);
        write_stringstream(output, streams.trajectory);
        write_stringstream(output, streams.data);
        write_stringstream(output, streams.stats);
    }

    void write_streams(std::ostream &output, waveform_streams &streams) {
        write_stringstream(output, streams.header);
        write_stringstream(output, streams.data);
    }
}

namespace Gadgetron::Core::Writers {

    void AcquisitionBucketWriter::serialize(
            std::ostream &stream,
            const Gadgetron::AcquisitionBucket&bucket
    ) {
        struct {
            acquisition_streams data, reference;
            waveform_streams waveforms;
        } streams;

        bucket_meta meta = initialize_bucket_meta(bucket);

        serialize_acquisitions(streams.data, bucket.data_, bucket.datastats_);
        serialize_acquisitions(streams.reference, bucket.ref_, bucket.refstats_);
        serialize_waveforms(streams.waveforms, bucket.waveform_);

        update_meta(meta.data, meta.data_stats, streams.data);
        update_meta(meta.reference, meta.reference_stats, streams.reference);
        update_meta(meta.waveforms, streams.waveforms);

        IO::write(stream, MessageID::GADGET_MESSAGE_BUCKET);
        IO::write(stream, meta);

        write_streams(stream, streams.data);
            write_streams(stream, streams.reference);
        write_streams(stream, streams.waveforms);
    }

    GADGETRON_WRITER_EXPORT(AcquisitionBucketWriter)
}



