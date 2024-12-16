#include "BucketToBufferGadget.h"
#include "hoNDArray_elemwise.h"
#include "hoNDArray_reductions.h"
#include <boost/algorithm/string.hpp>


using BufferKey =  Gadgetron::BucketToBufferGadget::BufferKey;


namespace std {
    template<>
    struct less<BufferKey>{
        bool operator()(const BufferKey& idx1, const BufferKey& idx2) const {
            return std::tie(idx1.average,idx1.slice,idx1.contrast,idx1.phase,idx1.repetition,idx1.set,idx1.segment) <
                std::tie(idx2.average,idx2.slice,idx2.contrast,idx2.phase,idx2.repetition,idx2.set,idx2.segment);
        }
    };

    template<> struct equal_to<BufferKey>{
        bool operator()(const BufferKey& idx1, const BufferKey& idx2) const {
            return idx1.average == idx2.average
                   && idx1.slice == idx2.slice && idx1.contrast == idx2.contrast && idx1.phase == idx2.phase
                   && idx1.repetition == idx2.repetition && idx1.set == idx2.set && idx1.segment == idx2.segment;
        }
    };
}
namespace Gadgetron {
    namespace {

        mrd::ReconAssembly& getReconAssembly(std::map<BufferKey, mrd::ReconData>& recon_data_buffers,
            const BufferKey& key, uint32_t espace) {

            // Look up the mrd::ReconBuffer entry corresponding to this encoding space
            // create if needed and set the fields of view and matrix size
            if (recon_data_buffers[key].buffers.size() < (espace + 1)) {
                recon_data_buffers[key].buffers.resize(espace + 1);
            }

            return recon_data_buffers[key].buffers[espace];
        }

        uint32_t getLimitSize(std::optional<mrd::LimitType> const& limit) {
            if (limit) {
                return limit->maximum - limit->minimum + 1;
            }
            return 1;
        }
    }

    void BucketToBufferGadget::process(Core::InputChannel<mrd::AcquisitionBucket>& input, Core::OutputChannel& out) {

        for (auto acq_bucket : input) {
            std::map<BufferKey, mrd::ReconData> recon_data_buffers;
            // GDEBUG_STREAM("BUCKET_SIZE " << acq_bucket.data.size() << " ESPACE " << acq_bucket.refstats.size());

            // Buffer the reference data
            for (auto& acq : acq_bucket.ref) {
                auto key              = getKey(acq.head.idx);
                uint32_t espace       = acq.head.encoding_space_ref.value_or(0);
                mrd::ReconAssembly& assembly = getReconAssembly(recon_data_buffers, key, espace);
                if (!assembly.ref) {
                    assembly.ref = makeDataBuffer(acq, header.encoding[espace], acq_bucket.refstats[espace], true);
                }

                add_acquisition(*assembly.ref, acq, header.encoding[espace], acq_bucket.refstats[espace], true);
            }

            // Buffer the bucketed Acquisitions
            for (auto& acq : acq_bucket.data) {
                auto key              = getKey(acq.head.idx);
                uint32_t espace       = acq.head.encoding_space_ref.value_or(0);
                mrd::ReconAssembly& assembly = getReconAssembly(recon_data_buffers, key, espace);
                if (assembly.data.data.empty()) {
                    assembly.data = makeDataBuffer(acq, header.encoding[espace], acq_bucket.datastats[espace], false);
                }

                add_acquisition(assembly.data, acq, header.encoding[espace], acq_bucket.datastats[espace], false);
            }

            // Send all the ReconData messages
            GDEBUG("End of bucket reached, sending out %d ReconData buffers\n", recon_data_buffers.size());

            for (auto& recon_data_buffer : recon_data_buffers) {
                if (acq_bucket.waveforms.empty())
                {
                    // GDEBUG_STREAM("Sending out ReconData buffers without waveforms ...");
                    out.push(recon_data_buffer.second);
                }
                else
                {
                    // GDEBUG_STREAM("Sending out ReconData buffers with waveforms ...");
                    out.push(recon_data_buffer.second, acq_bucket.waveforms);
                }
            }
        }
    }

    namespace {
        void clear(BucketToBufferGadget::Dimension dim, BufferKey& idx) {
            switch (dim) {

            case BucketToBufferGadget::Dimension::average: idx.average = 0; break;
            case BucketToBufferGadget::Dimension::contrast: idx.contrast = 0; break;
            case BucketToBufferGadget::Dimension::phase: idx.phase = 0; break;
            case BucketToBufferGadget::Dimension::repetition: idx.repetition = 0; break;
            case BucketToBufferGadget::Dimension::set: idx.set = 0; break;
            case BucketToBufferGadget::Dimension::segment: idx.segment = 0; break;
            case BucketToBufferGadget::Dimension::slice: break;
            case BucketToBufferGadget::Dimension::none: break;
            default: throw std::runtime_error("Invalid enum encountered");
            }
        }

        size_t getDimensionKey(BucketToBufferGadget::Dimension dim, const mrd::EncodingCounters& idx) {
            switch (dim) {

            case BucketToBufferGadget::Dimension::average: return idx.average.value_or(0);
            case BucketToBufferGadget::Dimension::contrast: return idx.contrast.value_or(0);
            case BucketToBufferGadget::Dimension::phase: return idx.phase.value_or(0);
            case BucketToBufferGadget::Dimension::repetition: return idx.repetition.value_or(0);
            case BucketToBufferGadget::Dimension::set: return idx.set.value_or(0);
            case BucketToBufferGadget::Dimension::segment: return idx.segment.value_or(0);
            case BucketToBufferGadget::Dimension::slice: return 0;
            case BucketToBufferGadget::Dimension::none: return 0;
            default: throw std::runtime_error("Invalid enum encountered");
            }
        }
    }

    BufferKey BucketToBufferGadget::getKey(const mrd::EncodingCounters& idx) const {
        BufferKey key(idx);
        clear(N_dimension, key);
        clear(S_dimension, key);
        if (!split_slices)
            key.slice = 0;
        if (ignore_segment)
            key.segment = 0;
        return key;
    }

    namespace {
        uint32_t getSizeFromDimension(BucketToBufferGadget::Dimension dimension, const mrd::EncodingLimitsType& stats) {
            switch (dimension) {
            case BucketToBufferGadget::Dimension::phase:
                if (stats.phase) {
                    return stats.phase->maximum - stats.phase->minimum + 1;
                }
            case BucketToBufferGadget::Dimension::contrast:
                if (stats.contrast) {
                    return stats.contrast->maximum - stats.contrast->minimum + 1;
                }
            case BucketToBufferGadget::Dimension::repetition:
                if (stats.repetition) {
                    return stats.repetition->maximum - stats.repetition->minimum + 1;
                }
            case BucketToBufferGadget::Dimension::set:
                if (stats.set) {
                    return stats.set->maximum - stats.set->minimum + 1;
                }
            case BucketToBufferGadget::Dimension::segment: // TODO: Is this an intentional fallthrough? See 3b643b814b3b9a88c7dbc48879471ce718c7ab56
            case BucketToBufferGadget::Dimension::average:
                if (stats.average) {
                    return stats.average->maximum - stats.average->minimum + 1;
                }
            case BucketToBufferGadget::Dimension::slice:
                if (stats.slice) {
                    return stats.slice->maximum - stats.slice->minimum + 1;
                }
            case BucketToBufferGadget::Dimension::none:;
                return 1;
            default: throw std::runtime_error("Illegal enum value.");
            }

            return 1;
        }
    }

    mrd::ReconBuffer BucketToBufferGadget::makeDataBuffer(const mrd::Acquisition& acq,
        mrd::EncodingType encoding, const mrd::EncodingLimitsType& stats, bool forref) const
    {
        // Allocate the reference data array
        // 7D,  fixed order [E0, E1, E2, CHA, N, S, LOC]
        // 11D, fixed order [E0, E1, E2, CHA, SLC, PHS, CON, REP, SET, SEG, AVE]
        const uint32_t NE0 = getNE0(acq, encoding);
        uint32_t NE1 = getNE1(encoding, stats, forref);
        uint32_t NE2 = getNE2(encoding, stats, forref);
        size_t NCHA = acq.Coils();
        uint32_t NLOC = getNLOC(encoding, stats);
        uint32_t NN = getSizeFromDimension(N_dimension, stats);
        uint32_t NS = getSizeFromDimension(S_dimension, stats);

        GDEBUG_CONDITION_STREAM(verbose, "Data dimensions [RO E1 E2 CHA N S SLC] : ["
                                             << NE0 << " " << NE1 << " " << NE2 << " " << NCHA << " " << NN << " " << NS
                                             << " " << NLOC << "]");

        mrd::ReconBuffer buffer;

        // Allocate the array for the data
        buffer.data = hoNDArray<std::complex<float>>(NE0, NE1, NE2, NCHA, NN, NS, NLOC);
        clear(&buffer.data);

        // Allocate the array for the headers
        buffer.headers = hoNDArray<mrd::AcquisitionHeader>(NE1, NE2, NN, NS, NLOC);

        // Allocate the array for the trajectories
        if (acq.TrajectoryDimensions() > 0 && acq.TrajectorySamples() > 0) {
            auto basis = acq.TrajectoryDimensions();
            auto samples = acq.TrajectorySamples();
            buffer.trajectory = hoNDArray<float>(samples, basis, NE1, NE2, NN, NS, NLOC);
            clear(&buffer.trajectory);
        }

        // Add the sampling description
        buffer.sampling = createSamplingDescription(encoding, stats, acq, forref);

        return buffer;
    }

    uint32_t BucketToBufferGadget::getNLOC(
        const mrd::EncodingType& encoding, const mrd::EncodingLimitsType& stats) const {
        uint32_t NLOC;
        if (split_slices) {
            NLOC = 1;
        } else {
            if (encoding.encoding_limits.slice.has_value()) {
                NLOC = encoding.encoding_limits.slice->maximum - encoding.encoding_limits.slice->minimum + 1;
            } else {
                NLOC = 1;
            }

            // if the AcquisitionAccumulateTriggerGadget sort by SLC, then the stats should be used to determine NLOC
            // size_t NLOC_received = stats.slice ? stats.slice->maximum - stats.slice->minimum + 1 : 1;
            size_t NLOC_received = getLimitSize(stats.slice); // TODO: Clean up
            if (NLOC_received < NLOC) {
                NLOC = NLOC_received;
            }
        }
        return NLOC;
    }
    uint32_t BucketToBufferGadget::getNE2(
        const mrd::EncodingType& encoding, const mrd::EncodingLimitsType& stats, bool forref) const {
        uint32_t NE2;

        /** TODO: This is ugly... */

        if (encoding.trajectory == mrd::Trajectory::kCartesian || encoding.trajectory == mrd::Trajectory::kEpi) {
            if (encoding.parallel_imaging) {
                if (forref && encoding.parallel_imaging->calibration_mode &&
                        (encoding.parallel_imaging->calibration_mode.value() == mrd::CalibrationMode::kSeparate ||
                        encoding.parallel_imaging->calibration_mode.value() == mrd::CalibrationMode::kExternal)) {
                    NE2 = encoding.encoding_limits.kspace_encoding_step_2->maximum
                          - encoding.encoding_limits.kspace_encoding_step_2->minimum + 1;
                } else {
                    NE2 = encoding.encoded_space.matrix_size.z;
                }
            } else {
                if (encoding.encoding_limits.kspace_encoding_step_2.has_value()) {
                    NE2 = encoding.encoding_limits.kspace_encoding_step_2->maximum
                          - encoding.encoding_limits.kspace_encoding_step_2->minimum + 1;
                } else {
                    NE2 = encoding.encoded_space.matrix_size.z;
                }
            }
        } else {
            if (encoding.encoding_limits.kspace_encoding_step_2.has_value()) {
                NE2 = encoding.encoding_limits.kspace_encoding_step_2->maximum
                      - encoding.encoding_limits.kspace_encoding_step_2->minimum + 1;
            } else {
                NE2 = getLimitSize(stats.kspace_encoding_step_2);
            }
        }
        return NE2;
    }
    uint32_t BucketToBufferGadget::getNE1(
        const mrd::EncodingType& encoding, const mrd::EncodingLimitsType& stats, bool forref) const {
        uint32_t NE1;

        /** TODO: This is also ugly... */

        if (encoding.trajectory == mrd::Trajectory::kCartesian || encoding.trajectory == mrd::Trajectory::kEpi) {
            if (encoding.parallel_imaging) {
                if (forref && encoding.parallel_imaging->calibration_mode &&
                        (encoding.parallel_imaging->calibration_mode.value() == mrd::CalibrationMode::kSeparate ||
                        encoding.parallel_imaging->calibration_mode.value() == mrd::CalibrationMode::kExternal)) {
                    NE1 = stats.kspace_encoding_step_1 ? stats.kspace_encoding_step_1->maximum - stats.kspace_encoding_step_1->minimum + 1 : 1;
                } else {
                    NE1 = encoding.encoded_space.matrix_size.y;
                }
            } else {
                if (encoding.encoding_limits.kspace_encoding_step_1.has_value()) {
                    NE1 = encoding.encoding_limits.kspace_encoding_step_1->maximum
                          - encoding.encoding_limits.kspace_encoding_step_1->minimum + 1;
                } else {
                    NE1 = encoding.encoded_space.matrix_size.y;
                }
            }
        } else {
            if (encoding.encoding_limits.kspace_encoding_step_1.has_value()) {
                NE1 = encoding.encoding_limits.kspace_encoding_step_1->maximum
                      - encoding.encoding_limits.kspace_encoding_step_1->minimum + 1;
            } else {
                NE1 = getLimitSize(stats.kspace_encoding_step_1);
            }
        }
        return NE1;
    }
    uint32_t BucketToBufferGadget::getNE0(
        const mrd::Acquisition& acq, const mrd::EncodingType& encoding) const {
        uint32_t NE0;
        if (encoding.trajectory == mrd::Trajectory::kCartesian || encoding.trajectory == mrd::Trajectory::kEpi) {
            // if separate or external calibration mode, using the acq length for NE0
            if (encoding.parallel_imaging) {
                NE0 = acq.Samples();
            } else {
                NE0 = acq.Samples() - acq.head.discard_pre.value_or(0) - acq.head.discard_post.value_or(0);
            }
        } else {
            NE0 = acq.Samples() - acq.head.discard_pre.value_or(0) - acq.head.discard_post.value_or(0);
        }
        return NE0;
    }

    mrd::SamplingDescription BucketToBufferGadget::createSamplingDescription(const mrd::EncodingType& encoding,
        const mrd::EncodingLimitsType& stats, const mrd::Acquisition& acq, bool forref) const
    {
        mrd::SamplingDescription sampling;
        sampling.encoded_fov    = encoding.encoded_space.field_of_view_mm;
        sampling.encoded_matrix = encoding.encoded_space.matrix_size;
        sampling.recon_fov      = encoding.recon_space.field_of_view_mm;
        sampling.recon_matrix   = encoding.recon_space.matrix_size;

        sampling.sampling_limits.kspace_encoding_step_0.minimum = 0;
        sampling.sampling_limits.kspace_encoding_step_0.maximum = acq.Samples() - 1;
        sampling.sampling_limits.kspace_encoding_step_0.center = acq.Samples() / 2;

        sampling.sampling_limits.kspace_encoding_step_1.minimum = encoding.encoding_limits.kspace_encoding_step_1->minimum;
        sampling.sampling_limits.kspace_encoding_step_1.maximum = encoding.encoding_limits.kspace_encoding_step_1->maximum;
        sampling.sampling_limits.kspace_encoding_step_1.center = encoding.encoding_limits.kspace_encoding_step_1->center;

        sampling.sampling_limits.kspace_encoding_step_2.minimum = encoding.encoding_limits.kspace_encoding_step_2->minimum;
        sampling.sampling_limits.kspace_encoding_step_2.maximum = encoding.encoding_limits.kspace_encoding_step_2->maximum;
        sampling.sampling_limits.kspace_encoding_step_2.center = encoding.encoding_limits.kspace_encoding_step_2->center;

        if (verbose) {
            GDEBUG_STREAM("Encoding space : " << acq.head.encoding_space_ref.value_or(0) << " - "
                          << int(encoding.trajectory) << " - FOV : [ " << encoding.encoded_space.field_of_view_mm.x << " "
                          << encoding.encoded_space.field_of_view_mm.y << " " << encoding.encoded_space.field_of_view_mm.z
                          << " ] "
                          << " - Matris size : [ " << encoding.encoded_space.matrix_size.x << " "
                          << encoding.encoded_space.matrix_size.y << " " << encoding.encoded_space.matrix_size.z << " ] ");

            GDEBUG_STREAM("Sampling limits : "
                          << "- RO : [ " << sampling.sampling_limits.kspace_encoding_step_0.minimum << " "
                          << sampling.sampling_limits.kspace_encoding_step_0.center << " " << sampling.sampling_limits.kspace_encoding_step_0.maximum
                          << " ] - E1 : [ " << sampling.sampling_limits.kspace_encoding_step_1.minimum << " "
                          << sampling.sampling_limits.kspace_encoding_step_1.center << " " << sampling.sampling_limits.kspace_encoding_step_1.maximum
                          << " ] - E2 : [ " << sampling.sampling_limits.kspace_encoding_step_2.minimum << " "
                          << sampling.sampling_limits.kspace_encoding_step_2.center << " " << sampling.sampling_limits.kspace_encoding_step_2.maximum << " ]");
        }

        // For cartesian trajectories, assume that any oversampling has been removed.
        if (encoding.trajectory == mrd::Trajectory::kCartesian) {
            sampling.encoded_fov.x    = encoding.recon_space.field_of_view_mm.x;
            sampling.encoded_matrix.x = encoding.recon_space.matrix_size.x;
        } else {
            sampling.encoded_fov.x    = encoding.encoded_space.field_of_view_mm.x;
            sampling.encoded_matrix.x = encoding.encoded_space.matrix_size.x;
        }

        // For cartesian trajectories, assume that any oversampling has been removed.
        if (((encoding.trajectory == mrd::Trajectory::kCartesian)) || (encoding.trajectory == mrd::Trajectory::kEpi)) {
            sampling.sampling_limits.kspace_encoding_step_0.minimum = acq.head.discard_pre.value_or(0);
            sampling.sampling_limits.kspace_encoding_step_0.maximum = acq.Samples() - acq.head.discard_post.value_or(0) - 1;
            sampling.sampling_limits.kspace_encoding_step_0.center  = acq.Samples() / 2;
        } else {
            sampling.sampling_limits.kspace_encoding_step_0.minimum = 0;
            sampling.sampling_limits.kspace_encoding_step_0.maximum = encoding.encoded_space.matrix_size.x - 1;
            sampling.sampling_limits.kspace_encoding_step_0.center  = encoding.encoded_space.matrix_size.x / 2;
        }

        // if the scan is cartesian
        if (((encoding.trajectory == mrd::Trajectory::kCartesian) &&
                (!forref ||
                    (forref && encoding.parallel_imaging && encoding.parallel_imaging->calibration_mode &&
                    (encoding.parallel_imaging->calibration_mode.value() == mrd::CalibrationMode::kEmbedded))))
            || ((encoding.trajectory == mrd::Trajectory::kEpi) && !forref)) {

            int32_t space_matrix_offset_E1 = 0;
            if (encoding.encoding_limits.kspace_encoding_step_1.has_value()) {
                space_matrix_offset_E1 = (int32_t)encoding.encoded_space.matrix_size.y / 2
                                         - (int32_t)encoding.encoding_limits.kspace_encoding_step_1->center;
            }

            int32_t space_matrix_offset_E2 = 0;
            if (encoding.encoding_limits.kspace_encoding_step_2.has_value() && encoding.encoded_space.matrix_size.z > 1) {
                space_matrix_offset_E2 = (int32_t)encoding.encoded_space.matrix_size.z / 2
                                         - (int32_t)encoding.encoding_limits.kspace_encoding_step_2->center;
            }

            // E1
            sampling.sampling_limits.kspace_encoding_step_1.minimum
                = encoding.encoding_limits.kspace_encoding_step_1->minimum + space_matrix_offset_E1;
            sampling.sampling_limits.kspace_encoding_step_1.maximum
                = encoding.encoding_limits.kspace_encoding_step_1->maximum + space_matrix_offset_E1;
            sampling.sampling_limits.kspace_encoding_step_1.center = sampling.encoded_matrix.y / 2;

            GADGET_CHECK_THROW(sampling.sampling_limits.kspace_encoding_step_1.minimum < encoding.encoded_space.matrix_size.y);
            GADGET_CHECK_THROW(sampling.sampling_limits.kspace_encoding_step_1.maximum >= sampling.sampling_limits.kspace_encoding_step_1.minimum);
            GADGET_CHECK_THROW(sampling.sampling_limits.kspace_encoding_step_1.center >= sampling.sampling_limits.kspace_encoding_step_1.minimum);
            GADGET_CHECK_THROW(sampling.sampling_limits.kspace_encoding_step_1.center <= sampling.sampling_limits.kspace_encoding_step_1.maximum);

            // E2
            sampling.sampling_limits.kspace_encoding_step_2.minimum
                = encoding.encoding_limits.kspace_encoding_step_2->minimum + space_matrix_offset_E2;
            sampling.sampling_limits.kspace_encoding_step_2.maximum
                = encoding.encoding_limits.kspace_encoding_step_2->maximum + space_matrix_offset_E2;
            sampling.sampling_limits.kspace_encoding_step_2.center = sampling.encoded_matrix.z / 2;

            GADGET_CHECK_THROW(sampling.sampling_limits.kspace_encoding_step_2.minimum < encoding.encoded_space.matrix_size.y);
            GADGET_CHECK_THROW(sampling.sampling_limits.kspace_encoding_step_2.maximum >= sampling.sampling_limits.kspace_encoding_step_2.minimum);
            GADGET_CHECK_THROW(sampling.sampling_limits.kspace_encoding_step_2.center >= sampling.sampling_limits.kspace_encoding_step_2.minimum);
            GADGET_CHECK_THROW(sampling.sampling_limits.kspace_encoding_step_2.center <= sampling.sampling_limits.kspace_encoding_step_2.maximum);
        } else {
            sampling.sampling_limits.kspace_encoding_step_1.minimum    = encoding.encoding_limits.kspace_encoding_step_1->minimum;
            sampling.sampling_limits.kspace_encoding_step_1.maximum    = encoding.encoding_limits.kspace_encoding_step_1->maximum;
            sampling.sampling_limits.kspace_encoding_step_1.center = encoding.encoding_limits.kspace_encoding_step_1->center;

            sampling.sampling_limits.kspace_encoding_step_2.minimum    = encoding.encoding_limits.kspace_encoding_step_2->minimum;
            sampling.sampling_limits.kspace_encoding_step_2.maximum    = encoding.encoding_limits.kspace_encoding_step_2->maximum;
            sampling.sampling_limits.kspace_encoding_step_2.center = encoding.encoding_limits.kspace_encoding_step_2->center;
        }

        if (verbose) {
            GDEBUG_STREAM("Encoding space : "
                          << int(encoding.trajectory) << " - FOV : [ " << encoding.encoded_space.field_of_view_mm.x << " "
                          << encoding.encoded_space.field_of_view_mm.y << " " << encoding.encoded_space.field_of_view_mm.z
                          << " ] "
                          << " - Matrix size : [ " << encoding.encoded_space.matrix_size.x << " "
                          << encoding.encoded_space.matrix_size.y << " " << encoding.encoded_space.matrix_size.z << " ] ");

            GDEBUG_STREAM("Sampling limits : "
                          << "- kspace_encoding_step_0 : [ " << sampling.sampling_limits.kspace_encoding_step_0.minimum << " "
                          << sampling.sampling_limits.kspace_encoding_step_0.center << " " << sampling.sampling_limits.kspace_encoding_step_0.maximum
                          << " ] - kspace_encoding_step_1 : [ " << sampling.sampling_limits.kspace_encoding_step_1.minimum << " "
                          << sampling.sampling_limits.kspace_encoding_step_1.center << " " << sampling.sampling_limits.kspace_encoding_step_1.maximum
                          << " ] - kspace_encoding_step_2 : [ " << sampling.sampling_limits.kspace_encoding_step_2.minimum << " "
                          << sampling.sampling_limits.kspace_encoding_step_2.center << " " << sampling.sampling_limits.kspace_encoding_step_2.maximum << " ]");
        }
        return sampling;
    }

    void BucketToBufferGadget::add_acquisition(mrd::ReconBuffer& dataBuffer, const mrd::Acquisition& acq,
        mrd::EncodingType encoding, const mrd::EncodingLimitsType& stats, bool forref) {

        uint16_t NE0  = (uint16_t)dataBuffer.data.get_size(0);
        uint16_t NE1  = (uint16_t)dataBuffer.data.get_size(1);
        uint16_t NE2  = (uint16_t)dataBuffer.data.get_size(2);
        uint16_t NCHA = (uint16_t)dataBuffer.data.get_size(3);
        uint16_t NN   = (uint16_t)dataBuffer.data.get_size(4);
        uint16_t NS   = (uint16_t)dataBuffer.data.get_size(5);
        uint16_t NLOC = (uint16_t)dataBuffer.data.get_size(6);

        const size_t slice_loc = split_slices || NLOC == 1 ? 0 : acq.head.idx.slice.value_or(0);

        // Stuff the data
        uint32_t npts_to_copy = acq.Samples() - acq.head.discard_pre.value_or(0) - acq.head.discard_post.value_or(0);
        long long offset;
        if (encoding.trajectory == mrd::Trajectory::kCartesian || encoding.trajectory == mrd::Trajectory::kEpi) {
            if ((acq.Samples() == NE0)
                && ((acq.head.center_sample == acq.Samples() / 2) || acq.head.center_sample >= acq.Samples()))
            {
                // acq has been corrected for center, e.g. by asymmetric handling
                offset = acq.head.discard_pre.value_or(0);
            } else {
                offset = (long long)dataBuffer.sampling.sampling_limits.kspace_encoding_step_0.center - (long long)acq.head.center_sample.value_or(0);
            }
        } else {
            // TODO what about EPI with asymmetric readouts?
            // TODO any other sort of trajectory?
            offset = 0;
        }

        long long roffset = (long long)NE0 - npts_to_copy - offset;

        if ((offset < 0) | (roffset < 0)) {
            throw std::runtime_error("Acquired reference data does not fit into the reference data buffer.\n");
        }

        uint32_t NUsed = (uint32_t)getDimensionKey(N_dimension, acq.head.idx);
        if (NUsed >= NN)
            NUsed = NN - 1;

        uint32_t SUsed = (uint32_t)getDimensionKey(S_dimension, acq.head.idx);
        if (SUsed >= NS)
            SUsed = NS - 1;

        int32_t e1 = (int32_t)acq.head.idx.kspace_encode_step_1.value_or(0);
        int32_t e2 = (int32_t)acq.head.idx.kspace_encode_step_2.value_or(0);

        bool is_cartesian_sampling = (encoding.trajectory == mrd::Trajectory::kCartesian);
        bool is_epi_sampling       = (encoding.trajectory == mrd::Trajectory::kEpi);
        if (is_cartesian_sampling || is_epi_sampling) {
            if (!forref || (forref &&
                    encoding.parallel_imaging &&
                    encoding.parallel_imaging->calibration_mode &&
                    (encoding.parallel_imaging->calibration_mode.value() == mrd::CalibrationMode::kEmbedded)))
            {
                // compute the center offset for E1 and E2
                int32_t space_matrix_offset_E1 = 0;
                if (encoding.encoding_limits.kspace_encoding_step_1.has_value()) {
                    space_matrix_offset_E1 = (int32_t)encoding.encoded_space.matrix_size.y / 2
                                             - (int32_t)encoding.encoding_limits.kspace_encoding_step_1->center;
                }

                int32_t space_matrix_offset_E2 = 0;
                if (encoding.encoding_limits.kspace_encoding_step_2.has_value()
                    && encoding.encoded_space.matrix_size.z > 1) {
                    space_matrix_offset_E2 = (int32_t)encoding.encoded_space.matrix_size.z / 2
                                             - (int32_t)encoding.encoding_limits.kspace_encoding_step_2->center;
                }

                // compute the used e1 and e2 indices and make sure they are in the valid range
                e1 = (int32_t)acq.head.idx.kspace_encode_step_1.value_or(0) + space_matrix_offset_E1;
                e2 = (int32_t)acq.head.idx.kspace_encode_step_2.value_or(0) + space_matrix_offset_E2;
            }

            // for external or separate mode, it is possible the starting numbers of ref lines are not zero, therefore
            // it is needed to subtract the staring ref line number because the ref array size is set up by the actual
            // number of lines acquired only assumption for external or separate ref line mode is that all ref lines are
            // numbered sequentially the acquisition order of ref line can be arbitrary
            if (forref &&
                    encoding.parallel_imaging && encoding.parallel_imaging->calibration_mode &&
                    ((encoding.parallel_imaging->calibration_mode.value() == mrd::CalibrationMode::kSeparate)
                        || (encoding.parallel_imaging->calibration_mode.value() == mrd::CalibrationMode::kExternal)))
            {
                if (stats.kspace_encoding_step_1 && stats.kspace_encoding_step_1->minimum > 0) {
                    e1 = acq.head.idx.kspace_encode_step_1.value_or(0) - stats.kspace_encoding_step_1->minimum;
                }

                if (stats.kspace_encoding_step_2 && stats.kspace_encoding_step_2->minimum > 0) {
                    e2 = acq.head.idx.kspace_encode_step_2.value_or(0) - stats.kspace_encoding_step_2->minimum;
                }
            }

            if (e1 < 0 || e1 >= (int32_t)NE1) {
                // if the incoming line is outside the encoding limits, something is wrong
                GADGET_CHECK_THROW(
                    acq.head.idx.kspace_encode_step_1.value_or(0) >= encoding.encoding_limits.kspace_encoding_step_1->minimum
                    && acq.head.idx.kspace_encode_step_1.value_or(0) <= encoding.encoding_limits.kspace_encoding_step_1->maximum);

                // if the incoming line is inside encoding limits but outside the encoded matrix, do not include the data
                GWARN_STREAM(
                    "incoming readout "
                    << acq.head.scan_counter.value_or(0)
                    << " is inside the encoding limits, but outside the encoded matrix for kspace_encode_step_1 : "
                    << e1 << " out of " << NE1);
                return;
            }

            if (e2 < 0 || e2 >= (int32_t)NE2) {
                GADGET_CHECK_THROW(
                    acq.head.idx.kspace_encode_step_2.value_or(0) >= encoding.encoding_limits.kspace_encoding_step_2->minimum
                    && acq.head.idx.kspace_encode_step_2.value_or(0) <= encoding.encoding_limits.kspace_encoding_step_2->maximum);

                GWARN_STREAM(
                    "incoming readout "
                    << acq.head.scan_counter.value_or(0)
                    << " is inside the encoding limits, but outside the encoded matrix for kspace_encode_step_2 : "
                    << e2 << " out of " << NE2);
                return;
            }
        }

        // Stuff the data
        std::complex<float>* pData = &dataBuffer.data(offset, e1, e2, 0, NUsed, SUsed, slice_loc);

        for (size_t cha = 0; cha < NCHA; cha++) {
            auto dataptr = pData + cha * NE0 * NE1 * NE2;
            auto fromptr = &acq.data(acq.head.discard_pre.value_or(0), cha);
            std::copy(fromptr, fromptr + npts_to_copy, dataptr);
        }

        // Stuff the header
        dataBuffer.headers(e1, e2, NUsed, SUsed, slice_loc) = acq.head;

        if (acq.TrajectoryDimensions() > 0 && acq.TrajectorySamples() > 0) {
            // Stuff the trajectory
            float* trajptr = &dataBuffer.trajectory(offset, 0, e1, e2, NUsed, SUsed, slice_loc);
            auto* fromptr  = &acq.trajectory(acq.head.discard_pre.value_or(0), 0);
            std::copy(fromptr, fromptr + npts_to_copy * acq.TrajectoryDimensions(), trajptr);
        }
    }


    BucketToBufferGadget::BucketToBufferGadget(const Core::Context& context, const Core::GadgetProperties& props)
        : ChannelGadget(context, props), header{ context.header } {}

    namespace {
        using Dimension = BucketToBufferGadget::Dimension;
        const std::map<std::string, BucketToBufferGadget::Dimension> dimension_from_name
            = { { "average", Dimension::average }, { "contrast", Dimension::contrast }, { "phase", Dimension::phase },
                { "repetition", Dimension::repetition }, { "set", Dimension::set }, { "segment", Dimension::segment },
                { "slice", Dimension::slice }, { "", Dimension::none }, { "none", Dimension::none }
              };
    }

    void from_string(const std::string& str, BucketToBufferGadget::Dimension& dim) {
        auto lower = str;
        boost::to_lower(lower);
        dim = dimension_from_name.at(lower);
    }

    GADGETRON_GADGET_EXPORT(BucketToBufferGadget)

} // namespace Gadgetron
