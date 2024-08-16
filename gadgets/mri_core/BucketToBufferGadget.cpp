#include "BucketToBufferGadget.h"
#include "hoNDArray_elemwise.h"
#include "hoNDArray_reductions.h"
#include "mri_core_data.h"
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

        IsmrmrdReconBit& getRBit(std::map<BufferKey, IsmrmrdReconData>& recon_data_buffers,
            const BufferKey& key, uint16_t espace) {

            // Look up the DataBuffered entry corresponding to this encoding space
            // create if needed and set the fields of view and matrix size
            if (recon_data_buffers[key].rbit_.size() < (espace + 1)) {
                recon_data_buffers[key].rbit_.resize(espace + 1);
            }

            return recon_data_buffers[key].rbit_[espace];
        }

    }

    void BucketToBufferGadget::process(Core::InputChannel<AcquisitionBucket>& input, Core::OutputChannel& out) {

        for (auto acq_bucket : input) {
            std::map<BufferKey, IsmrmrdReconData> recon_data_buffers;
            GDEBUG_STREAM("BUCKET_SIZE " << acq_bucket.data_.size() << " ESPACE " << acq_bucket.refstats_.size());
            // Iterate over the reference data of the bucket
            for (auto& acq : acq_bucket.ref_) {
                // Get a reference to the header for this acquisition

                const auto& acqhdr    = std::get<ISMRMRD::AcquisitionHeader>(acq);
                auto key              = getKey(acqhdr.idx);
                uint16_t espace       = acqhdr.encoding_space_ref;
                IsmrmrdReconBit& rbit = getRBit(recon_data_buffers, key, espace);
                if (!rbit.ref_) {
                    rbit.ref_ = makeDataBuffer(acqhdr, header.encoding[espace], acq_bucket.refstats_[espace], true);
                    rbit.ref_->sampling_ = createSamplingDescription(
                        header.encoding[espace], acq_bucket.refstats_[espace], acqhdr, true);
                }

                add_acquisition(*rbit.ref_, acq, header.encoding[espace], acq_bucket.refstats_[espace], true);

                // Stuff the data, header and trajectory into this data buffer
            }

            // Iterate over the reference data of the bucket
            for (auto& acq : acq_bucket.data_) {
                // Get a reference to the header for this acquisition

                const auto& acqhdr    = std::get<ISMRMRD::AcquisitionHeader>(acq);
                auto key              = getKey(acqhdr.idx);
                uint16_t espace       = acqhdr.encoding_space_ref;
                IsmrmrdReconBit& rbit = getRBit(recon_data_buffers, key, espace);
                if (rbit.data_.data_.empty()) {
                    rbit.data_ = makeDataBuffer(acqhdr, header.encoding[espace], acq_bucket.datastats_[espace], false);
                    rbit.data_.sampling_ = createSamplingDescription(
                        header.encoding[espace], acq_bucket.datastats_[espace], acqhdr, false);
                }

                add_acquisition(rbit.data_, acq, header.encoding[espace], acq_bucket.datastats_[espace], false);

                // Stuff the data, header and trajectory into this data buffer
            }

            // Send all the ReconData messages
            GDEBUG("End of bucket reached, sending out %d ReconData buffers\n", recon_data_buffers.size());

            for (auto& recon_data_buffer : recon_data_buffers) {
                if (acq_bucket.waveform_.empty())
                {
                    GDEBUG_STREAM("Sending out ReconData buffers without waveforms ...");
                    out.push(recon_data_buffer.second);
                }
                else
                {
                    GDEBUG_STREAM("Sending out ReconData buffers with waveforms ...");
                    out.push(recon_data_buffer.second, acq_bucket.waveform_);
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

        size_t getDimensionKey(BucketToBufferGadget::Dimension dim, const ISMRMRD::EncodingCounters& idx) {
            switch (dim) {

            case BucketToBufferGadget::Dimension::average: return idx.average;
            case BucketToBufferGadget::Dimension::contrast: return idx.contrast;
            case BucketToBufferGadget::Dimension::phase: return idx.phase;
            case BucketToBufferGadget::Dimension::repetition: return idx.repetition;
            case BucketToBufferGadget::Dimension::set: return idx.set;
            case BucketToBufferGadget::Dimension::segment: return idx.segment;
            case BucketToBufferGadget::Dimension::slice: return 0;
            case BucketToBufferGadget::Dimension::none: return 0;
            default: throw std::runtime_error("Invalid enum encountered");
            }
        }
    }

    BufferKey BucketToBufferGadget::getKey(const ISMRMRD::EncodingCounters& idx) const {
        BufferKey key = idx;
        clear(N_dimension, key);
        clear(S_dimension, key);
        if (!split_slices)
            key.slice = 0;
        if (ignore_segment)
            key.segment = 0;
        return key;
    }

    namespace {
        uint16_t getSizeFromDimension(BucketToBufferGadget::Dimension dimension, const AcquisitionBucketStats& stats) {
            switch (dimension) {
            case BucketToBufferGadget::Dimension::phase: return *stats.phase.rbegin() - *stats.phase.begin() + 1;
            case BucketToBufferGadget::Dimension::contrast:
                return *stats.contrast.rbegin() - *stats.contrast.begin() + 1;
            case BucketToBufferGadget::Dimension::repetition:
                return *stats.repetition.rbegin() - *stats.repetition.begin() + 1;
            case BucketToBufferGadget::Dimension::set: return *stats.set.rbegin() - *stats.set.begin() + 1;
            case BucketToBufferGadget::Dimension::segment:
            case BucketToBufferGadget::Dimension::average: return *stats.average.rbegin() - *stats.average.begin() + 1;
            case BucketToBufferGadget::Dimension::slice: return *stats.slice.rbegin() - *stats.slice.begin() + 1;
            case BucketToBufferGadget::Dimension::none:; return 1;
            default: throw std::runtime_error("Illegal enum value.");
            }
        }
    }

    IsmrmrdDataBuffered BucketToBufferGadget::makeDataBuffer(const ISMRMRD::AcquisitionHeader& acqhdr,
        ISMRMRD::Encoding encoding, const AcquisitionBucketStats& stats, bool forref) const {
        IsmrmrdDataBuffered buffer;

        // Allocate the reference data array
        // 7D,  fixed order [E0, E1, E2, CHA, N, S, LOC]
        // 11D, fixed order [E0, E1, E2, CHA, SLC, PHS, CON, REP, SET, SEG, AVE]
        const uint16_t NE0 = getNE0(acqhdr, encoding);

        uint16_t NE1 = getNE1(encoding, stats, forref);

        uint16_t NE2 = getNE2(encoding, stats, forref);

        uint16_t NCHA = acqhdr.active_channels;

        uint16_t NLOC = getNLOC(encoding, stats);

        uint16_t NN = getSizeFromDimension(N_dimension, stats);

        uint16_t NS = getSizeFromDimension(S_dimension, stats);


        GDEBUG_CONDITION_STREAM(verbose, "Data dimensions [RO E1 E2 CHA N S SLC] : ["
                                             << NE0 << " " << NE1 << " " << NE2 << " " << NCHA << " " << NN << " " << NS
                                             << " " << NLOC << "]");

        // Allocate the array for the data
        buffer.data_ = hoNDArray<std::complex<float>>(NE0, NE1, NE2, NCHA, NN, NS, NLOC);
        clear(&buffer.data_);

        // Allocate the array for the headers
        buffer.headers_ = hoNDArray<ISMRMRD::AcquisitionHeader>(NE1, NE2, NN, NS, NLOC);

        // Allocate the array for the trajectories
        uint16_t TRAJDIM = acqhdr.trajectory_dimensions;
        if (TRAJDIM > 0) {
            buffer.trajectory_ = hoNDArray<float>(TRAJDIM, NE0, NE1, NE2, NN, NS, NLOC);
            clear(*buffer.trajectory_);
        }
        return buffer;
    }

    uint16_t BucketToBufferGadget::getNLOC(
        const ISMRMRD::Encoding& encoding, const AcquisitionBucketStats& stats) const {
        uint16_t NLOC;
        if (split_slices) {
            NLOC = 1;
        } else {
            if (encoding.encodingLimits.slice.is_present()) {
                NLOC = encoding.encodingLimits.slice->maximum - encoding.encodingLimits.slice->minimum + 1;
            } else {
                NLOC = 1;
            }

            // if the AcquisitionAccumulateTriggerGadget sort by SLC, then the stats should be used to determine NLOC
            size_t NLOC_received = *stats.slice.rbegin() - *stats.slice.begin() + 1;
            if (NLOC_received < NLOC) {
                NLOC = NLOC_received;
            }
        }
        return NLOC;
    }

    uint16_t BucketToBufferGadget::getNE2(
        const ISMRMRD::Encoding& encoding, const AcquisitionBucketStats& stats, bool forref) const {
        uint16_t NE2;
        if (((encoding.trajectory == ISMRMRD::TrajectoryType::CARTESIAN))
            || (encoding.trajectory == ISMRMRD::TrajectoryType::EPI)) {
            if (encoding.parallelImaging) {
                if (forref
                    && (encoding.parallelImaging.get().calibrationMode.get() == "separate"
                        || encoding.parallelImaging.get().calibrationMode.get() == "external")) {
                    NE2 = encoding.encodingLimits.kspace_encoding_step_2->maximum
                          - encoding.encodingLimits.kspace_encoding_step_2->minimum + 1;
                } else {
                    NE2 = encoding.encodedSpace.matrixSize.z;
                }
            } else {
                if (encoding.encodingLimits.kspace_encoding_step_2.is_present()) {
                    NE2 = encoding.encodingLimits.kspace_encoding_step_2->maximum
                          - encoding.encodingLimits.kspace_encoding_step_2->minimum + 1;
                } else {
                    NE2 = encoding.encodedSpace.matrixSize.z;
                }
            }
        } else {
            if (encoding.encodingLimits.kspace_encoding_step_2.is_present()) {
                NE2 = encoding.encodingLimits.kspace_encoding_step_2->maximum
                      - encoding.encodingLimits.kspace_encoding_step_2->minimum + 1;
            } else {
                NE2 = *stats.kspace_encode_step_2.rbegin() - *stats.kspace_encode_step_2.begin() + 1;
            }
        }
        return NE2;
    }

    uint16_t BucketToBufferGadget::getNE1(
        const ISMRMRD::Encoding& encoding, const AcquisitionBucketStats& stats, bool forref) const {
        uint16_t NE1;
        if (((encoding.trajectory == ISMRMRD::TrajectoryType::CARTESIAN))
            || (encoding.trajectory == ISMRMRD::TrajectoryType::EPI)) {
            if (encoding.parallelImaging) {
                if (forref
                    && (encoding.parallelImaging.get().calibrationMode.get() == "separate"
                        || encoding.parallelImaging.get().calibrationMode.get() == "external")) {
                    NE1 = *stats.kspace_encode_step_1.rbegin() - *stats.kspace_encode_step_1.begin() + 1;
                } else {
                    NE1 = encoding.encodedSpace.matrixSize.y;
                }
            } else {
                if (encoding.encodingLimits.kspace_encoding_step_1.is_present()) {
                    NE1 = encoding.encodingLimits.kspace_encoding_step_1->maximum
                          - encoding.encodingLimits.kspace_encoding_step_1->minimum + 1;
                } else {
                    NE1 = encoding.encodedSpace.matrixSize.y;
                }
            }
        } else {
            if (encoding.encodingLimits.kspace_encoding_step_1.is_present()) {
                NE1 = encoding.encodingLimits.kspace_encoding_step_1->maximum
                      - encoding.encodingLimits.kspace_encoding_step_1->minimum + 1;
            } else {
                NE1 = *stats.kspace_encode_step_1.rbegin() - *stats.kspace_encode_step_1.begin() + 1;
            }
        }
        return NE1;
    }

    uint16_t BucketToBufferGadget::getNE0(
        const ISMRMRD::AcquisitionHeader& acqhdr, const ISMRMRD::Encoding& encoding) const {
        uint16_t NE0;
        if (((encoding.trajectory == ISMRMRD::TrajectoryType::CARTESIAN))
            || (encoding.trajectory == ISMRMRD::TrajectoryType::EPI)) {
            // if separate or external calibration mode, using the acq length for NE0
            if (encoding.parallelImaging) {
                NE0 = acqhdr.number_of_samples;
            } else {
                NE0 = acqhdr.number_of_samples - acqhdr.discard_pre - acqhdr.discard_post;
            }
        } else {
            NE0 = acqhdr.number_of_samples - acqhdr.discard_pre - acqhdr.discard_post;
        }
        return NE0;
    }
    namespace {
        template <class DIMSTRUCT> auto xyz_to_vector(const DIMSTRUCT& dimstruct) {
            std::array<decltype(dimstruct.x), 3> result = { dimstruct.x, dimstruct.y, dimstruct.z };
            return result;
        }
    }

    SamplingDescription BucketToBufferGadget::createSamplingDescription(const ISMRMRD::Encoding& encoding,
        const AcquisitionBucketStats& stats, const ISMRMRD::AcquisitionHeader& acqhdr, bool forref) const {
        auto sampling            = SamplingDescription();
        sampling.encoded_FOV_    = xyz_to_vector(encoding.encodedSpace.fieldOfView_mm);
        sampling.encoded_matrix_ = xyz_to_vector(encoding.encodedSpace.matrixSize);
        sampling.recon_FOV_      = xyz_to_vector(encoding.reconSpace.fieldOfView_mm);
        sampling.recon_matrix_   = xyz_to_vector(encoding.reconSpace.matrixSize);

        sampling.sampling_limits_[0].min_ = 0;
        sampling.sampling_limits_[0].max_ = acqhdr.number_of_samples-1;
        sampling.sampling_limits_[0].center_ = acqhdr.number_of_samples/2;

        sampling.sampling_limits_[1].min_ = encoding.encodingLimits.kspace_encoding_step_1->minimum;
        sampling.sampling_limits_[1].max_ = encoding.encodingLimits.kspace_encoding_step_1->maximum;
        sampling.sampling_limits_[1].center_ = encoding.encodingLimits.kspace_encoding_step_1->center;

        sampling.sampling_limits_[2].min_ = encoding.encodingLimits.kspace_encoding_step_2->minimum;
        sampling.sampling_limits_[2].max_ = encoding.encodingLimits.kspace_encoding_step_2->maximum;
        sampling.sampling_limits_[2].center_ = encoding.encodingLimits.kspace_encoding_step_2->center;

        if (verbose) {
            GDEBUG_STREAM("Encoding space : " << acqhdr.encoding_space_ref << " - "
                          << int(encoding.trajectory) << " - FOV : [ " << encoding.encodedSpace.fieldOfView_mm.x << " "
                          << encoding.encodedSpace.fieldOfView_mm.y << " " << encoding.encodedSpace.fieldOfView_mm.z
                          << " ] "
                          << " - Matris size : [ " << encoding.encodedSpace.matrixSize.x << " "
                          << encoding.encodedSpace.matrixSize.y << " " << encoding.encodedSpace.matrixSize.z << " ] ");

            GDEBUG_STREAM("Sampling limits : "
                          << "- RO : [ " << sampling.sampling_limits_[0].min_ << " "
                          << sampling.sampling_limits_[0].center_ << " " << sampling.sampling_limits_[0].max_
                          << " ] - E1 : [ " << sampling.sampling_limits_[1].min_ << " "
                          << sampling.sampling_limits_[1].center_ << " " << sampling.sampling_limits_[1].max_
                          << " ] - E2 : [ " << sampling.sampling_limits_[2].min_ << " "
                          << sampling.sampling_limits_[2].center_ << " " << sampling.sampling_limits_[2].max_ << " ]");
        }

        // For cartesian trajectories, assume that any oversampling has been removed.
        if (encoding.trajectory == ISMRMRD::TrajectoryType::CARTESIAN) {
            sampling.encoded_FOV_[0]    = encoding.reconSpace.fieldOfView_mm.x;
            sampling.encoded_matrix_[0] = encoding.reconSpace.matrixSize.x;
        } else {
            sampling.encoded_FOV_[0]    = encoding.encodedSpace.fieldOfView_mm.x;
            sampling.encoded_matrix_[0] = encoding.encodedSpace.matrixSize.x;
        }

        // For cartesian trajectories, assume that any oversampling has been removed.
        if (((encoding.trajectory == ISMRMRD::TrajectoryType::CARTESIAN))
            || (encoding.trajectory == ISMRMRD::TrajectoryType::EPI)) {
            sampling.sampling_limits_[0].min_    = acqhdr.discard_pre;
            sampling.sampling_limits_[0].max_    = acqhdr.number_of_samples - acqhdr.discard_post - 1;
            sampling.sampling_limits_[0].center_ = acqhdr.number_of_samples / 2;
        } else {
            sampling.sampling_limits_[0].min_    = 0;
            sampling.sampling_limits_[0].max_    = encoding.encodedSpace.matrixSize.x - 1;
            sampling.sampling_limits_[0].center_ = encoding.encodedSpace.matrixSize.x / 2;
        }

        // if the scan is cartesian
        if (((encoding.trajectory == ISMRMRD::TrajectoryType::CARTESIAN)
                && (!forref || (forref && (encoding.parallelImaging.get().calibrationMode.get() == "embedded"))))
            || ((encoding.trajectory == ISMRMRD::TrajectoryType::EPI) && !forref)) {
            int16_t space_matrix_offset_E1 = 0;
            if (encoding.encodingLimits.kspace_encoding_step_1.is_present()) {
                space_matrix_offset_E1 = (int16_t)encoding.encodedSpace.matrixSize.y / 2
                                         - (int16_t)encoding.encodingLimits.kspace_encoding_step_1->center;
            }

            int16_t space_matrix_offset_E2 = 0;
            if (encoding.encodingLimits.kspace_encoding_step_2.is_present() && encoding.encodedSpace.matrixSize.z > 1) {
                space_matrix_offset_E2 = (int16_t)encoding.encodedSpace.matrixSize.z / 2
                                         - (int16_t)encoding.encodingLimits.kspace_encoding_step_2->center;
            }

            // E1
            sampling.sampling_limits_[1].min_
                = encoding.encodingLimits.kspace_encoding_step_1->minimum + space_matrix_offset_E1;
            sampling.sampling_limits_[1].max_
                = encoding.encodingLimits.kspace_encoding_step_1->maximum + space_matrix_offset_E1;
            sampling.sampling_limits_[1].center_ = sampling.encoded_matrix_[1] / 2;

            GADGET_CHECK_THROW(sampling.sampling_limits_[1].min_ < encoding.encodedSpace.matrixSize.y);
            GADGET_CHECK_THROW(sampling.sampling_limits_[1].max_ >= sampling.sampling_limits_[1].min_);
            GADGET_CHECK_THROW(sampling.sampling_limits_[1].center_ >= sampling.sampling_limits_[1].min_);
            GADGET_CHECK_THROW(sampling.sampling_limits_[1].center_ <= sampling.sampling_limits_[1].max_);

            // E2
            sampling.sampling_limits_[2].min_
                = encoding.encodingLimits.kspace_encoding_step_2->minimum + space_matrix_offset_E2;
            sampling.sampling_limits_[2].max_
                = encoding.encodingLimits.kspace_encoding_step_2->maximum + space_matrix_offset_E2;
            sampling.sampling_limits_[2].center_ = sampling.encoded_matrix_[2] / 2;

            GADGET_CHECK_THROW(sampling.sampling_limits_[2].min_ < encoding.encodedSpace.matrixSize.y);
            GADGET_CHECK_THROW(sampling.sampling_limits_[2].max_ >= sampling.sampling_limits_[2].min_);
            GADGET_CHECK_THROW(sampling.sampling_limits_[2].center_ >= sampling.sampling_limits_[2].min_);
            GADGET_CHECK_THROW(sampling.sampling_limits_[2].center_ <= sampling.sampling_limits_[2].max_);
        } else {
            sampling.sampling_limits_[1].min_    = encoding.encodingLimits.kspace_encoding_step_1->minimum;
            sampling.sampling_limits_[1].max_    = encoding.encodingLimits.kspace_encoding_step_1->maximum;
            sampling.sampling_limits_[1].center_ = encoding.encodingLimits.kspace_encoding_step_1->center;

            sampling.sampling_limits_[2].min_    = encoding.encodingLimits.kspace_encoding_step_2->minimum;
            sampling.sampling_limits_[2].max_    = encoding.encodingLimits.kspace_encoding_step_2->maximum;
            sampling.sampling_limits_[2].center_ = encoding.encodingLimits.kspace_encoding_step_2->center;
        }

        return sampling;
    }

    void BucketToBufferGadget::add_acquisition(IsmrmrdDataBuffered& dataBuffer, const Core::Acquisition& acq,
        ISMRMRD::Encoding encoding, const AcquisitionBucketStats& stats, bool forref) {

        // The acquisition header and data
        const auto& acqhdr  = std::get<ISMRMRD::AcquisitionHeader>(acq);
        const auto& acqdata = std::get<hoNDArray<std::complex<float>>>(acq);
        // we make one for the trajectory down below if we need it

        uint16_t NE0  = (uint16_t)dataBuffer.data_.get_size(0);
        uint16_t NE1  = (uint16_t)dataBuffer.data_.get_size(1);
        uint16_t NE2  = (uint16_t)dataBuffer.data_.get_size(2);
        uint16_t NCHA = (uint16_t)dataBuffer.data_.get_size(3);
        uint16_t NN   = (uint16_t)dataBuffer.data_.get_size(4);
        uint16_t NS   = (uint16_t)dataBuffer.data_.get_size(5);
        uint16_t NLOC = (uint16_t)dataBuffer.data_.get_size(6);

        const size_t slice_loc = split_slices || NLOC == 1 ? 0 : acqhdr.idx.slice;

        // Stuff the data
        uint16_t npts_to_copy = acqhdr.number_of_samples - acqhdr.discard_pre - acqhdr.discard_post;
        long long offset;
        if (encoding.trajectory == ISMRMRD::TrajectoryType::CARTESIAN
            || encoding.trajectory == ISMRMRD::TrajectoryType::EPI) {
            if ((acqhdr.number_of_samples == dataBuffer.data_.get_size(0))
                && ((acqhdr.center_sample == acqhdr.number_of_samples / 2) || acqhdr.center_sample >= acqhdr.number_of_samples)) // acq has been corrected for center , e.g. by asymmetric handling
            {
                offset = acqhdr.discard_pre;
            } else {
                offset = (long long)dataBuffer.sampling_.sampling_limits_[0].center_ - (long long)acqhdr.center_sample;
            }
        } else {
            // TODO what about EPI with asymmetric readouts?
            // TODO any other sort of trajectory?
            offset = 0;
        }

        long long roffset = (long long)dataBuffer.data_.get_size(0) - npts_to_copy - offset;

        if ((offset < 0) | (roffset < 0)) {
            throw std::runtime_error("Acquired reference data does not fit into the reference data buffer.\n");
        }

        uint16_t NUsed = (uint16_t)getDimensionKey(N_dimension, acqhdr.idx);
        if (NUsed >= NN)
            NUsed = NN - 1;

        uint16_t SUsed = (uint16_t)getDimensionKey(S_dimension, acqhdr.idx);
        if (SUsed >= NS)
            SUsed = NS - 1;

        int16_t e1 = (int16_t)acqhdr.idx.kspace_encode_step_1;
        int16_t e2 = (int16_t)acqhdr.idx.kspace_encode_step_2;

        bool is_cartesian_sampling = (encoding.trajectory == ISMRMRD::TrajectoryType::CARTESIAN);
        bool is_epi_sampling       = (encoding.trajectory == ISMRMRD::TrajectoryType::EPI);
        if (is_cartesian_sampling || is_epi_sampling) {
            if (!forref || (forref && (encoding.parallelImaging.get().calibrationMode.get() == "embedded"))) {
                // compute the center offset for E1 and E2
                int16_t space_matrix_offset_E1 = 0;
                if (encoding.encodingLimits.kspace_encoding_step_1.is_present()) {
                    space_matrix_offset_E1 = (int16_t)encoding.encodedSpace.matrixSize.y / 2
                                             - (int16_t)encoding.encodingLimits.kspace_encoding_step_1->center;
                }

                int16_t space_matrix_offset_E2 = 0;
                if (encoding.encodingLimits.kspace_encoding_step_2.is_present()
                    && encoding.encodedSpace.matrixSize.z > 1) {
                    space_matrix_offset_E2 = (int16_t)encoding.encodedSpace.matrixSize.z / 2
                                             - (int16_t)encoding.encodingLimits.kspace_encoding_step_2->center;
                }

                // compute the used e1 and e2 indices and make sure they are in the valid range
                e1 = (int16_t)acqhdr.idx.kspace_encode_step_1 + space_matrix_offset_E1;
                e2 = (int16_t)acqhdr.idx.kspace_encode_step_2 + space_matrix_offset_E2;
            }

            // for external or separate mode, it is possible the starting numbers of ref lines are not zero, therefore
            // it is needed to subtract the staring ref line number because the ref array size is set up by the actual
            // number of lines acquired only assumption for external or separate ref line mode is that all ref lines are
            // numbered sequentially the acquisition order of ref line can be arbitrary
            if (forref
                && ((encoding.parallelImaging.get().calibrationMode.get() == "separate")
                    || (encoding.parallelImaging.get().calibrationMode.get() == "external"))) {
                if (*stats.kspace_encode_step_1.begin() > 0) {
                    e1 = acqhdr.idx.kspace_encode_step_1 - *stats.kspace_encode_step_1.begin();
                }

                if (*stats.kspace_encode_step_2.begin() > 0) {
                    e2 = acqhdr.idx.kspace_encode_step_2 - *stats.kspace_encode_step_2.begin();
                }
            }

            if (e1 < 0 || e1 >= (int16_t)NE1) {
                // if the incoming line is outside the encoding limits, something is wrong
                GADGET_CHECK_THROW(
                    acqhdr.idx.kspace_encode_step_1 >= encoding.encodingLimits.kspace_encoding_step_1->minimum
                    && acqhdr.idx.kspace_encode_step_1 <= encoding.encodingLimits.kspace_encoding_step_1->maximum);

                // if the incoming line is inside encoding limits but outside the encoded matrix, do not include the
                // data
                GWARN_STREAM(
                    "incoming readout "
                    << acqhdr.scan_counter
                    << " is inside the encoding limits, but outside the encoded matrix for kspace_encode_step_1 : "
                    << e1 << " out of " << NE1);
                return;
            }

            if (e2 < 0 || e2 >= (int16_t)NE2) {
                GADGET_CHECK_THROW(
                    acqhdr.idx.kspace_encode_step_2 >= encoding.encodingLimits.kspace_encoding_step_2->minimum
                    && acqhdr.idx.kspace_encode_step_2 <= encoding.encodingLimits.kspace_encoding_step_2->maximum);

                GWARN_STREAM(
                    "incoming readout "
                    << acqhdr.scan_counter
                    << " is inside the encoding limits, but outside the encoded matrix for kspace_encode_step_2 : "
                    << e2 << " out of " << NE2);
                return;
            }
        }

        std::complex<float>* pData = &dataBuffer.data_(offset, e1, e2, 0, NUsed, SUsed, slice_loc);

        for (uint16_t cha = 0; cha < NCHA; cha++) {
            auto dataptr = pData + cha * NE0 * NE1 * NE2;
            auto fromptr = &acqdata(acqhdr.discard_pre, cha);
            std::copy(fromptr, fromptr + npts_to_copy, dataptr);
        }

        dataBuffer.headers_(e1, e2, NUsed, SUsed, slice_loc) = acqhdr;

        if (acqhdr.trajectory_dimensions > 0) {

            const auto& acqtraj = *std::get<Core::optional<hoNDArray<float>>>(acq); // TODO do we need to check this?

            float* trajptr = &(*dataBuffer.trajectory_)(0, offset, e1, e2, NUsed, SUsed, slice_loc);
            auto* fromptr  = &acqtraj(0, acqhdr.discard_pre);
            std::copy(fromptr, fromptr + npts_to_copy * acqhdr.trajectory_dimensions, trajptr);
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

}
