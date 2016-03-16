
#include "GenericReconBucketToBufferGadget.h"
#include "GadgetIsmrmrdReadWrite.h"
#include "mri_core_data.h"
#include "hoNDArray_elemwise.h"
#include "hoNDArray_reductions.h"

namespace Gadgetron {

    GenericReconBucketToBufferGadget::GenericReconBucketToBufferGadget()
    {
    }

    GenericReconBucketToBufferGadget::~GenericReconBucketToBufferGadget()
    {
    }

    int GenericReconBucketToBufferGadget::process_config(ACE_Message_Block* mb)
    {
        return BaseClass::process_config(mb);
    }

    void GenericReconBucketToBufferGadget::allocateDataArrays(IsmrmrdDataBuffered & dataBuffer, ISMRMRD::AcquisitionHeader & acqhdr, ISMRMRD::Encoding encoding, IsmrmrdAcquisitionBucketStats & stats, bool forref)
    {
        if (dataBuffer.data_.get_number_of_elements() == 0)
        {
            uint16_t NE0;
            if (encoding.trajectory.compare("cartesian") == 0)
            {
                // if seperate or external calibration mode, using the acq length for NE0
                if (encoding.parallelImaging)
                {
                    NE0 = acqhdr.number_of_samples;
                }
                else
                {
                    NE0 = acqhdr.number_of_samples - acqhdr.discard_pre - acqhdr.discard_post;
                }
            }
            else
            {
                NE0 = acqhdr.number_of_samples - acqhdr.discard_pre - acqhdr.discard_post;
            }

            uint16_t NE1;
            if (encoding.trajectory.compare("cartesian") == 0)
            {
                if (encoding.parallelImaging)
                {
                    if (forref && (encoding.parallelImaging.get().calibrationMode.get() == "separate"
                        || encoding.parallelImaging.get().calibrationMode.get() == "external"))
                    {
                        NE1 = *stats.kspace_encode_step_1.rbegin() - *stats.kspace_encode_step_1.begin() + 1;
                    }
                    else
                    {
                        // make sure squared pixel recon
                        float dx = encoding.encodedSpace.fieldOfView_mm.x / encoding.encodedSpace.matrixSize.x;

                        NE1 = (size_t)std::floor(encoding.encodedSpace.fieldOfView_mm.y / dx + 0.5);
                        if (NE1 < encoding.encodedSpace.matrixSize.y)
                        {
                            NE1 = encoding.encodedSpace.matrixSize.y;
                        }
                    }
                }
                else
                {
                    if (encoding.encodingLimits.kspace_encoding_step_1.is_present())
                    {
                        NE1 = encoding.encodingLimits.kspace_encoding_step_1->maximum - encoding.encodingLimits.kspace_encoding_step_1->minimum + 1;
                    }
                    else
                    {
                        NE1 = encoding.encodedSpace.matrixSize.y;
                    }
                }
            }
            else {
                if (encoding.encodingLimits.kspace_encoding_step_1.is_present()) {
                    NE1 = encoding.encodingLimits.kspace_encoding_step_1->maximum - encoding.encodingLimits.kspace_encoding_step_1->minimum + 1;
                }
                else {
                    NE1 = *stats.kspace_encode_step_1.rbegin() - *stats.kspace_encode_step_1.begin() + 1;
                }
            }

            uint16_t NE2;
            if (encoding.trajectory.compare("cartesian") == 0)
            {
                if (encoding.parallelImaging)
                {
                    if (forref && (encoding.parallelImaging.get().calibrationMode.get() == "separate" || encoding.parallelImaging.get().calibrationMode.get() == "external"))
                    {
                        NE2 = encoding.encodingLimits.kspace_encoding_step_2->maximum - encoding.encodingLimits.kspace_encoding_step_2->minimum + 1;
                    }
                    else
                    {
                        NE2 = encoding.encodedSpace.matrixSize.z;
                    }
                }
                else
                {
                    if (encoding.encodingLimits.kspace_encoding_step_2.is_present())
                    {
                        NE2 = encoding.encodingLimits.kspace_encoding_step_2->maximum - encoding.encodingLimits.kspace_encoding_step_2->minimum + 1;
                    }
                    else
                    {
                        NE2 = encoding.encodedSpace.matrixSize.z;
                    }
                }
            }
            else {
                if (encoding.encodingLimits.kspace_encoding_step_2.is_present())
                {
                    NE2 = encoding.encodingLimits.kspace_encoding_step_2->maximum - encoding.encodingLimits.kspace_encoding_step_2->minimum + 1;
                }
                else
                {
                    NE2 = *stats.kspace_encode_step_2.rbegin() - *stats.kspace_encode_step_2.begin() + 1;
                }
            }

            uint16_t NCHA = acqhdr.active_channels;

            uint16_t NLOC;
            if (split_slices_)
            {
                NLOC = 1;
            }
            else
            {
                if (encoding.encodingLimits.slice.is_present())
                {
                    NLOC = encoding.encodingLimits.slice->maximum - encoding.encodingLimits.slice->minimum + 1;
                }
                else
                {
                    NLOC = *stats.slice.rbegin() - *stats.slice.begin() + 1;
                }
            }

            uint16_t NN;
            switch (N_)
            {
            case PHASE:
                NN = *stats.phase.rbegin() - *stats.phase.begin() + 1;
                break;
            case CONTRAST:
                NN = *stats.contrast.rbegin() - *stats.contrast.begin() + 1;
                break;
            case REPETITION:
                NN = *stats.repetition.rbegin() - *stats.repetition.begin() + 1;
                break;
            case SET:
                NN = *stats.set.rbegin() - *stats.set.begin() + 1;
                break;
            case SEGMENT:
                NN = *stats.segment.rbegin() - *stats.segment.begin() + 1;
                break;
            case AVERAGE:
                NN = *stats.average.rbegin() - *stats.average.begin() + 1;
                break;
            case SLICE:
                NN = *stats.slice.rbegin() - *stats.slice.begin() + 1;
                break;
            default:
                NN = 1;
            }

            uint16_t NS;
            switch (S_)
            {
            case PHASE:
                NS = *stats.phase.rbegin() - *stats.phase.begin() + 1;
                break;
            case CONTRAST:
                NS = *stats.contrast.rbegin() - *stats.contrast.begin() + 1;
                break;
            case REPETITION:
                NS = *stats.repetition.rbegin() - *stats.repetition.begin() + 1;
                break;
            case SET:
                NS = *stats.set.rbegin() - *stats.set.begin() + 1;
                break;
            case SEGMENT:
                NS = *stats.segment.rbegin() - *stats.segment.begin() + 1;
                break;
            case AVERAGE:
                NS = *stats.average.rbegin() - *stats.average.begin() + 1;
                break;
            case SLICE:
                NS = *stats.slice.rbegin() - *stats.slice.begin() + 1;
                break;
            default:
                NS = 1;
            }

            //Allocate the array for the data
            dataBuffer.data_.create(NE0, NE1, NE2, NCHA, NN, NS, NLOC);
            clear(&dataBuffer.data_);

            //Allocate the array for the headers
            dataBuffer.headers_.create(NE1, NE2, NN, NS, NLOC);

            //Allocate the array for the trajectories
            uint16_t TRAJDIM = acqhdr.trajectory_dimensions;
            if (TRAJDIM > 0)
            {
                dataBuffer.trajectory_ = hoNDArray<float>(TRAJDIM, NE0, NE1, NE2, NN, NS, NLOC);
                clear(dataBuffer.trajectory_.get_ptr());
            }
        }
    }

    void GenericReconBucketToBufferGadget::fillSamplingDescription(SamplingDescription & sampling, ISMRMRD::Encoding & encoding, IsmrmrdAcquisitionBucketStats & stats, ISMRMRD::AcquisitionHeader& acqhdr, bool forref)
    {
        if (encoding.trajectory.compare("cartesian") == 0)
        {
            sampling.encoded_FOV_[0] = encoding.reconSpace.fieldOfView_mm.x;
            sampling.encoded_matrix_[0] = encoding.reconSpace.matrixSize.x;
        }
        else
        {
            sampling.encoded_FOV_[0] = encoding.encodedSpace.fieldOfView_mm.x;
            sampling.encoded_matrix_[0] = encoding.encodedSpace.matrixSize.x;
        }

        sampling.encoded_FOV_[1] = encoding.encodedSpace.fieldOfView_mm.y;
        sampling.encoded_FOV_[2] = encoding.encodedSpace.fieldOfView_mm.z;

        sampling.encoded_matrix_[1] = encoding.encodedSpace.matrixSize.y;
        sampling.encoded_matrix_[2] = encoding.encodedSpace.matrixSize.z;

        sampling.recon_FOV_[0] = encoding.reconSpace.fieldOfView_mm.x;
        sampling.recon_FOV_[1] = encoding.reconSpace.fieldOfView_mm.y;
        sampling.recon_FOV_[2] = encoding.reconSpace.fieldOfView_mm.z;

        sampling.recon_matrix_[0] = encoding.reconSpace.matrixSize.x;
        sampling.recon_matrix_[1] = encoding.reconSpace.matrixSize.y;
        sampling.recon_matrix_[2] = encoding.reconSpace.matrixSize.z;

        if (encoding.trajectory.compare("cartesian") == 0)
        {
            sampling.sampling_limits_[0].min_ = acqhdr.discard_pre;
            sampling.sampling_limits_[0].max_ = acqhdr.number_of_samples - acqhdr.discard_post - 1;
            sampling.sampling_limits_[0].center_ = acqhdr.number_of_samples / 2;
        }
        else
        {
            sampling.sampling_limits_[0].min_ = 0;
            sampling.sampling_limits_[0].max_ = encoding.encodedSpace.matrixSize.x - 1;
            sampling.sampling_limits_[0].center_ = encoding.encodedSpace.matrixSize.x / 2;
        }

        if (!forref || (forref && (encoding.parallelImaging.get().calibrationMode.get() == "embedded")))
        {
            int16_t space_matrix_offset_E1 = 0;
            if (encoding.encodingLimits.kspace_encoding_step_1.is_present())
            {
                space_matrix_offset_E1 = (int16_t)encoding.encodedSpace.matrixSize.y / 2 - (int16_t)encoding.encodingLimits.kspace_encoding_step_1->center;
            }

            int16_t space_matrix_offset_E2 = 0;
            if (encoding.encodingLimits.kspace_encoding_step_2.is_present() && encoding.encodedSpace.matrixSize.z > 1)
            {
                space_matrix_offset_E2 = (int16_t)encoding.encodedSpace.matrixSize.z / 2 - (int16_t)encoding.encodingLimits.kspace_encoding_step_2->center;
            }

            {
                sampling.sampling_limits_[1].min_ = encoding.encodingLimits.kspace_encoding_step_1->minimum + space_matrix_offset_E1;
                sampling.sampling_limits_[1].max_ = encoding.encodingLimits.kspace_encoding_step_1->maximum + space_matrix_offset_E1;
                sampling.sampling_limits_[1].center_ = sampling.encoded_matrix_[1] / 2;

                if (sampling.sampling_limits_[1].min_ < 0) sampling.sampling_limits_[1].min_ = 0;
                if (sampling.sampling_limits_[1].min_ >= encoding.encodedSpace.matrixSize.y) sampling.sampling_limits_[1].min_ = encoding.encodedSpace.matrixSize.y - 1;

                if (sampling.sampling_limits_[1].max_ < 0) sampling.sampling_limits_[1].max_ = 0;
                if (sampling.sampling_limits_[1].max_ >= encoding.encodedSpace.matrixSize.y) sampling.sampling_limits_[1].max_ = encoding.encodedSpace.matrixSize.y - 1;

                if (sampling.sampling_limits_[1].min_ > sampling.sampling_limits_[1].max_) sampling.sampling_limits_[1].min_ = sampling.sampling_limits_[1].max_;

                if (sampling.sampling_limits_[1].center_ < sampling.sampling_limits_[1].min_) sampling.sampling_limits_[1].center_ = sampling.sampling_limits_[1].min_;
                if (sampling.sampling_limits_[1].center_ > sampling.sampling_limits_[1].max_) sampling.sampling_limits_[1].center_ = sampling.sampling_limits_[1].max_;
            }

            {
                sampling.sampling_limits_[2].min_ = encoding.encodingLimits.kspace_encoding_step_2->minimum + space_matrix_offset_E2;
                sampling.sampling_limits_[2].max_ = encoding.encodingLimits.kspace_encoding_step_2->maximum + space_matrix_offset_E2;
                sampling.sampling_limits_[2].center_ = sampling.encoded_matrix_[2] / 2;

                if (sampling.sampling_limits_[2].min_ < 0) sampling.sampling_limits_[2].min_ = 0;
                if (sampling.sampling_limits_[2].min_ >= encoding.encodedSpace.matrixSize.z) sampling.sampling_limits_[2].min_ = encoding.encodedSpace.matrixSize.z - 1;

                if (sampling.sampling_limits_[2].max_ < 0) sampling.sampling_limits_[2].max_ = 0;
                if (sampling.sampling_limits_[2].max_ >= encoding.encodedSpace.matrixSize.z) sampling.sampling_limits_[2].max_ = encoding.encodedSpace.matrixSize.z - 1;

                if (sampling.sampling_limits_[2].min_ > sampling.sampling_limits_[2].max_) sampling.sampling_limits_[2].min_ = sampling.sampling_limits_[2].max_;

                if (sampling.sampling_limits_[2].center_ < sampling.sampling_limits_[2].min_) sampling.sampling_limits_[2].center_ = sampling.sampling_limits_[2].min_;
                if (sampling.sampling_limits_[2].center_ > sampling.sampling_limits_[2].max_) sampling.sampling_limits_[2].center_ = sampling.sampling_limits_[2].max_;
            }
        }
        else
        {
            sampling.sampling_limits_[1].min_ = encoding.encodingLimits.kspace_encoding_step_1->minimum;
            sampling.sampling_limits_[1].max_ = encoding.encodingLimits.kspace_encoding_step_1->maximum;
            sampling.sampling_limits_[1].center_ = encoding.encodingLimits.kspace_encoding_step_1->center;

            sampling.sampling_limits_[2].min_ = encoding.encodingLimits.kspace_encoding_step_2->minimum;
            sampling.sampling_limits_[2].max_ = encoding.encodingLimits.kspace_encoding_step_2->maximum;
            sampling.sampling_limits_[2].center_ = encoding.encodingLimits.kspace_encoding_step_2->center;
        }
    }

    void GenericReconBucketToBufferGadget::stuff(std::vector<IsmrmrdAcquisitionData>::iterator it, IsmrmrdDataBuffered & dataBuffer, ISMRMRD::Encoding encoding, bool forref)
    {

        ISMRMRD::AcquisitionHeader & acqhdr = *it->head_->getObjectPtr();
        hoNDArray< std::complex<float> > & acqdata = *it->data_->getObjectPtr();

        size_t slice_loc;
        if (split_slices_)
        {
            slice_loc = 0;
        }
        else
        {
            slice_loc = acqhdr.idx.slice;
        }

        //Stuff the data
        uint16_t npts_to_copy = acqhdr.number_of_samples - acqhdr.discard_pre - acqhdr.discard_post;
        long long offset;
        if (encoding.trajectory.compare("cartesian") == 0)
        {
            if ((acqhdr.number_of_samples == dataBuffer.data_.get_size(0)) && (acqhdr.center_sample == acqhdr.number_of_samples / 2)) // acq has been corrected for center , e.g. by asymmetric handling
            {
                offset = acqhdr.discard_pre;
            }
            else
            {
                offset = (long long)dataBuffer.sampling_.sampling_limits_[0].center_ - (long long)acqhdr.center_sample;
            }
        }
        else
        {
            offset = 0;
        }
        long long roffset = (long long)dataBuffer.data_.get_size(0) - npts_to_copy - offset;

        if ((offset < 0) | (roffset < 0))
        {
            throw std::runtime_error("Acquired reference data does not fit into the reference data buffer.\n");
        }

        std::complex<float> *dataptr;
        uint16_t NE1 = (uint16_t)dataBuffer.data_.get_size(1);
        uint16_t NE2 = (uint16_t)dataBuffer.data_.get_size(2);
        uint16_t NCHA = (uint16_t)dataBuffer.data_.get_size(3);
        uint16_t NN = (uint16_t)dataBuffer.data_.get_size(4);
        uint16_t NS = (uint16_t)dataBuffer.data_.get_size(5);

        uint16_t NUsed = (uint16_t)getN(acqhdr.idx);
        if (NUsed >= NN) NUsed = NN - 1;

        uint16_t SUsed = (uint16_t)getS(acqhdr.idx);
        if (SUsed >= NS) SUsed = NS - 1;

        int16_t e1 = (int16_t)acqhdr.idx.kspace_encode_step_1;
        int16_t e2 = (int16_t)acqhdr.idx.kspace_encode_step_2;

        if (!forref || (forref && (encoding.parallelImaging.get().calibrationMode.get() == "embedded")))
        {
            // compute the center offset for E1 and E2
            int16_t space_matrix_offset_E1 = 0;
            if (encoding.encodingLimits.kspace_encoding_step_1.is_present())
            {
                space_matrix_offset_E1 = (int16_t)encoding.encodedSpace.matrixSize.y / 2 - (int16_t)encoding.encodingLimits.kspace_encoding_step_1->center;
            }

            int16_t space_matrix_offset_E2 = 0;
            if (encoding.encodingLimits.kspace_encoding_step_2.is_present() && encoding.encodedSpace.matrixSize.z > 1)
            {
                space_matrix_offset_E2 = (int16_t)encoding.encodedSpace.matrixSize.z / 2 - (int16_t)encoding.encodingLimits.kspace_encoding_step_2->center;
            }

            // compute the used e1 and e2 indices and make sure they are in the valid range
            e1 = (int16_t)acqhdr.idx.kspace_encode_step_1 + space_matrix_offset_E1;
            e2 = (int16_t)acqhdr.idx.kspace_encode_step_2 + space_matrix_offset_E2;
        }

        if (e1 < 0) e1 = 0;
        if (e1 >= NE1) e1 = NE1 - 1;

        if (e2 < 0) e2 = 0;
        if (e2 >= NE2) e2 = NE2 - 1;

        for (uint16_t cha = 0; cha < NCHA; cha++)
        {
            dataptr = &dataBuffer.data_(offset, e1, e2, cha, NUsed, SUsed, slice_loc);

            memcpy(dataptr, &acqdata(acqhdr.discard_pre, cha), sizeof(std::complex<float>)*npts_to_copy);
        }

        dataBuffer.headers_(e1, e2, NUsed, SUsed, slice_loc) = acqhdr;

        if (acqhdr.trajectory_dimensions > 0)
        {

            hoNDArray< float > & acqtraj = *it->traj_->getObjectPtr();  // TODO do we need to check this?

            float * trajptr;

            trajptr = &(*dataBuffer.trajectory_)(0, offset, e1, e2, NUsed, SUsed, slice_loc);

            memcpy(trajptr, &acqtraj(0, acqhdr.discard_pre), sizeof(float)*npts_to_copy*acqhdr.trajectory_dimensions);

        }
    }

    GADGET_FACTORY_DECLARE(GenericReconBucketToBufferGadget)
}
