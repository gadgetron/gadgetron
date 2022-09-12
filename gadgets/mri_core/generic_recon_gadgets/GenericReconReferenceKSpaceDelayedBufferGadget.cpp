
#include "GenericReconReferenceKSpaceDelayedBufferGadget.h"
#include <iomanip>

#include "hoNDArray_reductions.h"
#include "mri_core_def.h"

namespace Gadgetron {

    GenericReconReferenceKSpaceDelayedBufferGadget::GenericReconReferenceKSpaceDelayedBufferGadget() : BaseClass()
    {
    }

    GenericReconReferenceKSpaceDelayedBufferGadget::~GenericReconReferenceKSpaceDelayedBufferGadget()
    {
    }

    int GenericReconReferenceKSpaceDelayedBufferGadget::process_config(ACE_Message_Block* mb)
    {
        GADGET_CHECK_RETURN(BaseClass::process_config(mb) == GADGET_OK, GADGET_FAIL);

        ISMRMRD::IsmrmrdHeader h;
        try
        {
            deserialize(mb->rd_ptr(), h);
        }
        catch (...)
        {
            GDEBUG("Error parsing ISMRMRD Header");
        }

        if (!h.acquisitionSystemInformation)
        {
            GDEBUG("acquisitionSystemInformation not found in header. Bailing out");
            return GADGET_FAIL;
        }

        // -------------------------------------------------

        size_t NE = h.encoding.size();
        num_encoding_spaces_ = NE;
        GDEBUG_CONDITION_STREAM(verbose.value(), "Number of encoding spaces: " << NE);

        calib_mode_.resize(NE, ISMRMRD_noacceleration);
        SLC_.resize(NE, 0);

        ref_buf_.resize(NE);

        imaging_data_arrived_.resize(NE);

        for (size_t e = 0; e < h.encoding.size(); e++)
        {
            ISMRMRD::EncodingSpace e_space = h.encoding[e].encodedSpace;
            ISMRMRD::EncodingSpace r_space = h.encoding[e].reconSpace;
            ISMRMRD::EncodingLimits e_limits = h.encoding[e].encodingLimits;

            if (!h.encoding[e].parallelImaging)
            {
                GDEBUG_STREAM("Parallel Imaging section not found in header");
                calib_mode_[e] = ISMRMRD_noacceleration;
            }
            else
            {

                ISMRMRD::ParallelImaging p_imaging = *h.encoding[0].parallelImaging;
                std::string calib = *p_imaging.calibrationMode;

                bool separate = (calib.compare("separate") == 0);
                bool embedded = (calib.compare("embedded") == 0);
                bool external = (calib.compare("external") == 0);
                bool interleaved = (calib.compare("interleaved") == 0);
                bool other = (calib.compare("other") == 0);

                calib_mode_[e] = Gadgetron::ISMRMRD_noacceleration;
                if (p_imaging.accelerationFactor.kspace_encoding_step_1 > 1 || p_imaging.accelerationFactor.kspace_encoding_step_2 > 1)
                {
                    if (interleaved)
                        calib_mode_[e] = Gadgetron::ISMRMRD_interleaved;
                    else if (embedded)
                        calib_mode_[e] = Gadgetron::ISMRMRD_embedded;
                    else if (separate)
                        calib_mode_[e] = Gadgetron::ISMRMRD_separate;
                    else if (external)
                        calib_mode_[e] = Gadgetron::ISMRMRD_external;
                    else if (other)
                        calib_mode_[e] = Gadgetron::ISMRMRD_other;
                }
            }

            SLC_[e] = e_limits.slice.get().maximum + 1;
            imaging_data_arrived_[e].resize(SLC_[e], false);

            ref_buf_[e].resize(SLC_[e]);

            GDEBUG_STREAM("Econding space " << e << " has " << SLC_[e] << " slices, with acceleration mode being " << calib_mode_[e]);
        }

        return GADGET_OK;
    }

    int GenericReconReferenceKSpaceDelayedBufferGadget::process(Gadgetron::GadgetContainerMessage< ISMRMRD::AcquisitionHeader >* m1)
    {
        if (perform_timing.value()) { gt_timer_.start("GenericReconReferenceKSpaceDelayedBufferGadget::process"); }

        process_called_times_++;

        ISMRMRD::AcquisitionHeader* acq_header = m1->getObjectPtr();

        uint16_t espace = acq_header->encoding_space_ref;
        size_t slc = acq_header->idx.slice;

        // if acceleration mode is not separate, pass data along
        if(calib_mode_[espace] != Gadgetron::ISMRMRD_separate)
        {
            if (this->next()->putq(m1) == -1)
            {
                GERROR("GenericReconReferenceKSpaceDelayedBufferGadget::process, passing data on to next gadget");
                return GADGET_FAIL;
            }

            return GADGET_OK;
        }

        // whether this read out line is ref
        if (ISMRMRD::FlagBit(ISMRMRD::ISMRMRD_ACQ_IS_PARALLEL_CALIBRATION).isSet(acq_header->flags) ||
            ISMRMRD::FlagBit(ISMRMRD::ISMRMRD_ACQ_IS_PARALLEL_CALIBRATION_AND_IMAGING).isSet(acq_header->flags))
        {
            if(this->imaging_data_arrived_[espace][slc])
            {
                // forward this ref data
                if (this->next()->putq(m1) == -1)
                {
                    GERROR("GenericReconReferenceKSpaceDelayedBufferGadget::process, passing data on to next gadget");
                    return GADGET_FAIL;
                }

                return GADGET_OK;
            }
            else
            {
                // buffer ref data
                ref_buf_[espace][slc].push_back(m1);
            }
        }
        else 
        {
            // if an imaging line arrives
            if (!(
                ISMRMRD::FlagBit(ISMRMRD::ISMRMRD_ACQ_IS_PARALLEL_CALIBRATION).isSet(m1->getObjectPtr()->flags) ||
                ISMRMRD::FlagBit(ISMRMRD::ISMRMRD_ACQ_IS_PHASECORR_DATA).isSet(m1->getObjectPtr()->flags)
                ))
            {
                this->imaging_data_arrived_[espace][slc] = true;

                // send out the ref buffer
                size_t num_acq = ref_buf_[espace][slc].size();
                for (size_t n=0; n<num_acq; n++)
                {
                    GadgetContainerMessage< ISMRMRD::AcquisitionHeader >* tmpm1 = ref_buf_[espace][slc][n];
                    if (tmpm1 != NULL)
                    {
                        if (this->next()->putq(tmpm1) == -1)
                        {
                            GERROR("GenericReconReferenceKSpaceDelayedBufferGadget::process, passing ref buffer line on to next gadget");
                            tmpm1->release();
                            return GADGET_FAIL;
                        }
                    }

                    ref_buf_[espace][slc][n] = NULL;
                }

                ref_buf_[espace][slc].clear();

                // send out this imaging line
                if (this->next()->putq(m1) == -1)
                {
                    GERROR("GenericReconReferenceKSpaceDelayedBufferGadget::process, passing imaging line on to next gadget");
                    return GADGET_FAIL;
                }

                return GADGET_OK;
            }
        }

        return GADGET_OK;
    }

    int GenericReconReferenceKSpaceDelayedBufferGadget::close(unsigned long flags)
    {
        GDEBUG_CONDITION_STREAM(true, "GenericReconReferenceKSpaceDelayedBufferGadget - close(flags) : " << flags);

        if (BaseClass::close(flags) != GADGET_OK) return GADGET_FAIL;

        if (flags != 0)
        {
            size_t e, slc;
            for (e = 0; e < ref_buf_.size(); e++)
            {
                if (calib_mode_[e] != Gadgetron::ISMRMRD_separate)
                    continue;

                for (slc = 0; slc < SLC_[e]; slc++)
                {
                    size_t num_acq = ref_buf_[e][slc].size();
                    for (size_t n = 0; n<num_acq; n++)
                    {
                        GadgetContainerMessage< ISMRMRD::AcquisitionHeader >* tmpm1 = ref_buf_[e][slc][n];
                        if (tmpm1 != NULL)
                        {
                            if (this->next()->putq(tmpm1) == -1)
                            {
                                GERROR("GenericReconReferenceKSpaceDelayedBufferGadget::process, passing ref buffer line on to next gadget in close(...) ... ");
                                tmpm1->release();
                                return GADGET_FAIL;
                            }
                        }
                    }

                    ref_buf_[e][slc].clear();
                }
            }
        }

        return GADGET_OK;
    }

    GADGET_FACTORY_DECLARE(GenericReconReferenceKSpaceDelayedBufferGadget)
}
