#include "GtPlusAccumulatorWorkOrderTriggerGadget.h"
#include "GtPlusReconGadgetUtil.h"

using namespace Gadgetron::gtPlus;

namespace Gadgetron
{

GtPlusAccumulatorWorkOrderTriggerGadget::GtPlusAccumulatorWorkOrderTriggerGadget() : 
                                            image_counter_(0), image_series_(100), first_kspace_scan_(true), 
                                            triggered_in_close_(false), triggered_in_process_(false), triggered_in_process_last_acq_(false), 
                                            triggered_in_process_by_numOfKSpace_triggerDim1_(false), 
                                            prev_dim1_(-1), curr_dim1_(-1), 
                                            prev_dim2_(-1), curr_dim2_(-1), 
                                            count_dim1_(0), 
                                            last_acq_arrived_(false), 
                                            verboseMode_(false), 
                                            other_kspace_matching_Dim_(DIM_NONE)
{
    space_matrix_offset_E1_ = 0;
    space_matrix_offset_E2_ = 0;

    gtPlusISMRMRDReconUtil<ValueType>().clearAcquisitionHeaderISMRMRD(prev_acq_header_);
    memset(&meas_max_idx_ref_, 0, sizeof(ISMRMRD::EncodingCounters));

    ind_time_stamp_.resize(GT_DIM_NUM, 0);

    embedded_ref_lines_E1_ = 0;
    embedded_ref_lines_E2_ = 0;

    min_sampled_E1_ = 0;
    min_sampled_E2_ = 0;

    center_line_E1_ = 0;
    center_line_E2_ = 0;

    max_sampled_E1_ = 0;
    max_sampled_E2_ = 0;

    timeStampResolution_ = 0.0025f;
}

GtPlusAccumulatorWorkOrderTriggerGadget::~GtPlusAccumulatorWorkOrderTriggerGadget()
{

}

// extract necessary configuration information from the xml
int GtPlusAccumulatorWorkOrderTriggerGadget::process_config(ACE_Message_Block* mb)
{
    // gadget parameters
    image_series_ = image_series.value();

    noacceleration_triggerDim1_ = gtPlus_util_.getISMRMRDDimFromName(noacceleration_triggerDim1.value());
    noacceleration_triggerDim2_ = gtPlus_util_.getISMRMRDDimFromName(noacceleration_triggerDim2.value());
    noacceleration_numOfKSpace_triggerDim1_ = noacceleration_numOfKSpace_triggerDim1.value();

    interleaved_triggerDim1_ = gtPlus_util_.getISMRMRDDimFromName(interleaved_triggerDim1.value());
    interleaved_triggerDim2_ = gtPlus_util_.getISMRMRDDimFromName(interleaved_triggerDim2.value());
    interleaved_numOfKSpace_triggerDim1_ = interleaved_numOfKSpace_triggerDim1.value();

    embedded_triggerDim1_ = gtPlus_util_.getISMRMRDDimFromName(embedded_triggerDim1.value());
    embedded_triggerDim2_ = gtPlus_util_.getISMRMRDDimFromName(embedded_triggerDim2.value());
    embedded_numOfKSpace_triggerDim1_ = embedded_numOfKSpace_triggerDim1.value();

    separate_triggerDim1_ = gtPlus_util_.getISMRMRDDimFromName(separate_triggerDim1.value());
    separate_triggerDim2_ = gtPlus_util_.getISMRMRDDimFromName(separate_triggerDim2.value());
    separate_numOfKSpace_triggerDim1_ = separate_numOfKSpace_triggerDim1.value();

    other_kspace_matching_Dim_ = gtPlus_util_.getISMRMRDDimFromName(other_kspace_matching_Dim.value());

    verboseMode_ = verboseMode.value();

    timeStampResolution_ = timeStampResolution.value();
    if ( timeStampResolution_ < FLT_EPSILON ) timeStampResolution_ = 0.0025f;
    GDEBUG_CONDITION_STREAM(verboseMode_, "timeStampResolution_ is " << timeStampResolution_);

    // ---------------------------------------------------------------------------------------------------------
    // pass the xml file
    ISMRMRD::IsmrmrdHeader h;
    try {
      deserialize(mb->rd_ptr(),h);
    } catch (...) {
      GDEBUG("Error parsing ISMRMRD Header");
      throw;
      return GADGET_FAIL;
    }


    // This only supports two encoding spaces where the recon_space is the same size
    // e.g. Parallel imaging reference scan collected with GRE and data with EPI
    if (h.encoding.size() > 2)
    {
        GDEBUG("Number of encoding spaces: %d\n", h.encoding.size());
        GDEBUG("This GtPlusAccumulatorWorkOrderTriggerGadget only supports two encoding space\n");
        return GADGET_FAIL;
    } 
    else if (h.encoding.size() == 2)
    {
        if (! ((h.encoding[0].reconSpace.matrixSize.x == h.encoding[1].reconSpace.matrixSize.x) && 
            (h.encoding[0].reconSpace.matrixSize.y == h.encoding[1].reconSpace.matrixSize.y) && 
            (h.encoding[0].reconSpace.matrixSize.z == h.encoding[1].reconSpace.matrixSize.z) && 
            (h.encoding[0].reconSpace.fieldOfView_mm.x == h.encoding[1].reconSpace.fieldOfView_mm.x) &&
            (h.encoding[0].reconSpace.fieldOfView_mm.y == h.encoding[1].reconSpace.fieldOfView_mm.y) &&
            (h.encoding[0].reconSpace.fieldOfView_mm.z == h.encoding[1].reconSpace.fieldOfView_mm.z)) )
        {
            GDEBUG("Number of encoding spaces: %d\n", h.encoding.size());
            GDEBUG("This GtPlusAccumulatorWorkOrderTriggerGadget only supports two encoding spaces with identical recon spaces.\n");
            return GADGET_FAIL;
        }
    }

    // find out the PAT mode
    if (!h.encoding[0].parallelImaging)
    {
      GDEBUG("Parallel Imaging section not found in header");
      return GADGET_FAIL;
    }

    ISMRMRD::ParallelImaging p_imaging = *h.encoding[0].parallelImaging;

    workOrder_.acceFactorE1_ = (double)(p_imaging.accelerationFactor.kspace_encoding_step_1);
    workOrder_.acceFactorE2_ = (double)(p_imaging.accelerationFactor.kspace_encoding_step_2);

    GDEBUG_CONDITION_STREAM(verboseMode_, "acceFactorE1_ is " << workOrder_.acceFactorE1_);
    GDEBUG_CONDITION_STREAM(verboseMode_, "acceFactorE2_ is " << workOrder_.acceFactorE2_);

    workOrder_.InterleaveDim_ = Gadgetron::DIM_NONE;

    if ( !p_imaging.calibrationMode.is_present() )
    {
        GDEBUG("Parallel Imaging calibrationMode not found in header");
        return GADGET_FAIL;
    }

    std::string calib = *p_imaging.calibrationMode;
    if ( calib.compare("interleaved") == 0 )
    {
        workOrder_.CalibMode_ = Gadgetron::ISMRMRD_interleaved;
        GDEBUG_CONDITION_STREAM(verboseMode_, "Calibration mode is interleaved");

        if ( p_imaging.interleavingDimension )
        {
            if ( p_imaging.interleavingDimension->compare("phase") == 0 )
            {
                workOrder_.InterleaveDim_ = Gadgetron::DIM_Phase;
            }
            else if ( p_imaging.interleavingDimension->compare("repetition") == 0 )
            {
                workOrder_.InterleaveDim_ = Gadgetron::DIM_Repetition;
            }
            else if ( p_imaging.interleavingDimension->compare("average") == 0 )
            {
                workOrder_.InterleaveDim_ = Gadgetron::DIM_Average;
            }
            else if ( p_imaging.interleavingDimension->compare("contrast") == 0 )
            {
                workOrder_.InterleaveDim_ = Gadgetron::DIM_Contrast;
            }
            else if ( p_imaging.interleavingDimension->compare("other") == 0 )
            {
                workOrder_.InterleaveDim_ = Gadgetron::DIM_other1;
            }
            else
            {
                GDEBUG("Unknown interleaving dimension. Bailing out");
                return GADGET_FAIL;
            }
            GDEBUG_CONDITION_STREAM(verboseMode_, "InterleaveDim is " << gtPlus_util_.getISMRMRDDimName(workOrder_.InterleaveDim_));
        }
    }
    else if ( calib.compare("embedded") == 0 )
    {
        workOrder_.CalibMode_ = Gadgetron::ISMRMRD_embedded;
        GDEBUG_CONDITION_STREAM(verboseMode_, "Calibration mode is embedded");
    }
    else if ( calib.compare("separate") == 0 )
    {
        workOrder_.CalibMode_ = Gadgetron::ISMRMRD_separate;
        GDEBUG_CONDITION_STREAM(verboseMode_, "Calibration mode is separate");
    }
    else if ( calib.compare("external") == 0 )
    {
        workOrder_.CalibMode_ = Gadgetron::ISMRMRD_external;
    }
    else if ( (calib.compare("other") == 0) && workOrder_.acceFactorE1_==1 && workOrder_.acceFactorE2_==1 )
    {
        workOrder_.CalibMode_ = Gadgetron::ISMRMRD_noacceleration;
        workOrder_.acceFactorE1_=1;
    }
    else if ( (calib.compare("other") == 0) &&  (workOrder_.acceFactorE1_>1 || workOrder_.acceFactorE2_>1) )
    {
        workOrder_.CalibMode_ = Gadgetron::ISMRMRD_interleaved;
        workOrder_.acceFactorE1_=2;
        workOrder_.InterleaveDim_ = Gadgetron::DIM_Phase;
    }
    else
    {
        GDEBUG("Failed to process parallel imaging calibration mode");
        return GADGET_FAIL;
    }
    
    // ---------------------------------------------------------------------------------------------------------

    // find out the encoding space 

    findMatrixSizeEncoding(h, matrix_size_encoding_);
    findFOVEncoding(h, field_of_view_encoding_);

    findMatrixSizeRecon(h, matrix_size_recon_);
    findFOVRecon(h, field_of_view_recon_);

    GDEBUG_CONDITION_STREAM(verboseMode_, "Encoding matrix size: " << matrix_size_encoding_[0] << " " << matrix_size_encoding_[1] << " " << matrix_size_encoding_[2]);
    GDEBUG_CONDITION_STREAM(verboseMode_, "Encoding field_of_view : " << field_of_view_encoding_[0] << " " << field_of_view_encoding_[1] << " " << field_of_view_encoding_[2]);
    GDEBUG_CONDITION_STREAM(verboseMode_, "Recon matrix size : " << matrix_size_recon_[0] << " " << matrix_size_recon_[1] << " " << matrix_size_recon_[2]);
    GDEBUG_CONDITION_STREAM(verboseMode_, "Recon field_of_view :  " << field_of_view_recon_[0] << " " << field_of_view_recon_[1] << " " << field_of_view_recon_[2]);

    // ---------------------------------------------------------------------------------------------------------
    // handle partial fourier

    ISMRMRD::EncodingSpace e_space = h.encoding[0].encodedSpace;
    ISMRMRD::EncodingLimits e_limits = h.encoding[0].encodingLimits;

    workOrder_.kSpaceMaxEncode1_ = matrix_size_encoding_[1]-1;
    GDEBUG_CONDITION_STREAM(verboseMode_, "matrix size kSpaceMaxEncode1_ is " << workOrder_.kSpaceMaxEncode1_);

    workOrder_.kSpaceMaxEncode2_ = matrix_size_encoding_[2]-1;
    GDEBUG_CONDITION_STREAM(verboseMode_, "matrix size kSpaceMaxEncode2_ is " << workOrder_.kSpaceMaxEncode2_);

    space_size_[1] = workOrder_.kSpaceMaxEncode1_+1;
    space_size_[2] = workOrder_.kSpaceMaxEncode2_+1;

    if ( (!e_limits.kspace_encoding_step_1) || (!e_limits.kspace_encoding_step_2))
    {
        GDEBUG("kspace_encoding_step_1 and kspace_encoding_step_2 limits are required. Not found. Bailing out.");
        return GADGET_FAIL;
    }

    min_sampled_E1_ = e_limits.kspace_encoding_step_1->minimum;
    min_sampled_E2_ = e_limits.kspace_encoding_step_2->minimum;

    GDEBUG_CONDITION_STREAM(verboseMode_, "min_sampled_E1_ is " << min_sampled_E1_);
    GDEBUG_CONDITION_STREAM(verboseMode_, "min_sampled_E2_ is " << min_sampled_E2_);

    max_sampled_E1_ = e_limits.kspace_encoding_step_1->maximum;
    max_sampled_E2_ = e_limits.kspace_encoding_step_2->maximum;

    GDEBUG_CONDITION_STREAM(verboseMode_, "max_sampled_E1_ is " << max_sampled_E1_);
    GDEBUG_CONDITION_STREAM(verboseMode_, "max_sampled_E2_ is " << max_sampled_E2_);

    center_line_E1_ = e_limits.kspace_encoding_step_1->center;
    center_line_E2_ = e_limits.kspace_encoding_step_2->center;

    GDEBUG_CONDITION_STREAM(verboseMode_, "center_line_E1_ is " << center_line_E1_);
    GDEBUG_CONDITION_STREAM(verboseMode_, "center_line_E2_ is " << center_line_E2_);

    workOrder_.kSpaceCenterEncode1_ = center_line_E1_;
    GDEBUG_CONDITION_STREAM(verboseMode_, "kSpaceCenterEncode1_ is " << workOrder_.kSpaceCenterEncode1_);

    workOrder_.kSpaceCenterEncode2_ = center_line_E2_;
    GDEBUG_CONDITION_STREAM(verboseMode_, "kSpaceCenterEncode2_ is " << workOrder_.kSpaceCenterEncode2_);

    // ---------------------------------------------------------------------------------------------------------
    // handle retro-gating
    if (h.userParameters)
    {
        for (std::vector<ISMRMRD::UserParameterLong>::const_iterator  i = h.userParameters->userParameterLong.begin (); i != h.userParameters->userParameterLong.end(); ++i)
        {
            if (i->name == "RetroGatedImages")
            {
                workOrder_.retro_gated_images_ = i->value;
            }
            else if ( i->name == "RetroGatedSegmentSize")
            {
                workOrder_.retro_gated_segment_size_ = i->value;
            }
            else if ( i->name == "EmbeddedRefLinesE1")
            {
                embedded_ref_lines_E1_ = i->value;
            }
            else if ( i->name == "EmbeddedRefLinesE2")
            {
                embedded_ref_lines_E2_ = i->value;
            }
        }
    }

    // ---------------------------------------------------------------------------------------------------------
    // encoding limits

    if ( std::abs(2*field_of_view_recon_[0]-field_of_view_encoding_[0]) < 1.0 )
    {
        meas_max_ro_ = e_space.matrixSize.x/2;
    }
    else
    {
        meas_max_ro_ = e_space.matrixSize.x;
    }
    space_size_[0] = meas_max_ro_;

    meas_max_idx_.kspace_encode_step_1 = (uint16_t)matrix_size_encoding_[1]-1;

    meas_max_idx_.set = (e_limits.set && (e_limits.set->maximum>0)) ? e_limits.set->maximum : 0;
    meas_max_idx_.phase = (e_limits.phase && (e_limits.phase->maximum>0)) ? e_limits.phase->maximum : 0;

    // if it is retro-gating
    if ( workOrder_.retro_gated_images_ > 0 )
    {
        meas_max_idx_.phase = (uint16_t)(workOrder_.retro_gated_images_ - 1);
    }

    meas_max_idx_.kspace_encode_step_2 = (uint16_t)matrix_size_encoding_[2]-1;

    meas_max_idx_.contrast = (e_limits.contrast && (e_limits.contrast->maximum > 0)) ? e_limits.contrast->maximum : 0;

    meas_max_idx_.slice = (e_limits.slice && (e_limits.slice->maximum > 0)) ? e_limits.slice->maximum : 0;

    meas_max_idx_.repetition = e_limits.repetition ? e_limits.repetition->maximum : 0;

    meas_max_idx_.average = e_limits.average ? e_limits.average->maximum : 0;

    meas_max_idx_.segment = 0;

    return GADGET_OK;
}

int GtPlusAccumulatorWorkOrderTriggerGadget::process(GadgetContainerMessage<ISMRMRD::AcquisitionHeader>* m1, GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2)
{
    // logic to control whether to store kspace and ref data
    bool bIsKSpace, bIsRef, bIsNoise, bIsPhaseCorr, bIsReflect, bIsOther, bIsNavigator, bIsRTFeedback, bIsHPFeedback, bIsDummyScan;
    if ( !checkStatus(m1->getObjectPtr()->flags, m1->getObjectPtr()->number_of_samples, 
            bIsKSpace, bIsRef, bIsNoise, bIsPhaseCorr, bIsReflect, bIsOther,
            bIsNavigator, bIsRTFeedback, bIsHPFeedback, bIsDummyScan) )
    {
        GDEBUG("Failed check readout status\n");
        return GADGET_FAIL;
    }

    size_t scan_counter = m1->getObjectPtr()->scan_counter;

    if ( scan_counter%1000 == 0 )
    {
        GDEBUG_CONDITION_STREAM(verboseMode_, "--> receive scan : " << scan_counter);
    }

    // combine the segmentes
    //if ( workOrder_.retro_gated_images_ == 0 )
    //{
        m1->getObjectPtr()->idx.segment = 0;
    //}

    if ( (bIsNavigator || bIsRTFeedback || bIsHPFeedback || bIsDummyScan) && !bIsKSpace && !bIsRef )
    {
        m1->release();
        return GADGET_OK;
    }

    if ( !bIsRTFeedback && bIsKSpace && first_kspace_scan_ && m1->getObjectPtr()->center_sample>0 )
    {
        if ( (workOrder_.start_RO_<0) && (workOrder_.end_RO_<0) )
        {
	    workOrder_.start_RO_ = m1->getObjectPtr()->discard_pre;
	    workOrder_.end_RO_ = m1->getObjectPtr()->number_of_samples - m1->getObjectPtr()->discard_post - 1;

            GDEBUG_CONDITION_STREAM(verboseMode_, "start_RO : " << workOrder_.start_RO_);
            GDEBUG_CONDITION_STREAM(verboseMode_, "end_RO : " << workOrder_.end_RO_);

            workOrder_.kSpaceCenterRO_ = m1->getObjectPtr()->center_sample;
            workOrder_.kSpaceMaxRO_ = m1->getObjectPtr()->number_of_samples;
        }

        // if partial fourier or asymmetric echo is used, correct the kSpaceCenter
        //if ( std::abs( (long long)(space_size_[1])-(long long)max_sampled_E1_) > workOrder_.acceFactorE1_ )
        //{
        //    GDEBUG_CONDITION_STREAM(verboseMode_, "Partial fourier along E1 ... ");

        //    // if ( (m1->getObjectPtr()->idx.user[5]>0) && (std::abs( (long long)m1->getObjectPtr()->idx.user[5] - (long long)space_size_[1]/2 )<2) )
        //    if ( (m1->getObjectPtr()->idx.user[5]>0) )
        //    {
        //        workOrder_.kSpaceCenterEncode1_ = m1->getObjectPtr()->idx.user[5];
        //    }

        //    if ( 2*workOrder_.kSpaceCenterEncode1_ >= (max_sampled_E1_+1) )
        //    {
        //        space_matrix_offset_E1_ = 0;

        //        workOrder_.start_E1_ = 0;
        //        workOrder_.end_E1_ = (int)max_sampled_E1_;
        //    }
        //    else
        //    {
        //        space_matrix_offset_E1_ = space_size_[1] - max_sampled_E1_ -1;

        //        workOrder_.start_E1_ = (int)space_matrix_offset_E1_;
        //        workOrder_.end_E1_ = (int)workOrder_.kSpaceMaxEncode1_;
        //    }

        //    workOrder_.kSpaceMaxEncode1_ = 2*workOrder_.kSpaceCenterEncode1_-1;
        //}
        //else
        //{
        //    space_matrix_offset_E1_ = 0;
        //}

        //if ( std::abs( (long long)space_size_[2] - (long long)max_sampled_E2_) > workOrder_.acceFactorE2_ )
        //{
        //    GDEBUG_CONDITION_STREAM(verboseMode_, "Partial fourier along E2 ... ");

        //    // if ( (m1->getObjectPtr()->idx.user[6]>0) && (std::abs( (long long)m1->getObjectPtr()->idx.user[6] - (long long)space_size_[2]/2 )<2) )
        //    if ( (m1->getObjectPtr()->idx.user[6]>0) )
        //    {
        //        workOrder_.kSpaceCenterEncode2_ = m1->getObjectPtr()->idx.user[6];
        //    }

        //    if ( 2*workOrder_.kSpaceCenterEncode2_ >= (max_sampled_E2_+1) )
        //    {
        //        space_matrix_offset_E2_ = 0;

        //        workOrder_.start_E2_ = 0;
        //        workOrder_.end_E2_ = (int)max_sampled_E2_;
        //    }
        //    else
        //    {
        //        space_matrix_offset_E2_ = space_size_[2] - max_sampled_E2_-1;

        //        workOrder_.start_E2_ = (int)space_matrix_offset_E2_;
        //        workOrder_.end_E2_ = (int)workOrder_.kSpaceMaxEncode2_;
        //    }

        //    workOrder_.kSpaceMaxEncode2_ = 2*workOrder_.kSpaceCenterEncode2_-1;
        //}
        //else
        //{
        //    space_matrix_offset_E2_ = 0;
        //}

        space_matrix_offset_E1_ = (matrix_size_encoding_[1] / 2 - center_line_E1_);
        GDEBUG_CONDITION_STREAM(verboseMode_, "Center line offset along E1 : " << space_matrix_offset_E1_);

        if (std::abs((long long)center_line_E1_ - (long long)((max_sampled_E1_ + 1) / 2)) > workOrder_.acceFactorE1_)
        {
            GDEBUG_CONDITION_STREAM(verboseMode_, "Partial fourier along E1 ... ");

            workOrder_.start_E1_ = (int)(min_sampled_E1_ + space_matrix_offset_E1_);
            workOrder_.end_E1_ = (int)(max_sampled_E1_ + space_matrix_offset_E1_);

            GADGET_CHECK_RETURN_FALSE(workOrder_.start_E1_ >= 0);
            GADGET_CHECK_RETURN_FALSE(workOrder_.start_E1_ < matrix_size_encoding_[1]);

            GADGET_CHECK_RETURN_FALSE(workOrder_.end_E1_ >= 0);
            GADGET_CHECK_RETURN_FALSE(workOrder_.end_E1_ < matrix_size_encoding_[1]);
        }

        // -------------------------------------------------------------------------

        space_matrix_offset_E2_ = (matrix_size_encoding_[2] / 2 - center_line_E2_);
        GDEBUG_CONDITION_STREAM(verboseMode_, "Center line offset along E2 : " << space_matrix_offset_E2_);

        if (std::abs((long long)center_line_E2_ - (long long)((max_sampled_E2_ + 1) / 2)) > workOrder_.acceFactorE2_)
        {
            GDEBUG_CONDITION_STREAM(verboseMode_, "Partial fourier along E2 ... ");

            workOrder_.start_E2_ = (int)(min_sampled_E2_ + space_matrix_offset_E2_);
            workOrder_.end_E2_ = (int)(max_sampled_E2_ + space_matrix_offset_E2_);

            GADGET_CHECK_RETURN_FALSE(workOrder_.start_E2_ >= 0);
            GADGET_CHECK_RETURN_FALSE(workOrder_.start_E2_ < matrix_size_encoding_[2]);

            GADGET_CHECK_RETURN_FALSE(workOrder_.end_E2_ >= 0);
            GADGET_CHECK_RETURN_FALSE(workOrder_.end_E2_ < matrix_size_encoding_[2]);
        }

        first_kspace_scan_ = false;
    }

    // store kspace read out
    if ( bIsKSpace )
    {
        if ( !storeImageData(m1, m2, bIsReflect) )
        {
            GDEBUG("Failed check readout status\n");
            return GADGET_FAIL;
        }
    }

    // store ref read out
    if ( bIsRef && (workOrder_.CalibMode_ != Gadgetron::ISMRMRD_interleaved) )
    {
        if ( !storeRefData(m1, m2, bIsReflect) )
        {
            GDEBUG("Failed check readout status\n");
            return GADGET_FAIL;
        }
    }

    // store phaseCorr read out
    if ( bIsPhaseCorr )
    {
        ISMRMRD::AcquisitionHeader* pMDH = m1->getObjectPtr();
        hoNDArray< ValueType >* pRefLine = m2->getObjectPtr();

        ReadOutBuffer item;
        item.acqHead_ = *pMDH;
        item.data_ = *pRefLine;
        item.isReflect_ = bIsReflect;
        phaseCorrBuffer_.push_back(item);
    }

    // store noise read out
    if ( bIsNoise )
    {
        ISMRMRD::AcquisitionHeader* pMDH = m1->getObjectPtr();
        hoNDArray< ValueType >* pRefLine = m2->getObjectPtr();

        ReadOutBuffer item;
        item.acqHead_ = *pMDH;
        item.data_ = *pRefLine;
        item.isReflect_ = bIsReflect;
        noiseBuffer_.push_back(item);
    }

    // store other read out
    if ( bIsOther )
    {
        ISMRMRD::AcquisitionHeader* pMDH = m1->getObjectPtr();
        hoNDArray< ValueType >* pRefLine = m2->getObjectPtr();

        if ( other_kspace_matching_Dim_ != DIM_NONE )
        {
            if ( prev_acq_header_.measurement_uid != 0 )
            {
                size_t v = getDimValue(prev_acq_header_, other_kspace_matching_Dim_);
                setDimValue(*pMDH, other_kspace_matching_Dim_, v+1);
            }
        }

        ReadOutBuffer item;
        item.acqHead_ = *pMDH;
        item.data_ = *pRefLine;
        item.isReflect_ = bIsReflect;
        otherBuffer_.push_back(item);
    }

    // perform triggering
    if ( !triggerWorkOrder(m1, false, bIsKSpace) )
    {
        GDEBUG("Failed triggerWorkOrder(m1)\n");
        return GADGET_FAIL;
    }

    m1->release();
    return GADGET_OK;
}

bool GtPlusAccumulatorWorkOrderTriggerGadget::needTriggerWorkOrderAllInClose()
{
    // already triggered for last acquisition
    if ( triggered_in_process_last_acq_ ) return false;

    // if never triggered in process(...) and the last acqusition does arrive
    // if the last acquisition does not arrive, the user may has cancel the scan
    if ( !triggered_in_process_ && !triggered_in_process_last_acq_ && last_acq_arrived_ ) return true;

    if ( workOrder_.CalibMode_ == ISMRMRD_interleaved )
    {
        return ((interleaved_triggerDim1_==DIM_NONE)&&(interleaved_triggerDim2_==DIM_NONE));
    }
    else if ( workOrder_.CalibMode_ == ISMRMRD_embedded )
    {
        return ((embedded_triggerDim1_==DIM_NONE)&&(embedded_triggerDim2_==DIM_NONE));
    }
    else if ( (workOrder_.CalibMode_ == ISMRMRD_separate) 
            || (workOrder_.CalibMode_ == ISMRMRD_external) )
    {
        return ((separate_triggerDim1_==DIM_NONE)&&(separate_triggerDim2_==DIM_NONE));
    }
    else if ( workOrder_.CalibMode_ == ISMRMRD_noacceleration )
    {
        return ((noacceleration_triggerDim1_==DIM_NONE)&&(noacceleration_triggerDim2_==DIM_NONE));
    }
    else
    {
        GERROR_STREAM("Unsupported calibration mode : " << workOrder_.CalibMode_);
        return true;
    }

    return true;
}

bool GtPlusAccumulatorWorkOrderTriggerGadget::triggerWorkOrder(GadgetContainerMessage<ISMRMRD::AcquisitionHeader>* m1, bool inClose, bool isKSpace)
{
    if ( workOrder_.CalibMode_ == ISMRMRD_interleaved )
    {
        if ( inClose )
        {
            GADGET_CHECK_RETURN_FALSE(triggerWorkOrderLastCountInClose(interleaved_triggerDim1_, interleaved_triggerDim2_, interleaved_numOfKSpace_triggerDim1_));
        }
        else
        {
            if ( isKSpace )
            {
                GADGET_CHECK_RETURN_FALSE(triggerWorkOrder(m1, interleaved_triggerDim1_, interleaved_triggerDim2_, interleaved_numOfKSpace_triggerDim1_));
            }
        }
    }
    else if ( workOrder_.CalibMode_ == ISMRMRD_embedded )
    {
        if ( inClose )
        {
            GADGET_CHECK_RETURN_FALSE(triggerWorkOrderLastCountInClose(embedded_triggerDim1_, embedded_triggerDim2_, embedded_numOfKSpace_triggerDim1_));
        }
        else
        {
            if ( isKSpace )
            {
                GADGET_CHECK_RETURN_FALSE(triggerWorkOrder(m1, embedded_triggerDim1_, embedded_triggerDim2_, embedded_numOfKSpace_triggerDim1_));
            }
        }
    }
    else if ( (workOrder_.CalibMode_ == ISMRMRD_separate) 
            || (workOrder_.CalibMode_ == ISMRMRD_external) )
    {
        if ( inClose )
        {
            GADGET_CHECK_RETURN_FALSE(triggerWorkOrderLastCountInClose(separate_triggerDim1_, separate_triggerDim2_, separate_numOfKSpace_triggerDim1_));
        }
        else
        {
            if ( isKSpace )
            {
                GADGET_CHECK_RETURN_FALSE(triggerWorkOrder(m1, separate_triggerDim1_, separate_triggerDim2_, separate_numOfKSpace_triggerDim1_));
            }
        }
    }
    else if ( workOrder_.CalibMode_ == ISMRMRD_noacceleration )
    {
        if ( inClose )
        {
            GADGET_CHECK_RETURN_FALSE(triggerWorkOrderLastCountInClose(noacceleration_triggerDim1_, noacceleration_triggerDim2_, noacceleration_numOfKSpace_triggerDim1_));
        }
        else
        {
            if ( isKSpace )
            {
                GADGET_CHECK_RETURN_FALSE(triggerWorkOrder(m1, noacceleration_triggerDim1_, noacceleration_triggerDim2_, noacceleration_numOfKSpace_triggerDim1_));
            }
        }
    }
    else
    {
        GERROR_STREAM("Unsupported calibration mode : " << workOrder_.CalibMode_);
        return false;
    }

    return true;
}

bool GtPlusAccumulatorWorkOrderTriggerGadget::
resetTriggerStatus(GadgetContainerMessage< ISMRMRD::AcquisitionHeader >* m1)
{
    // return !gtPlusISMRMRDReconUtil<ValueType>().hasIdenticalGeometryISMRMRD(*(m1->getObjectPtr()), prev_acq_header_);
    return false;
}

bool GtPlusAccumulatorWorkOrderTriggerGadget::
triggerWorkOrder(GadgetContainerMessage<ISMRMRD::AcquisitionHeader>* m1,
            Gadgetron::ISMRMRDDIM& triggerDim1_, 
            Gadgetron::ISMRMRDDIM& triggerDim2_,
            int numOfKSpace_triggerDim1_)
{
    //bool is_first_acq_in_slice = ISMRMRD::FlagBit(ISMRMRD::ISMRMRD_ACQ_FIRST_IN_SLICE).isSet(m1->getObjectPtr()->flags);
    //if ( !is_first_acq_in_slice ) return true;

    bool is_last_acq = ( ((ISMRMRD::FlagBit(ISMRMRD::ISMRMRD_ACQ_LAST_IN_REPETITION).isSet(m1->getObjectPtr()->flags)) || (ISMRMRD::FlagBit(ISMRMRD::ISMRMRD_ACQ_LAST_IN_SLICE).isSet(m1->getObjectPtr()->flags)) ) 
                                && (m1->getObjectPtr()->idx.repetition==meas_max_idx_.repetition)
                                && (m1->getObjectPtr()->idx.slice==meas_max_idx_.slice)
                                && (m1->getObjectPtr()->idx.set==meas_max_idx_.set)
                                && (m1->getObjectPtr()->idx.contrast==meas_max_idx_.contrast)
                                && (m1->getObjectPtr()->idx.phase==meas_max_idx_.phase)
                                && (m1->getObjectPtr()->idx.average==meas_max_idx_.average) );

    // if retro gating, use the end of acq flag
    if ( !is_last_acq && (workOrder_.retro_gated_images_ > 0) )
    {
        is_last_acq = (ISMRMRD::FlagBit(ISMRMRD::ISMRMRD_ACQ_LAST_IN_MEASUREMENT).isSet(m1->getObjectPtr()->flags));
    }

    if ( is_last_acq ) last_acq_arrived_ = true;

    curr_dim1_ = getDimValue(*(m1->getObjectPtr()), triggerDim1_);
    curr_dim2_ = getDimValue(*(m1->getObjectPtr()), triggerDim2_);

    if ( is_last_acq && ( (triggerDim1_!=DIM_NONE) || (triggerDim2_!=DIM_NONE) ) )
    {
        GDEBUG_CONDITION_STREAM(true, "Last scan in measurement - " << gtPlusISMRMRDReconUtil<ValueType>().getISMRMRDDimName(triggerDim1_ ) << " = " << curr_dim1_ << " - " << gtPlusISMRMRDReconUtil<ValueType>().getISMRMRDDimName(triggerDim2_ ) << " = " << curr_dim2_);

        if ( curr_dim1_==0 && curr_dim2_== 0 )
        {
            GDEBUG_CONDITION_STREAM(true, "Last scan in measurement - not trigger ... ");
            return true;
        }

        triggered_in_process_last_acq_ = true;
        GDEBUG_CONDITION_STREAM(true, "Last scan in measurement - triggered_in_process_last_acq_ : " << triggered_in_process_last_acq_);

        if ( workOrder_.CalibMode_ == ISMRMRD_interleaved )
        {
            GADGET_CHECK_RETURN_FALSE(triggerWorkOrderLastCountInClose(interleaved_triggerDim1_, interleaved_triggerDim2_, interleaved_numOfKSpace_triggerDim1_));
        }
        else if ( workOrder_.CalibMode_ == ISMRMRD_embedded )
        {
            GADGET_CHECK_RETURN_FALSE(triggerWorkOrderLastCountInClose(embedded_triggerDim1_, embedded_triggerDim2_, embedded_numOfKSpace_triggerDim1_));
        }
        else if ( (workOrder_.CalibMode_ == ISMRMRD_separate) 
                || (workOrder_.CalibMode_ == ISMRMRD_external) )
        {
            GADGET_CHECK_RETURN_FALSE(triggerWorkOrderLastCountInClose(separate_triggerDim1_, separate_triggerDim2_, separate_numOfKSpace_triggerDim1_));
        }
        else if ( workOrder_.CalibMode_ == ISMRMRD_noacceleration )
        {
            GADGET_CHECK_RETURN_FALSE(triggerWorkOrderLastCountInClose(noacceleration_triggerDim1_, noacceleration_triggerDim2_, noacceleration_numOfKSpace_triggerDim1_));
        }
        else
        {
            triggered_in_process_last_acq_ = false;
            GERROR_STREAM("Unsupported calibration mode : " << workOrder_.CalibMode_);
            return false;
        }

        return true;
    }

    if ( prev_dim1_ == -1 )
    {
        prev_dim1_ = curr_dim1_;
        count_dim1_ = 0;
        GDEBUG_CONDITION_STREAM(verboseMode_, "Current Dim1 : " << gtPlusISMRMRDReconUtil<ValueType>().getISMRMRDDimName(triggerDim1_ ) << " = " << curr_dim1_);
    }

    if ( prev_dim2_ == -1 )
    {
        prev_dim2_ = curr_dim2_;
        count_dim1_ = 0;
        GDEBUG_CONDITION_STREAM(verboseMode_, "Current Dim2 : " << gtPlusISMRMRDReconUtil<ValueType>().getISMRMRDDimName(triggerDim2_ ) << " = " << curr_dim2_);
    }

    if ( prev_acq_header_.measurement_uid == 0 ) prev_acq_header_ = *(m1->getObjectPtr());

    bool workFlow_BufferKernel_ = false;
    bool workFlow_use_BufferedKernel_ = false;

    if ( prev_dim1_ != curr_dim1_ )
    {
        count_dim1_++;
        GDEBUG_CONDITION_STREAM(verboseMode_, "Current Dim1 : " << gtPlusISMRMRDReconUtil<ValueType>().getISMRMRDDimName(triggerDim1_ ) << " = " << curr_dim1_);
        GDEBUG_CONDITION_STREAM(verboseMode_, "Current Dim2 : " << gtPlusISMRMRDReconUtil<ValueType>().getISMRMRDDimName(triggerDim2_ ) << " = " << curr_dim2_);
        GDEBUG_CONDITION_STREAM(verboseMode_, "count_dim1_ : " << count_dim1_);
    }

    if ( (triggerDim1_==DIM_NONE) && (triggerDim2_==DIM_NONE) )
    {
        prev_dim1_ = curr_dim1_;
        prev_dim2_ = curr_dim2_;
        prev_acq_header_ = *(m1->getObjectPtr());
        return true;
    }

    int numOfAcquiredKSpaceForTriggerDim1 = numOfKSpace_triggerDim1_;
    if ( workOrder_.CalibMode_ == ISMRMRD_interleaved )
    {
        numOfAcquiredKSpaceForTriggerDim1 = (int)(numOfKSpace_triggerDim1_ * workOrder_.acceFactorE1_ * workOrder_.acceFactorE2_);
    }

    // trigger whenever the Dim2 is changed
    if (  triggerDim1_==DIM_NONE && triggerDim2_!=DIM_NONE  )
    {
        prev_dim1_ = curr_dim1_;
        prev_acq_header_ = *(m1->getObjectPtr());

        size_t prev_dim2_local_ = prev_dim2_;
        prev_dim2_ = curr_dim2_;

        if ( curr_dim2_!= prev_dim2_local_ )
        {
            GDEBUG_CONDITION_STREAM(verboseMode_, "Trigger Dim2 : " << gtPlusISMRMRDReconUtil<ValueType>().getISMRMRDDimName(triggerDim2_ ) << " = " << prev_dim2_local_);
            GADGET_CHECK_RETURN_FALSE(triggerByDimEqual(triggerDim2_, prev_dim2_local_, workFlow_BufferKernel_, workFlow_use_BufferedKernel_));
            triggered_in_process_ = true;
        }
    }

    // trigger whenever the Dim1 is changed
    if (  triggerDim1_!=DIM_NONE && triggerDim2_==DIM_NONE  )
    {
        prev_dim2_ = curr_dim2_;

        size_t prev_dim1_local_ = prev_dim1_;
        prev_dim1_ = curr_dim1_;

        if ( numOfKSpace_triggerDim1_ > 0 )
        {
            if ( curr_dim1_!= prev_dim1_local_ )
            {
                if ( resetTriggerStatus(m1) )
                {
                    count_dim1_ = 0;
                    GDEBUG_CONDITION_STREAM(verboseMode_, "Trigger Dim1 : " << gtPlusISMRMRDReconUtil<ValueType>().getISMRMRDDimName(triggerDim1_ ) << " = " << prev_dim1_local_);

                    workFlow_BufferKernel_ = false;
                    workFlow_use_BufferedKernel_ = true;
                    GADGET_CHECK_RETURN_FALSE(triggerByDimEqual(triggerDim1_, prev_dim1_local_, workFlow_BufferKernel_, workFlow_use_BufferedKernel_));
                    triggered_in_process_ = true;
                }

                if ( count_dim1_ == numOfAcquiredKSpaceForTriggerDim1 )
                {
                    GDEBUG_CONDITION_STREAM(verboseMode_, "Trigger Dim1 : " << gtPlusISMRMRDReconUtil<ValueType>().getISMRMRDDimName(triggerDim1_ ) << " = " << prev_dim1_local_);

                    workFlow_BufferKernel_ = true;
                    workFlow_use_BufferedKernel_ = false;
                    GADGET_CHECK_RETURN_FALSE(triggerByDimLessEqual(triggerDim1_, prev_dim1_local_, workFlow_BufferKernel_, workFlow_use_BufferedKernel_));
                    triggered_in_process_ = true;
                }
                else if ( count_dim1_ > numOfAcquiredKSpaceForTriggerDim1 )
                {
                    GDEBUG_CONDITION_STREAM(verboseMode_, "Trigger Dim1 : " << gtPlusISMRMRDReconUtil<ValueType>().getISMRMRDDimName(triggerDim1_ ) << " = " << prev_dim1_local_);

                    workFlow_BufferKernel_ = false;
                    workFlow_use_BufferedKernel_ = true;
                    GADGET_CHECK_RETURN_FALSE(triggerByDimEqual(triggerDim1_, prev_dim1_local_, workFlow_BufferKernel_, workFlow_use_BufferedKernel_));
                    triggered_in_process_ = true;
                }
            }

            prev_acq_header_ = *(m1->getObjectPtr());
        }
        else
        {
            if ( curr_dim1_!= prev_dim1_local_ )
            {
                GDEBUG_CONDITION_STREAM(verboseMode_, "Trigger Dim1 : " << gtPlusISMRMRDReconUtil<ValueType>().getISMRMRDDimName(triggerDim1_ ) << " = " << prev_dim1_local_);
                GADGET_CHECK_RETURN_FALSE(triggerByDimEqual(triggerDim1_, prev_dim1_local_, workFlow_BufferKernel_, workFlow_use_BufferedKernel_));
                triggered_in_process_ = true;
            }

            prev_acq_header_ = *(m1->getObjectPtr());
        }
    }

    if (  triggerDim1_!=DIM_NONE && triggerDim2_!=DIM_NONE  )
    {
        size_t prev_dim1_local_ = prev_dim1_;
        size_t prev_dim2_local_ = prev_dim2_;

        prev_dim1_ = curr_dim1_;
        prev_dim2_ = curr_dim2_;

        if ( numOfKSpace_triggerDim1_ > 0 )
        {
            if ( (curr_dim2_!=prev_dim2_local_) || resetTriggerStatus(m1) )
            {
                if ( count_dim1_ > numOfAcquiredKSpaceForTriggerDim1 )
                {
                    count_dim1_ = 0;
                    GDEBUG_CONDITION_STREAM(verboseMode_, "Trigger Dim1 : " << gtPlusISMRMRDReconUtil<ValueType>().getISMRMRDDimName(triggerDim1_ ) << " = " << prev_dim1_local_ 
                        << "; Dim2 : " << gtPlusISMRMRDReconUtil<ValueType>().getISMRMRDDimName(triggerDim2_ ) << " = " << prev_dim2_local_);

                    workFlow_BufferKernel_ = false;
                    workFlow_use_BufferedKernel_ = true;

                    GADGET_CHECK_RETURN_FALSE(triggerByDimEqual(triggerDim1_, prev_dim1_local_, triggerDim2_, prev_dim2_local_, workFlow_BufferKernel_, workFlow_use_BufferedKernel_));
                }

                if ( count_dim1_ <= numOfAcquiredKSpaceForTriggerDim1 && !triggered_in_process_by_numOfKSpace_triggerDim1_ ) // the trigger never happened
                {
                    count_dim1_ = 0;
                    GDEBUG_CONDITION_STREAM(verboseMode_, "Trigger Dim1 : " << gtPlusISMRMRDReconUtil<ValueType>().getISMRMRDDimName(triggerDim1_ ) << " = " << prev_dim1_local_ 
                        << "; Dim2 : " << gtPlusISMRMRDReconUtil<ValueType>().getISMRMRDDimName(triggerDim2_ ) << " = " << prev_dim2_local_);

                    workFlow_BufferKernel_ = false;
                    workFlow_use_BufferedKernel_ = false;

                    GADGET_CHECK_RETURN_FALSE(triggerByDim1LessEqualDim2Equal(triggerDim1_, prev_dim1_local_, triggerDim2_, prev_dim2_local_, workFlow_BufferKernel_, workFlow_use_BufferedKernel_));
                }

                triggered_in_process_ = true;
                triggered_in_process_by_numOfKSpace_triggerDim1_ = false; // reset this flag to be false for next dim2
            }

            if (curr_dim1_!=prev_dim1_local_)
            {
                if ( count_dim1_ == numOfAcquiredKSpaceForTriggerDim1 )
                {
                    GDEBUG_CONDITION_STREAM(verboseMode_, "Trigger Dim1 : " << gtPlusISMRMRDReconUtil<ValueType>().getISMRMRDDimName(triggerDim1_ ) << " = " << prev_dim1_local_ 
                        << "; Dim2 : " << gtPlusISMRMRDReconUtil<ValueType>().getISMRMRDDimName(triggerDim2_ ) << " = " << prev_dim2_local_);

                    workFlow_BufferKernel_ = true;
                    workFlow_use_BufferedKernel_ = false;
                    GADGET_CHECK_RETURN_FALSE(triggerByDim1LessEqualDim2Equal(triggerDim1_, prev_dim1_local_, triggerDim2_, prev_dim2_local_, workFlow_BufferKernel_, workFlow_use_BufferedKernel_));
                    triggered_in_process_ = true;
                    triggered_in_process_by_numOfKSpace_triggerDim1_ = true;
                }
                else if ( count_dim1_ > numOfAcquiredKSpaceForTriggerDim1 )
                {
                    GDEBUG_CONDITION_STREAM(verboseMode_, "Trigger Dim1 : " << gtPlusISMRMRDReconUtil<ValueType>().getISMRMRDDimName(triggerDim1_ ) << " = " << prev_dim1_local_ 
                        << "; Dim2 : " << gtPlusISMRMRDReconUtil<ValueType>().getISMRMRDDimName(triggerDim2_ ) << " = " << prev_dim2_local_);

                    workFlow_BufferKernel_ = false;
                    workFlow_use_BufferedKernel_ = true;
                    GADGET_CHECK_RETURN_FALSE(triggerByDimEqual(triggerDim1_, prev_dim1_local_, triggerDim2_, prev_dim2_local_, workFlow_BufferKernel_, workFlow_use_BufferedKernel_));
                    triggered_in_process_ = true;
                    triggered_in_process_by_numOfKSpace_triggerDim1_ = true;
                }
            }

            prev_acq_header_ = *(m1->getObjectPtr());
        }
        else
        {
            // trigger whenever the Dim2 is changed
            if ( curr_dim2_!= prev_dim2_local_ )
            {
                GDEBUG_CONDITION_STREAM(verboseMode_, "Trigger Dim2 : " << gtPlusISMRMRDReconUtil<ValueType>().getISMRMRDDimName(triggerDim2_ ) << " = " << prev_dim2_local_);
                GADGET_CHECK_RETURN_FALSE(triggerByDimEqual(triggerDim2_, prev_dim2_local_, workFlow_BufferKernel_, workFlow_use_BufferedKernel_));
                triggered_in_process_ = true;
            }

            prev_acq_header_ = *(m1->getObjectPtr());
        }
    }

    return true;
}

bool GtPlusAccumulatorWorkOrderTriggerGadget::
triggerWorkOrderLastCountInClose(Gadgetron::ISMRMRDDIM& triggerDim1_, Gadgetron::ISMRMRDDIM& triggerDim2_, int numOfKSpace_triggerDim1_)
{
    GDEBUG_CONDITION_STREAM(verboseMode_, "Current Dim1 InClose : " << gtPlusISMRMRDReconUtil<ValueType>().getISMRMRDDimName(triggerDim1_ ) << " = " << curr_dim1_);
    GDEBUG_CONDITION_STREAM(verboseMode_, "Current Dim2 InClose : " << gtPlusISMRMRDReconUtil<ValueType>().getISMRMRDDimName(triggerDim2_ ) << " = " << curr_dim2_);

    if ( prev_dim1_ != curr_dim1_ )
    {
        count_dim1_++;
    }

    bool workFlow_BufferKernel_ = false;
    bool workFlow_use_BufferedKernel_ = false;

    int numOfAcquiredKSpaceForTriggerDim1 = numOfKSpace_triggerDim1_;
    if ( workOrder_.CalibMode_ == ISMRMRD_interleaved )
    {
        numOfAcquiredKSpaceForTriggerDim1 = (int)(numOfKSpace_triggerDim1_ * workOrder_.acceFactorE1_ * workOrder_.acceFactorE2_);
    }

    size_t prev_dim1_local_ = prev_dim1_;
    size_t prev_dim2_local_ = prev_dim2_;

    prev_dim1_ = curr_dim1_;
    prev_dim2_ = curr_dim2_;

    // trigger whenever the Dim2 is changed
    if (  triggerDim1_==DIM_NONE && triggerDim2_!=DIM_NONE  )
    {
        GDEBUG_CONDITION_STREAM(verboseMode_, "Trigger Dim2 : " << gtPlusISMRMRDReconUtil<ValueType>().getISMRMRDDimName(triggerDim2_ ) << " = " << prev_dim2_local_);
        GADGET_CHECK_RETURN_FALSE(triggerByDimEqual(triggerDim2_, prev_dim2_local_, workFlow_BufferKernel_, workFlow_use_BufferedKernel_));
    }

    // trigger whenever the Dim1 is changed
    if (  triggerDim1_!=DIM_NONE && triggerDim2_==DIM_NONE  )
    {
        if ( numOfKSpace_triggerDim1_ > 0 )
        {
            if ( count_dim1_ <= numOfAcquiredKSpaceForTriggerDim1 )
            {
                GDEBUG_CONDITION_STREAM(verboseMode_, "Trigger Dim1 : " << gtPlusISMRMRDReconUtil<ValueType>().getISMRMRDDimName(triggerDim1_ ) << " <= " << prev_dim1_local_);
                workFlow_BufferKernel_ = true;
                workFlow_use_BufferedKernel_ = false;
                GADGET_CHECK_RETURN_FALSE(triggerByDimLessEqual(triggerDim1_, prev_dim1_local_, workFlow_BufferKernel_, workFlow_use_BufferedKernel_));
            }
            else if ( count_dim1_ > numOfAcquiredKSpaceForTriggerDim1 )
            {
                GDEBUG_CONDITION_STREAM(verboseMode_, "Trigger Dim1 : " << gtPlusISMRMRDReconUtil<ValueType>().getISMRMRDDimName(triggerDim1_ ) << " = " << prev_dim1_local_);
                workFlow_BufferKernel_ = false;
                workFlow_use_BufferedKernel_ = true;
                GADGET_CHECK_RETURN_FALSE(triggerByDimEqual(triggerDim1_, prev_dim1_local_, workFlow_BufferKernel_, workFlow_use_BufferedKernel_));
            }
        }
        else
        {
            GDEBUG_CONDITION_STREAM(verboseMode_, "Trigger Dim1 : " << gtPlusISMRMRDReconUtil<ValueType>().getISMRMRDDimName(triggerDim1_ ) << " = " << prev_dim1_local_);
            GADGET_CHECK_RETURN_FALSE(triggerByDimEqual(triggerDim1_, prev_dim1_local_, workFlow_BufferKernel_, workFlow_use_BufferedKernel_));
        }
    }

    if (  triggerDim1_!=DIM_NONE && triggerDim2_!=DIM_NONE  )
    {
        if ( numOfKSpace_triggerDim1_ > 0 )
        {
            if ( count_dim1_ <= numOfAcquiredKSpaceForTriggerDim1 ) // no more data will be available, so have to do the recon
            {
                GDEBUG_CONDITION_STREAM(verboseMode_, "Trigger Dim1 : " << gtPlusISMRMRDReconUtil<ValueType>().getISMRMRDDimName(triggerDim1_ ) << " <= " << prev_dim1_local_);
                workFlow_BufferKernel_ = true;
                workFlow_use_BufferedKernel_ = false;
                GADGET_CHECK_RETURN_FALSE(triggerByDim1LessEqualDim2Equal(triggerDim1_, prev_dim1_local_, triggerDim2_, prev_dim2_local_, workFlow_BufferKernel_, workFlow_use_BufferedKernel_));
            }
            else if ( count_dim1_ > numOfAcquiredKSpaceForTriggerDim1 )
            {
                GDEBUG_CONDITION_STREAM(verboseMode_, "Trigger Dim1 : " << gtPlusISMRMRDReconUtil<ValueType>().getISMRMRDDimName(triggerDim1_ ) << " = " << prev_dim1_local_);
                workFlow_BufferKernel_ = false;
                workFlow_use_BufferedKernel_ = true;
                GADGET_CHECK_RETURN_FALSE(triggerByDimEqual(triggerDim1_, prev_dim1_local_, triggerDim2_, prev_dim2_local_, workFlow_BufferKernel_, workFlow_use_BufferedKernel_));
            }
        }
        else
        {
            GDEBUG_CONDITION_STREAM(verboseMode_, "Trigger Dim2 : " << gtPlusISMRMRDReconUtil<ValueType>().getISMRMRDDimName(triggerDim2_ ) << " = " << prev_dim2_local_);
            GADGET_CHECK_RETURN_FALSE(triggerByDimEqual(triggerDim2_, prev_dim2_local_, workFlow_BufferKernel_, workFlow_use_BufferedKernel_));
        }
    }

    return true;
}

bool GtPlusAccumulatorWorkOrderTriggerGadget::checkStatus(uint64_t flag, int samples, 
    bool& bIsKSpace, bool& bIsRef, bool& bIsNoise, bool& bIsPhaseCorr, bool& bIsReflect, bool& bIsOther,
    bool& bIsNavigator, bool& bIsRTFeedback, bool& bIsHPFeedback, bool& bIsDummyScan)
{
    bIsNoise = ISMRMRD::FlagBit(ISMRMRD::ISMRMRD_ACQ_IS_NOISE_MEASUREMENT).isSet(flag);
    bool is_ref = ISMRMRD::FlagBit(ISMRMRD::ISMRMRD_ACQ_IS_PARALLEL_CALIBRATION).isSet(flag);
    bool is_ref_kspace = ISMRMRD::FlagBit(ISMRMRD::ISMRMRD_ACQ_IS_PARALLEL_CALIBRATION_AND_IMAGING).isSet(flag);
    bIsReflect = ISMRMRD::FlagBit(ISMRMRD::ISMRMRD_ACQ_IS_REVERSE).isSet(flag);
    bIsPhaseCorr = ISMRMRD::FlagBit(ISMRMRD::ISMRMRD_ACQ_IS_PHASECORR_DATA).isSet(flag);
    bIsNavigator = ISMRMRD::FlagBit(ISMRMRD::ISMRMRD_ACQ_IS_NAVIGATION_DATA).isSet(flag);
    bIsRTFeedback = ISMRMRD::FlagBit(ISMRMRD::ISMRMRD_ACQ_IS_RTFEEDBACK_DATA).isSet(flag);
    bIsHPFeedback = ISMRMRD::FlagBit(ISMRMRD::ISMRMRD_ACQ_IS_HPFEEDBACK_DATA).isSet(flag);
    bIsDummyScan = ISMRMRD::FlagBit(ISMRMRD::ISMRMRD_ACQ_IS_DUMMYSCAN_DATA).isSet(flag);

    bIsKSpace = false;
    bIsRef = false;
    bIsOther = false;

    if ( bIsNoise || bIsDummyScan )
    {
        return true;
    }

    if ( workOrder_.CalibMode_==ISMRMRD_noacceleration )
    {
        bIsKSpace = true;
        bIsRef = false;
    }

    // in interleaved mode, only store the image data
    if ( workOrder_.CalibMode_==ISMRMRD_interleaved )
    {
        bIsKSpace = true;
        bIsRef = false;
    }

    // in embedded, kspace stores only the undersampled lines
    // ref stores all lines used for references
    if ( workOrder_.CalibMode_==ISMRMRD_embedded )
    {
        if ( is_ref && !is_ref_kspace )
        {
            bIsKSpace = false;
            bIsRef = true;
        }

        if ( !is_ref && is_ref_kspace )
        {
            bIsKSpace = true;
            bIsRef = true;
        }

        if ( is_ref && is_ref_kspace )
        {
            bIsKSpace = true;
            bIsRef = true;
        }

        if ( !is_ref && !is_ref_kspace )
        {
            bIsKSpace = true;
            bIsRef = false;
        }
    }

    // in separate mode
    if ( workOrder_.CalibMode_==ISMRMRD_separate 
    || workOrder_.CalibMode_==ISMRMRD_external )
    {
        if ( is_ref )
        {
            bIsKSpace = false;
            bIsRef = true;
        }

        if ( !is_ref )
        {
            bIsKSpace = true;
            bIsRef = false;
        }
    }

    // store other data, e.g. AIF
    // only for tpat
    if ( !is_ref && !is_ref_kspace && (samples != meas_max_ro_) )
    {
        bIsOther = true;
        bIsKSpace = false;
        bIsRef = false;
    }

    return true;
}

bool GtPlusAccumulatorWorkOrderTriggerGadget::storeImageData(GadgetContainerMessage<ISMRMRD::AcquisitionHeader>* m1, GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2, bool isReflect)
{
    try
    {
        size_t ii;
        size_t samples =  m1->getObjectPtr()->number_of_samples;
        ISMRMRD::EncodingCounters idx = m1->getObjectPtr()->idx;

        /*if ( workOrder_.retro_gated_images_ == 0 )
        {*/
            idx.segment = 0; // combine the segments
        //}

        if ( workOrder_.data_.get_number_of_elements() <= 0 )
        {
            meas_max_channel_ = m1->getObjectPtr()->active_channels;

            size_t E1 = workOrder_.kSpaceMaxEncode1_+1;
            size_t E2 = workOrder_.kSpaceMaxEncode2_+1;
            if ( E2 == 0 ) E2 = 1;

            if ( E1 < matrix_size_encoding_[1] ) E1 = matrix_size_encoding_[1];
            if ( E2 < matrix_size_encoding_[2] ) E2 = matrix_size_encoding_[2];

            if ( samples > meas_max_ro_ ) meas_max_ro_ = samples;

            // find the loop counter boundary and allocate the buffer
            GDEBUG_CONDITION_STREAM(verboseMode_, "[RO E1 Cha Slice E2 Con Phase Rep Set Seg Ave] = [" 
                               << meas_max_ro_ 
                               << " " << E1 
                               << " " << meas_max_channel_ 
                               << " " << meas_max_idx_.slice+1 
                               << " " << E2 
                               << " " << meas_max_idx_.contrast+1 
                               << " " << meas_max_idx_.phase+1 
                               << " " << meas_max_idx_.repetition+1 
                               << " " << meas_max_idx_.set+1 
                               << " " << meas_max_idx_.segment+1 
                               << " " << meas_max_idx_.average+1 << "]");

            dimensions_.clear();
            dimensions_.push_back(meas_max_ro_);
            dimensions_.push_back(E1);
            dimensions_.push_back(meas_max_channel_);
            dimensions_.push_back(meas_max_idx_.slice+1);
            dimensions_.push_back(E2);
            dimensions_.push_back(meas_max_idx_.contrast+1);
            dimensions_.push_back(meas_max_idx_.phase+1);
            dimensions_.push_back(meas_max_idx_.repetition+1);
            dimensions_.push_back(meas_max_idx_.set+1);
            dimensions_.push_back(meas_max_idx_.segment+1);
            dimensions_.push_back(meas_max_idx_.average+1);

            size_t N = dimensions_.size();
            for ( ii=0; ii<N; ii++ )
            {
                GDEBUG_CONDITION_STREAM(verboseMode_, "dimensions_[" << ii << "] = " << dimensions_[ii]);
            }

            // allocate data buffer
            try
            {
                workOrder_.data_.create(&dimensions_);
                Gadgetron::clear(workOrder_.data_);

                std::vector<size_t> reflect_dimensions_(dimensions_);
                reflect_dimensions_[0] = 1;
                reflect_dimensions_[2] = 1;
                workOrder_.reflect_.create(&reflect_dimensions_);
                Gadgetron::clear(workOrder_.reflect_);

                std::vector<size_t> dim(dimensions_);
                dim[0] = 1;
                dim[2] = 1;
                workOrder_.time_stamp_.create(dim);
                Gadgetron::fill(workOrder_.time_stamp_, (real_value_type)(-1) );

                workOrder_.physio_time_stamp_.create(dim);
                Gadgetron::fill(workOrder_.physio_time_stamp_, (real_value_type)(-1) );
            }
            catch(...)
            {
                GDEBUG("Failed create buffer\n");
                return false;
            }

            // allocate message buffer
            size_t matrix_size[GT_DIM_NUM];
            for ( ii=0; ii<GT_DIM_NUM; ii++ )
            {
                matrix_size[ii] = dimensions_[ii];
            }

            if (!(messageImage_ = new GtPlusGadgetImageArray(matrix_size))) 
            {
                GDEBUG("Failed create buffer\n");
                return false;
            }
        }

        // if necessary, shift the E1/E2 indexes
        idx.kspace_encode_step_1 += space_matrix_offset_E1_;
        idx.kspace_encode_step_2 += space_matrix_offset_E2_;

        if ( idx.kspace_encode_step_1 >= dimensions_[1] )
        {
            return true;
        }

        size_t dataN = workOrder_.data_.get_number_of_elements();
        std::complex<float>* b = workOrder_.data_.begin();
        std::complex<float>* d = m2->getObjectPtr()->get_data_ptr();
        if (samples != static_cast<size_t>(dimensions_[0])) 
        {
            GDEBUG("Wrong number of samples received\n");
            return false;
        }

        //Copy the data for all the channels
        hoNDArray<std::complex<float> > reflectBuf;
        if ( isReflect )
        {
            reflectBuf.create(samples);
        }

        std::vector<size_t> pos(GT_DIM_NUM);
        for (size_t c = 0; c < m1->getObjectPtr()->active_channels; c++) 
        {
            pos[0] = 0;
            pos[1] = idx.kspace_encode_step_1;
            pos[2] = c;
            pos[3] = idx.slice;
            pos[4] = idx.kspace_encode_step_2;
            pos[5] = idx.contrast;
            pos[6] = idx.phase;
            pos[7] = idx.repetition;
            pos[8] = idx.set;
            pos[9] = idx.segment;
            pos[10] = idx.average;
            size_t offsetBuffer = workOrder_.data_.calculate_offset(pos);

            if ( offsetBuffer >= dataN )
            {
                break;
            }

            if ( offsetBuffer >= dataN )
            {
                break;
            }

            if ( isReflect )
            {
                for ( size_t s=0; s<samples; s++ )
                {
                    reflectBuf(samples-1-s) = d[c*samples+s];
                }

                memcpy(b+offsetBuffer, reflectBuf.begin(), sizeof(std::complex<float>)*samples);
            }
            else
            {
                memcpy(b+offsetBuffer, d+c*samples, sizeof(std::complex<float>)*samples);
            }

            pos[2] = 0;
            offsetBuffer = workOrder_.reflect_.calculate_offset(pos);
            workOrder_.reflect_.at(offsetBuffer) = isReflect;
        }

        if ( !fillImageInfo(m1, messageImage_, idx) )
        {
            GDEBUG("Failed in fillImageInfo(m1, messageImage_, idx)\n");
            return false;
        }
    }
    catch(...)
    {
        GDEBUG("Errors in GtPlusAccumulatorWorkOrderTriggerGadget::storeImageData(...) ... \n");
        return false;
    }

    return true;
}

bool GtPlusAccumulatorWorkOrderTriggerGadget::
storeRefData(GadgetContainerMessage<ISMRMRD::AcquisitionHeader>* m1, GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2, bool isReflect)
{
    try
    {
        size_t ii;
        size_t samples =  m1->getObjectPtr()->number_of_samples;
        ISMRMRD::EncodingCounters idx = m1->getObjectPtr()->idx;

        /*if ( workOrder_.retro_gated_images_ == 0 )
        {*/
            idx.segment = 0; // combine the segments
        //}

        if ( workOrder_.ref_.get_number_of_elements() <= 0 )
        {
            meas_max_channel_ = m1->getObjectPtr()->active_channels;

            size_t E1 = workOrder_.kSpaceMaxEncode1_+1;
            size_t E2 = workOrder_.kSpaceMaxEncode2_+1;
            if ( E2 == 0 ) E2 = 1;

            if ( E1 < matrix_size_encoding_[1] ) E1 = matrix_size_encoding_[1];
            if ( E2 < matrix_size_encoding_[2] ) E2 = matrix_size_encoding_[2];

            size_t RO = meas_max_ro_;

            if ( (samples < meas_max_ro_) 
                && (( workOrder_.CalibMode_==ISMRMRD_separate || workOrder_.CalibMode_==ISMRMRD_external )) )
            {
                RO = samples;
            }

            if ( RO < samples ) RO = samples;

            // find the loop counter boundary and allocate the buffer
            GDEBUG_CONDITION_STREAM(verboseMode_, "[RO E1 Cha Slice E2 Con Phase Rep Set Seg Ave] = [" 
                               << RO 
                               << " " << E1 
                               << " " << meas_max_channel_ 
                               << " " << meas_max_idx_.slice+1 
                               << " " << E2 
                               << " " << meas_max_idx_.contrast+1 
                               << " " << meas_max_idx_.phase+1 
                               << " " << meas_max_idx_.repetition+1 
                               << " " << meas_max_idx_.set+1 
                               << " " << meas_max_idx_.segment+1 
                               << " " << meas_max_idx_.average+1 << "]");

            dimensions_.clear();
            dimensions_.push_back(RO);
            dimensions_.push_back(E1);
            dimensions_.push_back(meas_max_channel_);
            dimensions_.push_back(meas_max_idx_.slice+1);
            dimensions_.push_back(E2);
            dimensions_.push_back(meas_max_idx_.contrast+1);
            dimensions_.push_back(meas_max_idx_.phase+1);
            dimensions_.push_back(meas_max_idx_.repetition+1);
            dimensions_.push_back(meas_max_idx_.set+1);
            dimensions_.push_back(meas_max_idx_.segment+1);
            dimensions_.push_back(meas_max_idx_.average+1);

            size_t N = dimensions_.size();
            for ( ii=0; ii<N; ii++ )
            {
                GDEBUG_CONDITION_STREAM(verboseMode_, "ref dimensions_[" << ii << "] = " << dimensions_[ii]);
            }

            // allocate data buffer
            try
            {
                workOrder_.ref_.create(&dimensions_);
                Gadgetron::clear(workOrder_.ref_);

                std::vector<size_t> reflect_dimensions_(dimensions_);
                reflect_dimensions_[0] = 1;
                reflect_dimensions_[2] = 1;
                workOrder_.reflect_ref_.create(&reflect_dimensions_);
                Gadgetron::clear(workOrder_.reflect_ref_);
            }
            catch(...)
            {
                GDEBUG("Failed create ref buffer\n");
                return false;
            }
        }

        // if necessary, shift the E1/E2 indexes
        if ( workOrder_.CalibMode_ == ISMRMRD_embedded )
        {
            if ( workOrder_.start_E1_ > 0 )
            {
                idx.kspace_encode_step_1 += workOrder_.start_E1_;
            }

            if ( workOrder_.start_E2_ > 0 )
            {
                idx.kspace_encode_step_2 += workOrder_.start_E2_;
            }
        }

        // for the seperate or external mode, store the maximal idx
        if ( (workOrder_.CalibMode_ == ISMRMRD_separate) || (workOrder_.CalibMode_ == ISMRMRD_external) )
        {
            if ( idx.kspace_encode_step_1 > meas_max_idx_ref_.kspace_encode_step_1 )    meas_max_idx_ref_.kspace_encode_step_1 = idx.kspace_encode_step_1;
            if ( idx.kspace_encode_step_2 > meas_max_idx_ref_.kspace_encode_step_2 )    meas_max_idx_ref_.kspace_encode_step_2 = idx.kspace_encode_step_2;
            if ( idx.average > meas_max_idx_ref_.average )                              meas_max_idx_ref_.average = idx.average;
            if ( idx.slice > meas_max_idx_ref_.slice )                                  meas_max_idx_ref_.slice = idx.slice;
            if ( idx.contrast > meas_max_idx_ref_.contrast )                            meas_max_idx_ref_.contrast = idx.contrast;
            if ( idx.phase > meas_max_idx_ref_.phase )                                  meas_max_idx_ref_.phase = idx.phase;
            if ( idx.repetition > meas_max_idx_ref_.repetition )                        meas_max_idx_ref_.repetition = idx.repetition;
            if ( idx.set > meas_max_idx_ref_.set )                                      meas_max_idx_ref_.set = idx.set;
            if ( idx.segment > meas_max_idx_ref_.segment )                              meas_max_idx_ref_.segment = idx.segment;
            if ( idx.average > meas_max_idx_ref_.average )                              meas_max_idx_ref_.average = idx.average;

            size_t ii;
            for ( ii=0; ii<ISMRMRD::ISMRMRD_USER_INTS; ii++ )
            {
                if ( idx.user[ii] > meas_max_idx_ref_.user[ii] ) meas_max_idx_ref_.user[ii] = idx.user[ii];
            }
        }

        size_t refN = workOrder_.ref_.get_number_of_elements();
        std::complex<float>* b = workOrder_.ref_.begin();
        std::complex<float>* d = m2->getObjectPtr()->get_data_ptr();
        if (samples != static_cast<size_t>(workOrder_.ref_.get_size(0))) 
        {
            GDEBUG("Wrong number of samples received\n");
            return false;
        }

        //Copy the data for all the channels
        hoNDArray<std::complex<float> > reflectBuf;
        if ( isReflect )
        {
            reflectBuf.create(samples);
        }

        std::vector<size_t> pos(GT_DIM_NUM);
        for (uint16_t c = 0; c < m1->getObjectPtr()->active_channels; c++) 
        {
            pos[0] = 0;
            pos[1] = idx.kspace_encode_step_1;
            pos[2] = c;
            pos[3] = idx.slice;
            pos[4] = idx.kspace_encode_step_2;
            pos[5] = idx.contrast;
            pos[6] = idx.phase;
            pos[7] = idx.repetition;
            pos[8] = idx.set;
            pos[9] = idx.segment;
            pos[10] = idx.average;

            size_t offsetBuffer = workOrder_.ref_.calculate_offset(pos);
            if ( offsetBuffer >= refN )
            {
                break;
            }

            if ( isReflect )
            {
                for ( size_t s=0; s<samples; s++ )
                {
                    reflectBuf(samples-1-s) = d[c*samples+s];
                }

                memcpy(b+offsetBuffer, reflectBuf.begin(), sizeof(std::complex<float>)*samples);
            }
            else
            {
                memcpy(b+offsetBuffer, d+c*samples, sizeof(std::complex<float>)*samples);
            }

            pos[2] = 0;
            offsetBuffer = workOrder_.reflect_ref_.calculate_offset(pos);
            workOrder_.reflect_ref_.at(offsetBuffer) = isReflect;
        }

        // if it is embedded mode, store the acquisition and physio time stamp
        if ( workOrder_.CalibMode_ == ISMRMRD_embedded )
        {
            ind_time_stamp_[0] = 0;
            ind_time_stamp_[1] = idx.kspace_encode_step_1;
            ind_time_stamp_[2] = 0;
            ind_time_stamp_[3] = idx.slice;
            ind_time_stamp_[4] = idx.kspace_encode_step_2;
            ind_time_stamp_[5] = idx.contrast;
            ind_time_stamp_[6] = idx.phase;
            ind_time_stamp_[7] = idx.repetition;
            ind_time_stamp_[8] = idx.set;
            ind_time_stamp_[9] = idx.segment;
            ind_time_stamp_[10] = idx.average;

            workOrder_.time_stamp_(ind_time_stamp_) = (real_value_type)(m1->getObjectPtr()->acquisition_time_stamp) * timeStampResolution_;
            workOrder_.physio_time_stamp_(ind_time_stamp_) = (real_value_type)(m1->getObjectPtr()->physiology_time_stamp[0]) * timeStampResolution_;
        }
    }
    catch(...)
    {
        GDEBUG("Errors in GtPlusAccumulatorWorkOrderTriggerGadget::storeRefData(...) ... \n");
        return false;
    }

    return true;
}

bool GtPlusAccumulatorWorkOrderTriggerGadget::
fillBuffer(ReadOutBufferType& readOutBuffer, BufferType& buf, ReflectBufferType& reflectBuf)
{
    try
    {
        // find the maximal dimension of all buffered ICE readout
        size_t numOfReadOuts = readOutBuffer.size();
        ISMRMRD::EncodingCounters max_idx;
        max_idx.kspace_encode_step_1 = 0;
        max_idx.average = 0;
        max_idx.slice = 0;
        max_idx.kspace_encode_step_2 = 0;
        max_idx.contrast = 0;
        max_idx.phase = 0;
        max_idx.repetition = 0;
        max_idx.set = 0;
        max_idx.segment = 0;
        max_idx.average = 0;
        size_t max_channel = 0;
        size_t max_col = 0;

        size_t a;
        for (a = 0; a < numOfReadOuts; a++) 
        {
            ISMRMRD::EncodingCounters idx = readOutBuffer[a].acqHead_.idx;

            if ( readOutBuffer[a].acqHead_.number_of_samples > max_col ) 
                max_col=readOutBuffer[a].acqHead_.number_of_samples;

            if ( idx.kspace_encode_step_1 > max_idx.kspace_encode_step_1 ) 
                max_idx.kspace_encode_step_1=idx.kspace_encode_step_1;

            if ( idx.slice > max_idx.slice ) 
                max_idx.slice = idx.slice;

            if ( idx.kspace_encode_step_2 > max_idx.kspace_encode_step_2 ) 
                max_idx.kspace_encode_step_2 = idx.kspace_encode_step_2;

            if ( idx.contrast > max_idx.contrast ) 
                max_idx.contrast = idx.contrast;

            if ( idx.phase > max_idx.phase ) 
                max_idx.phase = idx.phase;

            if ( idx.repetition > max_idx.repetition ) 
                max_idx.repetition = idx.repetition;

            if ( idx.set > max_idx.set ) 
                max_idx.set = idx.set;

            if ( idx.segment > max_idx.segment ) 
                max_idx.segment = idx.segment;

            if ( idx.average > max_idx.average ) 
                max_idx.average = idx.average;

            if ( readOutBuffer[a].acqHead_.active_channels > max_channel ) 
                max_channel = readOutBuffer[a].acqHead_.active_channels;
        }

        GDEBUG_CONDITION_STREAM(verboseMode_, "[RO E1 Cha Slice E2 Contrast Phase Rep Set Seg Ave] = [" 
                               << max_col 
                               << " " << max_idx.kspace_encode_step_1+1 
                               << " " << max_channel 
                               << " " << max_idx.slice+1 
                               << " " << max_idx.kspace_encode_step_2+1 
                               << " " << max_idx.contrast+1 
                               << " " << max_idx.phase+1 
                               << " " << max_idx.repetition+1 
                               << " " << max_idx.set+1 
                               << " " << max_idx.segment+1 
                               << " " << max_idx.average+1 << "]");

        // alloate buffer for data
        std::vector<size_t> dims(GT_DIM_NUM);
        dims[0] = max_col;
        dims[1] = max_idx.kspace_encode_step_1+1;
        dims[2] = max_channel;
        dims[3] = max_idx.slice+1;
        dims[4] = max_idx.kspace_encode_step_2+1;
        dims[5] = max_idx.contrast+1;
        dims[6] = max_idx.phase+1;
        dims[7] = max_idx.repetition+1;
        dims[8] = max_idx.set+1;
        dims[9] = max_idx.segment+1;
        dims[10] = max_idx.average+1;

        try
        {
            buf.create(&dims);
            Gadgetron::clear(buf);

            std::vector<size_t> reflect_dims(dims);
            reflect_dims[0] = 1;
            reflect_dims[2] = 1;
            reflectBuf.create(&reflect_dims);
            Gadgetron::clear(reflectBuf);
        }
        catch(...)
        {
            GDEBUG("Failed create buffer\n");
            return false;
        }

        std::complex<float>* b = buf.begin();

        // copy the data
        uint16_t c;
        std::vector<size_t> pos(GT_DIM_NUM);

        for ( a=0; a<numOfReadOuts; a++) 
        {
            ISMRMRD::EncodingCounters idx = readOutBuffer[a].acqHead_.idx;
            std::complex<float>* d = const_cast<std::complex<float>*>(readOutBuffer[a].data_.begin());

            for ( c=0; c<readOutBuffer[a].acqHead_.active_channels; c++) 
            {
                pos[0] = 0;
                pos[1] = idx.kspace_encode_step_1;
                pos[2] = c;
                pos[3] = idx.slice;
                pos[4] = idx.kspace_encode_step_2;
                pos[5] = idx.contrast;
                pos[6] = idx.phase;
                pos[7] = idx.repetition;
                pos[8] = idx.set;
                pos[9] = idx.segment;
                pos[10] = idx.average;
                long long offsetBuffer = buf.calculate_offset(pos);

                memcpy(b+offsetBuffer, d+c*readOutBuffer[a].acqHead_.number_of_samples, sizeof(std::complex<float>)*readOutBuffer[a].acqHead_.number_of_samples);

                pos[2] = 0;
                offsetBuffer = reflectBuf.calculate_offset(pos);
                reflectBuf.at(offsetBuffer) = readOutBuffer[a].isReflect_;
            }
        }
    }
    catch(...)
    {
        GDEBUG("Errors in GtPlusAccumulatorWorkOrderTriggerGadget::fillBuffer(...) ... \n");
        return false;
    }

    return true;
}

//XUE-TODO: Functions DO NOT return booleans in the Gadgetron
bool GtPlusAccumulatorWorkOrderTriggerGadget::fillImageInfo(GadgetContainerMessage<ISMRMRD::AcquisitionHeader>* m1, GtPlusGadgetImageArray* messageImage, const ISMRMRD::EncodingCounters& idx)
{

    try
    {
        // fill the message info
        size_t offset = messageImage->get_offset(idx.slice, idx.kspace_encode_step_2, idx.contrast, idx.phase, idx.repetition, idx.set, idx.segment, idx.average);

        if( (offset >= messageImage->max_num_of_images_)
            || (idx.slice>=messageImage->matrix_size[3])
            || (idx.kspace_encode_step_2>=messageImage->matrix_size[4])
            || (idx.contrast>=messageImage->matrix_size[5])
            || (idx.phase>=messageImage->matrix_size[6])
            || (idx.repetition>=messageImage->matrix_size[7])
            || (idx.set>=messageImage->matrix_size[8])
            || (idx.segment>=messageImage->matrix_size[9])
            || (idx.average>=messageImage->matrix_size[10]) )
        {
            GWARN_STREAM("Incoming image is over the boundary of buffer [SLC E2 CON PHS REP SET SEG AVE] = [ " 
                                                                            << idx.slice << " " << idx.kspace_encode_step_2 << " " 
                                                                            << idx.contrast << " " << idx.phase << " " 
                                                                            << idx.repetition << " " << idx.set << " " 
                                                                            << idx.segment << " " << idx.average << " ] ");
            return true;
        }

        if( offset >= messageImage->max_num_of_images_ )
        {
            GWARN_STREAM("Incoming image is over the boundary of buffer [SLC E2 CON PHS REP SET SEG AVE] = [ " 
                                                                            << idx.slice << " " << idx.kspace_encode_step_2 << " " 
                                                                            << idx.contrast << " " << idx.phase << " " << idx.repetition << " " 
                                                                            << idx.set << " " << idx.segment << " " << idx.average << " ] ");
            return true;
        }

        // if it is the first acq in a slice, fill in all information
        bool is_first_acq_in_slice = ISMRMRD::FlagBit(ISMRMRD::ISMRMRD_ACQ_FIRST_IN_SLICE).isSet(m1->getObjectPtr()->flags);

        /*if ( is_first_acq_in_slice 
            || ( messageImage->imageArray_[offset].version==0 
                    && messageImage->imageArray_[offset].flags==0 
                    && messageImage->imageArray_[offset].measurement_uid==0 ) )*/
        if ( messageImage->imageArray_[offset].version==0 
                    && messageImage->imageArray_[offset].flags==0 
                    && messageImage->imageArray_[offset].measurement_uid==0 )
        {
            GDEBUG_CONDITION_STREAM(verboseMode_, "--> buffer image header - offset = " << offset << " - [SLC E2 CON PHS REP SET SEG AVE] = [" 
                                                                      << idx.slice << " " 
                                                                      << idx.kspace_encode_step_2 << " " 
                                                                      << idx.contrast << " " 
                                                                      << idx.phase << " " 
                                                                      << idx.repetition << " " 
                                                                      << idx.set << " " 
                                                                      << idx.segment << " " 
                                                                      << idx.average << "]");

            messageImage->imageArray_[offset].version = m1->getObjectPtr()->version;
            messageImage->imageArray_[offset].flags = m1->getObjectPtr()->flags;
            messageImage->imageArray_[offset].measurement_uid = m1->getObjectPtr()->measurement_uid;

            //messageImage->imageArray_[offset].matrix_size[0] = dimensions_[0];
            //messageImage->imageArray_[offset].matrix_size[1] = dimensions_[1];
            //messageImage->imageArray_[offset].matrix_size[2] = dimensions_[2];

            messageImage->imageArray_[offset].set_matrix_size(0, dimensions_[0]);
            messageImage->imageArray_[offset].set_matrix_size(1, dimensions_[1]);
            messageImage->imageArray_[offset].set_matrix_size(2, dimensions_[2]);

            messageImage->imageArray_[offset].field_of_view[0] = field_of_view_recon_[0];
            messageImage->imageArray_[offset].field_of_view[1] = field_of_view_recon_[1];
            messageImage->imageArray_[offset].field_of_view[2] = field_of_view_recon_[2];

            messageImage->imageArray_[offset].channels = m1->getObjectPtr()->active_channels;

            messageImage->imageArray_[offset].position[0] = m1->getObjectPtr()->position[0];
            messageImage->imageArray_[offset].position[1] = m1->getObjectPtr()->position[1];
            messageImage->imageArray_[offset].position[2] = m1->getObjectPtr()->position[2];

            //messageImage->imageArray_[offset].quaternion[0] = m1->getObjectPtr()->quaternion[0];
            //messageImage->imageArray_[offset].quaternion[1] = m1->getObjectPtr()->quaternion[1];
            //messageImage->imageArray_[offset].quaternion[2] = m1->getObjectPtr()->quaternion[2];
            //messageImage->imageArray_[offset].quaternion[3] = m1->getObjectPtr()->quaternion[3];

            messageImage->imageArray_[offset].read_dir[0] = m1->getObjectPtr()->read_dir[0];
            messageImage->imageArray_[offset].read_dir[1] = m1->getObjectPtr()->read_dir[1];
            messageImage->imageArray_[offset].read_dir[2] = m1->getObjectPtr()->read_dir[2];

            messageImage->imageArray_[offset].phase_dir[0] = m1->getObjectPtr()->phase_dir[0];
            messageImage->imageArray_[offset].phase_dir[1] = m1->getObjectPtr()->phase_dir[1];
            messageImage->imageArray_[offset].phase_dir[2] = m1->getObjectPtr()->phase_dir[2];

            messageImage->imageArray_[offset].slice_dir[0] = m1->getObjectPtr()->slice_dir[0];
            messageImage->imageArray_[offset].slice_dir[1] = m1->getObjectPtr()->slice_dir[1];
            messageImage->imageArray_[offset].slice_dir[2] = m1->getObjectPtr()->slice_dir[2];

            messageImage->imageArray_[offset].patient_table_position[0] = m1->getObjectPtr()->patient_table_position[0];
            messageImage->imageArray_[offset].patient_table_position[1] = m1->getObjectPtr()->patient_table_position[1];
            messageImage->imageArray_[offset].patient_table_position[2] = m1->getObjectPtr()->patient_table_position[2];

            messageImage->imageArray_[offset].average = m1->getObjectPtr()->idx.average;
            messageImage->imageArray_[offset].slice = m1->getObjectPtr()->idx.slice;
            messageImage->imageArray_[offset].contrast = m1->getObjectPtr()->idx.contrast;
            messageImage->imageArray_[offset].phase = m1->getObjectPtr()->idx.phase;
            messageImage->imageArray_[offset].repetition = m1->getObjectPtr()->idx.repetition;
            messageImage->imageArray_[offset].set = m1->getObjectPtr()->idx.set;
            messageImage->imageArray_[offset].average = m1->getObjectPtr()->idx.average;

            messageImage->imageArray_[offset].acquisition_time_stamp = m1->getObjectPtr()->acquisition_time_stamp;

            messageImage->imageArray_[offset].physiology_time_stamp[0] = m1->getObjectPtr()->physiology_time_stamp[0];
            messageImage->imageArray_[offset].physiology_time_stamp[1] = m1->getObjectPtr()->physiology_time_stamp[1];
            messageImage->imageArray_[offset].physiology_time_stamp[2] = m1->getObjectPtr()->physiology_time_stamp[2];

            messageImage->imageArray_[offset].data_type = ISMRMRD::ISMRMRD_CXFLOAT;

            messageImage->imageArray_[offset].image_type = ISMRMRD::ISMRMRD_IMTYPE_MAGNITUDE;

            messageImage->imageArray_[offset].image_index = (uint16_t)(++image_counter_);
            messageImage->imageArray_[offset].image_series_index = (uint16_t)image_series_;

            // need to store the free user parameters
            memcpy(messageImage->imageArray_[offset].user_int, m1->getObjectPtr()->user_int, sizeof(int32_t)*8);
            memcpy(messageImage->imageArray_[offset].user_float, m1->getObjectPtr()->user_float, sizeof(float)*8);
        }

        // whether or not this acq is the first in a slice, we need to fill the TimeStamps and PMUTimeStamps
        if ( idx.kspace_encode_step_1 < messageImage->imageArray_[offset].time_stamps.size() )
        {
            messageImage->imageArray_[offset].time_stamps[idx.kspace_encode_step_1] = m1->getObjectPtr()->acquisition_time_stamp;
            messageImage->imageArray_[offset].pmu_time_stamps[idx.kspace_encode_step_1] = m1->getObjectPtr()->physiology_time_stamp[0];

            ind_time_stamp_[0] = 0;
            ind_time_stamp_[1] = idx.kspace_encode_step_1;
            ind_time_stamp_[2] = 0;
            ind_time_stamp_[3] = idx.slice;
            ind_time_stamp_[4] = idx.kspace_encode_step_2;
            ind_time_stamp_[5] = idx.contrast;
            ind_time_stamp_[6] = idx.phase;
            ind_time_stamp_[7] = idx.repetition;
            ind_time_stamp_[8] = idx.set;
            ind_time_stamp_[9] = idx.segment;
            ind_time_stamp_[10] = idx.average;

            workOrder_.time_stamp_(ind_time_stamp_) = (real_value_type)(m1->getObjectPtr()->acquisition_time_stamp) * timeStampResolution_;
            workOrder_.physio_time_stamp_(ind_time_stamp_) = (real_value_type)(m1->getObjectPtr()->physiology_time_stamp[0]) * timeStampResolution_;
        }
    }
    catch(...)
    {
        GDEBUG("Errors in GtPlusAccumulatorWorkOrderTriggerGadget::fillImageInfo(...) ... \n");
        return false;
    }

    return true;
}

size_t GtPlusAccumulatorWorkOrderTriggerGadget::
computeEncodedSizeE1(size_t centerE1, size_t maxE1)
{
    size_t E1;
    if ( (maxE1+1)%2 == 0 )
    {
        E1 = 2*centerE1;
    }
    else
    {
        E1 = 2*centerE1+1;
    }

    return E1;
}

size_t GtPlusAccumulatorWorkOrderTriggerGadget::
computeEncodedSizeE2(size_t centerE2, size_t maxE2)
{
    size_t E2;
    if ( (maxE2+1)%2 == 0 )
    {
        E2 = 2*centerE2;
    }
    else
    {
        E2 = 2*centerE2+1;
    }

    return E2;
}

bool GtPlusAccumulatorWorkOrderTriggerGadget::
triggerByDimEqual(Gadgetron::ISMRMRDDIM& triggerDim, size_t value, bool workFlow_BufferKernel_, bool workFlow_use_BufferedKernel_)
{
    try
    {
        GDEBUG_CONDITION_STREAM(verboseMode_, "GtPlusAccumulatorWorkOrderTriggerGadget - triggerByDimEqual(triggerDim, value) ... ");

        GadgetContainerMessage<GtPlusGadgetImageArray>* cm1 = new GadgetContainerMessage<GtPlusGadgetImageArray>();
        GadgetContainerMessage< WorkOrderType >* cm2 = new GadgetContainerMessage< WorkOrderType >();
        cm1->cont(cm2);

        workOrder_.duplicate(*cm2->getObjectPtr());
        cm2->getObjectPtr()->workFlow_BufferKernel_ = workFlow_BufferKernel_;
        cm2->getObjectPtr()->workFlow_use_BufferedKernel_ = workFlow_use_BufferedKernel_;

        bool lessEqual = false;

        // copy the image content
        GADGET_CHECK_RETURN_FALSE(gtPlusISMRMRDReconUtil<ValueType>().extractSubArrayForDim(workOrder_.data_, cm2->getObjectPtr()->data_, triggerDim, value, lessEqual));
        GADGET_CHECK_RETURN_FALSE(gtPlusISMRMRDReconUtil<unsigned short>().extractSubArrayForDim(workOrder_.reflect_, cm2->getObjectPtr()->reflect_, triggerDim, value, lessEqual));
        GADGET_CHECK_RETURN_FALSE(gtPlusISMRMRDReconUtil<real_value_type>().extractSubArrayForDim(workOrder_.time_stamp_, cm2->getObjectPtr()->time_stamp_, triggerDim, value, lessEqual));
        GADGET_CHECK_RETURN_FALSE(gtPlusISMRMRDReconUtil<real_value_type>().extractSubArrayForDim(workOrder_.physio_time_stamp_, cm2->getObjectPtr()->physio_time_stamp_, triggerDim, value, lessEqual));

        // copy the ref
        if ( workOrder_.ref_.get_number_of_elements()>0 
                && workOrder_.ref_.get_number_of_dimensions()==workOrder_.data_.get_number_of_dimensions() )
        {
            GADGET_CHECK_RETURN_FALSE(gtPlusISMRMRDReconUtil<ValueType>().extractSubArrayForDim(workOrder_.ref_, cm2->getObjectPtr()->ref_, triggerDim, value, lessEqual));
            GADGET_CHECK_RETURN_FALSE(gtPlusISMRMRDReconUtil<unsigned short>().extractSubArrayForDim(workOrder_.reflect_ref_, cm2->getObjectPtr()->reflect_ref_, triggerDim, value, lessEqual));

            // for seperate and external mode, further truncate the reference data
            if ( (workOrder_.CalibMode_ == ISMRMRD_separate) || (workOrder_.CalibMode_ == ISMRMRD_external) )
            {
                hoNDArray<ValueType> ref;
                hoNDArray<unsigned short> reflect_ref;

                GADGET_CHECK_RETURN_FALSE(gtPlusISMRMRDReconUtil<ValueType>().extractSubArrayForMaxEncodingCounters(cm2->getObjectPtr()->ref_, ref, meas_max_idx_ref_));
                GADGET_CHECK_RETURN_FALSE(gtPlusISMRMRDReconUtil<unsigned short>().extractSubArrayForMaxEncodingCounters(cm2->getObjectPtr()->reflect_ref_, reflect_ref, meas_max_idx_ref_));

                cm2->getObjectPtr()->ref_ = ref;
                cm2->getObjectPtr()->reflect_ref_ = reflect_ref;
            }
        }

        // copy the message image array
        GADGET_CHECK_RETURN_FALSE(messageImage_->extractGadgetImageArrayEqual(triggerDim, value, *(cm1->getObjectPtr()) ));

        if (!phaseCorrBuffer_.empty())
        {
            GDEBUG_CONDITION_STREAM(verboseMode_, "GtPlusAccumulatorWorkOrderTriggerGadget - phase correction signal found : " << phaseCorrBuffer_.size());

            if ( !fillBuffer(phaseCorrBuffer_, workOrder_.phaseCorr_, workOrder_.reflect_phaseCorr_) )
            {
                GDEBUG("fillBuffer(phaseCorrBuffer_) failed ... \n");
                cm1->release();
                return false;
            }

            cm2->getObjectPtr()->phaseCorr_ = workOrder_.phaseCorr_;
            cm2->getObjectPtr()->reflect_phaseCorr_ = workOrder_.reflect_phaseCorr_;
        }

        if (!noiseBuffer_.empty())
        {
            GDEBUG_CONDITION_STREAM(verboseMode_, "GtPlusAccumulatorGadget - noise signal found : " << noiseBuffer_.size());

            ReflectBufferType tmpBuf;
            if ( !fillBuffer(noiseBuffer_, workOrder_.noise_, tmpBuf) )
            {
                GDEBUG("fillBuffer(noiseBuffer_) failed ... \n");
                cm1->release();
                return false;
            }

            cm2->getObjectPtr()->noise_ = workOrder_.noise_;
        }

        if (!otherBuffer_.empty())
        {
            GDEBUG_CONDITION_STREAM(verboseMode_, "GtPlusAccumulatorGadget - other signal found : " << otherBuffer_.size());

            if ( !fillBuffer(otherBuffer_, workOrder_.other_, workOrder_.reflect_other_) )
            {
                GDEBUG("fillBuffer(otherBuffer_) failed ... \n");
                cm1->release();
                return false;
            }

            GADGET_CHECK_RETURN_FALSE(gtPlusISMRMRDReconUtil<ValueType>().extractSubArrayForDim(workOrder_.other_, cm2->getObjectPtr()->other_, triggerDim, value, lessEqual));

            cm2->getObjectPtr()->reflect_other_ = workOrder_.reflect_other_;
        }

        // send to next gadget
        if (this->next()->putq(cm1) < 0) 
        {
            return false;
        }
    }
    catch(...)
    {
        GERROR_STREAM("Errors in GtPlusAccumulatorWorkOrderTriggerGadget::triggerByDimEqual(triggerDim, value) ... ");
        return false;
    }

    return true;
}

bool GtPlusAccumulatorWorkOrderTriggerGadget::
triggerByDimLessEqual(Gadgetron::ISMRMRDDIM& triggerDim, size_t value, bool workFlow_BufferKernel_, bool workFlow_use_BufferedKernel_)
{
    try
    {
        GDEBUG_CONDITION_STREAM(verboseMode_, "GtPlusAccumulatorWorkOrderTriggerGadget - triggerByDimEqual(triggerDim, value) ... ");

        GadgetContainerMessage<GtPlusGadgetImageArray>* cm1 = new GadgetContainerMessage<GtPlusGadgetImageArray>();
        GadgetContainerMessage< WorkOrderType >* cm2 = new GadgetContainerMessage< WorkOrderType >();
        cm1->cont(cm2);

        workOrder_.duplicate(*cm2->getObjectPtr());
        cm2->getObjectPtr()->workFlow_BufferKernel_ = workFlow_BufferKernel_;
        cm2->getObjectPtr()->workFlow_use_BufferedKernel_ = workFlow_use_BufferedKernel_;

        bool lessEqual = true;

        // copy the image content
        GADGET_CHECK_RETURN_FALSE(gtPlusISMRMRDReconUtil<ValueType>().extractSubArrayForDim(workOrder_.data_, cm2->getObjectPtr()->data_, triggerDim, value, lessEqual));
        GADGET_CHECK_RETURN_FALSE(gtPlusISMRMRDReconUtil<unsigned short>().extractSubArrayForDim(workOrder_.reflect_, cm2->getObjectPtr()->reflect_, triggerDim, value, lessEqual));
        GADGET_CHECK_RETURN_FALSE(gtPlusISMRMRDReconUtil<real_value_type>().extractSubArrayForDim(workOrder_.time_stamp_, cm2->getObjectPtr()->time_stamp_, triggerDim, value, lessEqual));
        GADGET_CHECK_RETURN_FALSE(gtPlusISMRMRDReconUtil<real_value_type>().extractSubArrayForDim(workOrder_.physio_time_stamp_, cm2->getObjectPtr()->physio_time_stamp_, triggerDim, value, lessEqual));

        // copy the ref
        if ( workOrder_.ref_.get_number_of_elements()>0 
                && workOrder_.ref_.get_number_of_dimensions()==workOrder_.data_.get_number_of_dimensions() )
        {
            GADGET_CHECK_RETURN_FALSE(gtPlusISMRMRDReconUtil<ValueType>().extractSubArrayForDim(workOrder_.ref_, cm2->getObjectPtr()->ref_, triggerDim, value, lessEqual));
            GADGET_CHECK_RETURN_FALSE(gtPlusISMRMRDReconUtil<unsigned short>().extractSubArrayForDim(workOrder_.reflect_ref_, cm2->getObjectPtr()->reflect_ref_, triggerDim, value, lessEqual));

            // for seperate and external mode, further truncate the reference data
            if ( (workOrder_.CalibMode_ == ISMRMRD_separate) || (workOrder_.CalibMode_ == ISMRMRD_external) )
            {
                hoNDArray<ValueType> ref;
                hoNDArray<unsigned short> reflect_ref;

                GADGET_CHECK_RETURN_FALSE(gtPlusISMRMRDReconUtil<ValueType>().extractSubArrayForMaxEncodingCounters(cm2->getObjectPtr()->ref_, ref, meas_max_idx_ref_));
                GADGET_CHECK_RETURN_FALSE(gtPlusISMRMRDReconUtil<unsigned short>().extractSubArrayForMaxEncodingCounters(cm2->getObjectPtr()->reflect_ref_, reflect_ref, meas_max_idx_ref_));

                cm2->getObjectPtr()->ref_ = ref;
                cm2->getObjectPtr()->reflect_ref_ = reflect_ref;
            }
        }

        // copy the message image array
        GADGET_CHECK_RETURN_FALSE(messageImage_->extractGadgetImageArrayLessEqual(triggerDim, value, *(cm1->getObjectPtr()) ));

        if (!phaseCorrBuffer_.empty())
        {
            GDEBUG_CONDITION_STREAM(verboseMode_, "GtPlusAccumulatorWorkOrderTriggerGadget - phase correction signal found : " << phaseCorrBuffer_.size());

            if ( !fillBuffer(phaseCorrBuffer_, workOrder_.phaseCorr_, workOrder_.reflect_phaseCorr_) )
            {
                GDEBUG("fillBuffer(phaseCorrBuffer_) failed ... \n");
                cm1->release();
                return false;
            }

            cm2->getObjectPtr()->phaseCorr_ = workOrder_.phaseCorr_;
            cm2->getObjectPtr()->reflect_phaseCorr_ = workOrder_.reflect_phaseCorr_;
        }

        if (!noiseBuffer_.empty())
        {
            GDEBUG_CONDITION_STREAM(verboseMode_, "GtPlusAccumulatorGadget - noise signal found : " << noiseBuffer_.size());

            ReflectBufferType tmpBuf;
            if ( !fillBuffer(noiseBuffer_, workOrder_.noise_, tmpBuf) )
            {
                GDEBUG("fillBuffer(noiseBuffer_) failed ... \n");
                cm1->release();
                return false;
            }

            cm2->getObjectPtr()->noise_ = workOrder_.noise_;
        }

        if (!otherBuffer_.empty())
        {
            GDEBUG_CONDITION_STREAM(verboseMode_, "GtPlusAccumulatorGadget - other signal found : " << otherBuffer_.size());

            if ( !fillBuffer(otherBuffer_, workOrder_.other_, workOrder_.reflect_other_) )
            {
                GDEBUG("fillBuffer(otherBuffer_) failed ... \n");
                cm1->release();
                return false;
            }

            GADGET_CHECK_RETURN_FALSE(gtPlusISMRMRDReconUtil<ValueType>().extractSubArrayForDim(workOrder_.other_, cm2->getObjectPtr()->other_, triggerDim, value, lessEqual));

            cm2->getObjectPtr()->reflect_other_ = workOrder_.reflect_other_;
        }

        // send to next gadget
        if (this->next()->putq(cm1) < 0) 
        {
            return false;
        }
    }
    catch(...)
    {
        GERROR_STREAM("Errors in GtPlusAccumulatorWorkOrderTriggerGadget::triggerByDimLessEqual(triggerDim, value) ... ");
        return false;
    }

    return true;
}

bool GtPlusAccumulatorWorkOrderTriggerGadget::
triggerByDimEqual(Gadgetron::ISMRMRDDIM& triggerDim1, size_t value1, Gadgetron::ISMRMRDDIM& triggerDim2, size_t value2, bool workFlow_BufferKernel_, bool workFlow_use_BufferedKernel_)
{
    try
    {
        GDEBUG_CONDITION_STREAM(verboseMode_, "GtPlusAccumulatorWorkOrderTriggerGadget - triggerByDimEqual(triggerDim1, value1, triggerDim2, value2) ... ");

        GadgetContainerMessage<GtPlusGadgetImageArray>* cm1 = new GadgetContainerMessage<GtPlusGadgetImageArray>();
        GadgetContainerMessage< WorkOrderType >* cm2 = new GadgetContainerMessage< WorkOrderType >();
        cm1->cont(cm2);

        workOrder_.duplicate(*cm2->getObjectPtr());
        cm2->getObjectPtr()->workFlow_BufferKernel_ = workFlow_BufferKernel_;
        cm2->getObjectPtr()->workFlow_use_BufferedKernel_ = workFlow_use_BufferedKernel_;

        bool lessEqual = false;

        // copy the image content
        GADGET_CHECK_RETURN_FALSE(gtPlusISMRMRDReconUtil<ValueType>().extractSubArrayForDim(workOrder_.data_, cm2->getObjectPtr()->data_, triggerDim1, value1, triggerDim2, value2, lessEqual));
        GADGET_CHECK_RETURN_FALSE(gtPlusISMRMRDReconUtil<unsigned short>().extractSubArrayForDim(workOrder_.reflect_, cm2->getObjectPtr()->reflect_, triggerDim1, value1, triggerDim2, value2, lessEqual));
        GADGET_CHECK_RETURN_FALSE(gtPlusISMRMRDReconUtil<real_value_type>().extractSubArrayForDim(workOrder_.time_stamp_, cm2->getObjectPtr()->time_stamp_, triggerDim1, value1, triggerDim2, value2, lessEqual));
        GADGET_CHECK_RETURN_FALSE(gtPlusISMRMRDReconUtil<real_value_type>().extractSubArrayForDim(workOrder_.physio_time_stamp_, cm2->getObjectPtr()->physio_time_stamp_, triggerDim1, value1, triggerDim2, value2, lessEqual));

        // copy the ref
        if ( workOrder_.ref_.get_number_of_elements()>0 
                && workOrder_.ref_.get_number_of_dimensions()==workOrder_.data_.get_number_of_dimensions() )
        {
            GADGET_CHECK_RETURN_FALSE(gtPlusISMRMRDReconUtil<ValueType>().extractSubArrayForDim(workOrder_.ref_, cm2->getObjectPtr()->ref_, triggerDim1, value1, triggerDim2, value2, lessEqual));
            GADGET_CHECK_RETURN_FALSE(gtPlusISMRMRDReconUtil<unsigned short>().extractSubArrayForDim(workOrder_.reflect_ref_, cm2->getObjectPtr()->reflect_ref_, triggerDim1, value1, triggerDim2, value2, lessEqual));

            // for seperate and external mode, further truncate the reference data
            if ( (workOrder_.CalibMode_ == ISMRMRD_separate) || (workOrder_.CalibMode_ == ISMRMRD_external) )
            {
                hoNDArray<ValueType> ref;
                hoNDArray<unsigned short> reflect_ref;

                GADGET_CHECK_RETURN_FALSE(gtPlusISMRMRDReconUtil<ValueType>().extractSubArrayForMaxEncodingCounters(cm2->getObjectPtr()->ref_, ref, meas_max_idx_ref_));
                GADGET_CHECK_RETURN_FALSE(gtPlusISMRMRDReconUtil<unsigned short>().extractSubArrayForMaxEncodingCounters(cm2->getObjectPtr()->reflect_ref_, reflect_ref, meas_max_idx_ref_));

                cm2->getObjectPtr()->ref_ = ref;
                cm2->getObjectPtr()->reflect_ref_ = reflect_ref;
            }
        }

        // copy the message image array
        GADGET_CHECK_RETURN_FALSE(messageImage_->extractGadgetImageArrayEqual(triggerDim1, value1, triggerDim2, value2, *(cm1->getObjectPtr()) ));

        if (!phaseCorrBuffer_.empty())
        {
            GDEBUG_CONDITION_STREAM(verboseMode_, "GtPlusAccumulatorWorkOrderTriggerGadget - phase correction signal found : " << phaseCorrBuffer_.size());

            if ( !fillBuffer(phaseCorrBuffer_, workOrder_.phaseCorr_, workOrder_.reflect_phaseCorr_) )
            {
                GDEBUG("fillBuffer(phaseCorrBuffer_) failed ... \n");
                cm1->release();
                return false;
            }

            cm2->getObjectPtr()->phaseCorr_ = workOrder_.phaseCorr_;
            cm2->getObjectPtr()->reflect_phaseCorr_ = workOrder_.reflect_phaseCorr_;
        }

        if (!noiseBuffer_.empty())
        {
            GDEBUG_CONDITION_STREAM(verboseMode_, "GtPlusAccumulatorGadget - noise signal found : " << noiseBuffer_.size());

            ReflectBufferType tmpBuf;
            if ( !fillBuffer(noiseBuffer_, workOrder_.noise_, tmpBuf) )
            {
                GDEBUG("fillBuffer(noiseBuffer_) failed ... \n");
                cm1->release();
                return false;
            }

            cm2->getObjectPtr()->noise_ = workOrder_.noise_;
        }

        if (!otherBuffer_.empty())
        {
            GDEBUG_CONDITION_STREAM(verboseMode_, "GtPlusAccumulatorGadget - other signal found : " << otherBuffer_.size());

            if ( !fillBuffer(otherBuffer_, workOrder_.other_, workOrder_.reflect_other_) )
            {
                GDEBUG("fillBuffer(otherBuffer_) failed ... \n");
                cm1->release();
                return false;
            }

            GADGET_CHECK_RETURN_FALSE(gtPlusISMRMRDReconUtil<ValueType>().extractSubArrayForDim(workOrder_.other_, cm2->getObjectPtr()->other_, triggerDim1, value1, false));

            cm2->getObjectPtr()->reflect_other_ = workOrder_.reflect_other_;
        }

        // send to next gadget
        if (this->next()->putq(cm1) < 0) 
        {
            return false;
        }
    }
    catch(...)
    {
        GERROR_STREAM("Errors in GtPlusAccumulatorWorkOrderTriggerGadget::triggerByDimEqual(triggerDim1, value1, triggerDim2, value2) ... ");
        return false;
    }

    return true;
}


bool GtPlusAccumulatorWorkOrderTriggerGadget::
triggerByDim1LessEqualDim2Equal(Gadgetron::ISMRMRDDIM& triggerDim1, size_t value1, Gadgetron::ISMRMRDDIM& triggerDim2, size_t value2, bool workFlow_BufferKernel_, bool workFlow_use_BufferedKernel_)
{
    try
    {
        GDEBUG_CONDITION_STREAM(verboseMode_, "GtPlusAccumulatorWorkOrderTriggerGadget - triggerByDim1LessEqualDim2Equal(triggerDim1, value1, triggerDim2, value2) ... ");

        GadgetContainerMessage<GtPlusGadgetImageArray>* cm1 = new GadgetContainerMessage<GtPlusGadgetImageArray>();
        GadgetContainerMessage< WorkOrderType >* cm2 = new GadgetContainerMessage< WorkOrderType >();

        workOrder_.duplicate(*cm2->getObjectPtr());
        cm2->getObjectPtr()->workFlow_BufferKernel_ = workFlow_BufferKernel_;
        cm2->getObjectPtr()->workFlow_use_BufferedKernel_ = workFlow_use_BufferedKernel_;

        cm1->cont(cm2);

        // copy the image content
        GADGET_CHECK_RETURN_FALSE(gtPlusISMRMRDReconUtil<ValueType>().extractSubArrayForDim1LessEqualDim2Equal(workOrder_.data_, cm2->getObjectPtr()->data_, triggerDim1, value1, triggerDim2, value2));
        GADGET_CHECK_RETURN_FALSE(gtPlusISMRMRDReconUtil<unsigned short>().extractSubArrayForDim1LessEqualDim2Equal(workOrder_.reflect_, cm2->getObjectPtr()->reflect_, triggerDim1, value1, triggerDim2, value2));
        GADGET_CHECK_RETURN_FALSE(gtPlusISMRMRDReconUtil<real_value_type>().extractSubArrayForDim1LessEqualDim2Equal(workOrder_.time_stamp_, cm2->getObjectPtr()->time_stamp_, triggerDim1, value1, triggerDim2, value2));
        GADGET_CHECK_RETURN_FALSE(gtPlusISMRMRDReconUtil<real_value_type>().extractSubArrayForDim1LessEqualDim2Equal(workOrder_.physio_time_stamp_, cm2->getObjectPtr()->physio_time_stamp_, triggerDim1, value1, triggerDim2, value2));

        // copy the ref
        if ( workOrder_.ref_.get_number_of_elements()>0 
                && workOrder_.ref_.get_number_of_dimensions()==workOrder_.data_.get_number_of_dimensions() )
        {
            GADGET_CHECK_RETURN_FALSE(gtPlusISMRMRDReconUtil<ValueType>().extractSubArrayForDim1LessEqualDim2Equal(workOrder_.ref_, cm2->getObjectPtr()->ref_, triggerDim1, value1, triggerDim2, value2));
            GADGET_CHECK_RETURN_FALSE(gtPlusISMRMRDReconUtil<unsigned short>().extractSubArrayForDim1LessEqualDim2Equal(workOrder_.reflect_ref_, cm2->getObjectPtr()->reflect_ref_, triggerDim1, value1, triggerDim2, value2));

            // for seperate and external mode, further truncate the reference data
            if ( (workOrder_.CalibMode_ == ISMRMRD_separate) || (workOrder_.CalibMode_ == ISMRMRD_external) )
            {
                hoNDArray<ValueType> ref;
                hoNDArray<unsigned short> reflect_ref;

                GADGET_CHECK_RETURN_FALSE(gtPlusISMRMRDReconUtil<ValueType>().extractSubArrayForMaxEncodingCounters(cm2->getObjectPtr()->ref_, ref, meas_max_idx_ref_));
                GADGET_CHECK_RETURN_FALSE(gtPlusISMRMRDReconUtil<unsigned short>().extractSubArrayForMaxEncodingCounters(cm2->getObjectPtr()->reflect_ref_, reflect_ref, meas_max_idx_ref_));

                cm2->getObjectPtr()->ref_ = ref;
                cm2->getObjectPtr()->reflect_ref_ = reflect_ref;
            }
        }

        // copy the message image array
        GADGET_CHECK_RETURN_FALSE(messageImage_->extractGadgetImageArray_Dim1LessEqual_Dim2Equal(triggerDim1, value1, triggerDim2, value2, *(cm1->getObjectPtr()) ));

        if (!phaseCorrBuffer_.empty())
        {
            GDEBUG_CONDITION_STREAM(verboseMode_, "GtPlusAccumulatorWorkOrderTriggerGadget - phase correction signal found : " << phaseCorrBuffer_.size());

            if ( !fillBuffer(phaseCorrBuffer_, workOrder_.phaseCorr_, workOrder_.reflect_phaseCorr_) )
            {
                GDEBUG("fillBuffer(phaseCorrBuffer_) failed ... \n");
                cm1->release();
                return false;
            }

            cm2->getObjectPtr()->phaseCorr_ = workOrder_.phaseCorr_;
            cm2->getObjectPtr()->reflect_phaseCorr_ = workOrder_.reflect_phaseCorr_;
        }

        if (!noiseBuffer_.empty())
        {
            GDEBUG_CONDITION_STREAM(verboseMode_, "GtPlusAccumulatorGadget - noise signal found : " << noiseBuffer_.size());

            ReflectBufferType tmpBuf;
            if ( !fillBuffer(noiseBuffer_, workOrder_.noise_, tmpBuf) )
            {
                GDEBUG("fillBuffer(noiseBuffer_) failed ... \n");
                cm1->release();
                return false;
            }

            cm2->getObjectPtr()->noise_ = workOrder_.noise_;
        }

        if (!otherBuffer_.empty())
        {
            GDEBUG_CONDITION_STREAM(verboseMode_, "GtPlusAccumulatorGadget - other signal found : " << otherBuffer_.size());

            if ( !fillBuffer(otherBuffer_, workOrder_.other_, workOrder_.reflect_other_) )
            {
                GDEBUG("fillBuffer(otherBuffer_) failed ... \n");
                cm1->release();
                return false;
            }

            GADGET_CHECK_RETURN_FALSE(gtPlusISMRMRDReconUtil<ValueType>().extractSubArrayForDim(workOrder_.other_, cm2->getObjectPtr()->other_, triggerDim1, value1, true));

            cm2->getObjectPtr()->reflect_other_ = workOrder_.reflect_other_;
        }

        // send to next gadget
        if (this->next()->putq(cm1) < 0) 
        {
            return false;
        }
    }
    catch(...)
    {
        GERROR_STREAM("Errors in GtPlusAccumulatorWorkOrderTriggerGadget::triggerByDim1LessEqualDim2Equal(triggerDim1, value1, triggerDim2, value2) ... ");
        return false;
    }

    return true;
}

bool GtPlusAccumulatorWorkOrderTriggerGadget::triggerWorkOrderAllInClose()
{
    try
    {
        GDEBUG_CONDITION_STREAM(verboseMode_, "GtPlusAccumulatorWorkOrderTriggerGadget - triggerWorkOrderAllInClose ... ");

        GadgetContainerMessage<GtPlusGadgetImageArray>* cm1 = new GadgetContainerMessage<GtPlusGadgetImageArray>();
        GadgetContainerMessage< WorkOrderType >* cm2 = new GadgetContainerMessage< WorkOrderType >();

        workOrder_.duplicate(*cm2->getObjectPtr());
        cm2->getObjectPtr()->workFlow_BufferKernel_ = false;
        cm2->getObjectPtr()->workFlow_use_BufferedKernel_ = false;

        cm1->cont(cm2);

        // copy the image content
        cm2->getObjectPtr()->data_ = workOrder_.data_;
        cm2->getObjectPtr()->time_stamp_ = workOrder_.time_stamp_;
        cm2->getObjectPtr()->physio_time_stamp_ = workOrder_.physio_time_stamp_;
        cm2->getObjectPtr()->reflect_ = workOrder_.reflect_;

        // copy the ref
        cm2->getObjectPtr()->ref_ = workOrder_.ref_;
        cm2->getObjectPtr()->reflect_ref_ = workOrder_.reflect_ref_;

        // for seperate and external mode, further truncate the reference data
        if ( (workOrder_.CalibMode_ == ISMRMRD_separate) || (workOrder_.CalibMode_ == ISMRMRD_external) )
        {
            hoNDArray<ValueType> ref;
            hoNDArray<unsigned short> reflect_ref;

            GADGET_CHECK_RETURN_FALSE(gtPlusISMRMRDReconUtil<ValueType>().extractSubArrayForMaxEncodingCounters(cm2->getObjectPtr()->ref_, ref, meas_max_idx_ref_));
            GADGET_CHECK_RETURN_FALSE(gtPlusISMRMRDReconUtil<unsigned short>().extractSubArrayForMaxEncodingCounters(cm2->getObjectPtr()->reflect_ref_, reflect_ref, meas_max_idx_ref_));

            cm2->getObjectPtr()->ref_ = ref;
            cm2->getObjectPtr()->reflect_ref_ = reflect_ref;
        }

        // copy the message image array
        GADGET_CHECK_RETURN_FALSE(cm1->getObjectPtr()->copy(*messageImage_));

        if (!phaseCorrBuffer_.empty())
        {
            GDEBUG_CONDITION_STREAM(verboseMode_, "GtPlusAccumulatorWorkOrderTriggerGadget - phase correction signal found : " << phaseCorrBuffer_.size());

            if ( !fillBuffer(phaseCorrBuffer_, workOrder_.phaseCorr_, workOrder_.reflect_phaseCorr_) )
            {
                GDEBUG("fillBuffer(phaseCorrBuffer_) failed ... \n");
                cm1->release();
                return false;
            }

            cm2->getObjectPtr()->phaseCorr_ = workOrder_.phaseCorr_;
            cm2->getObjectPtr()->reflect_phaseCorr_ = workOrder_.reflect_phaseCorr_;
        }

        if (!noiseBuffer_.empty())
        {
            GDEBUG_CONDITION_STREAM(verboseMode_, "GtPlusAccumulatorGadget - noise signal found : " << noiseBuffer_.size());

            ReflectBufferType tmpBuf;
            if ( !fillBuffer(noiseBuffer_, workOrder_.noise_, tmpBuf) )
            {
                GDEBUG("fillBuffer(noiseBuffer_) failed ... \n");
                cm1->release();
                return false;
            }

            cm2->getObjectPtr()->noise_ = workOrder_.noise_;
        }

        if (!otherBuffer_.empty())
        {
            GDEBUG_CONDITION_STREAM(verboseMode_, "GtPlusAccumulatorGadget - other signal found : " << otherBuffer_.size());

            if ( !fillBuffer(otherBuffer_, workOrder_.other_, workOrder_.reflect_other_) )
            {
                GDEBUG("fillBuffer(otherBuffer_) failed ... \n");
                cm1->release();
                return false;
            }

            cm2->getObjectPtr()->other_ = workOrder_.other_;
            cm2->getObjectPtr()->reflect_other_ = workOrder_.reflect_other_;
        }

        // send to next gadget
        if (this->next()->putq(cm1) < 0) 
        {
            return false;
        }
    }
    catch(...)
    {
        GERROR_STREAM("Errors in GtPlusAccumulatorWorkOrderTriggerGadget::triggerWorkOrderAllInClose() ... ");
        return false;
    }

    return true;
}

size_t GtPlusAccumulatorWorkOrderTriggerGadget::
getDimValue(const ISMRMRD::AcquisitionHeader& acqHeader, Gadgetron::ISMRMRDDIM& dim)
{
    if ( dim == DIM_Encoding1 )             return acqHeader.idx.kspace_encode_step_1;
    if ( dim == DIM_Slice )                 return acqHeader.idx.slice;
    if ( dim == DIM_Encoding2 )             return acqHeader.idx.kspace_encode_step_2;
    if ( dim == DIM_Contrast )              return acqHeader.idx.contrast;
    if ( dim == DIM_Phase )                 return acqHeader.idx.phase;
    if ( dim == DIM_Repetition )            return acqHeader.idx.repetition;
    if ( dim == DIM_Set )                   return acqHeader.idx.set;
    if ( dim == DIM_Segment )               return acqHeader.idx.segment;
    if ( dim == DIM_Average )               return acqHeader.idx.average;

    return 0;
}

void GtPlusAccumulatorWorkOrderTriggerGadget::
setDimValue(ISMRMRD::AcquisitionHeader& acqHeader, Gadgetron::ISMRMRDDIM& dim, size_t value)
{
    if ( dim == DIM_Encoding1 ) acqHeader.idx.kspace_encode_step_1  = (uint16_t)value;
    if ( dim == DIM_Slice ) acqHeader.idx.slice                     = (uint16_t)value;
    if ( dim == DIM_Encoding2 ) acqHeader.idx.kspace_encode_step_2  = (uint16_t)value;
    if ( dim == DIM_Contrast ) acqHeader.idx.contrast               = (uint16_t)value;
    if ( dim == DIM_Phase ) acqHeader.idx.phase                     = (uint16_t)value;
    if ( dim == DIM_Repetition ) acqHeader.idx.repetition           = (uint16_t)value;
    if ( dim == DIM_Set ) acqHeader.idx.set                         = (uint16_t)value;
    if ( dim == DIM_Segment ) acqHeader.idx.segment                 = (uint16_t)value;
    if ( dim == DIM_Average ) acqHeader.idx.average                 = (uint16_t)value;

    return;
}

int GtPlusAccumulatorWorkOrderTriggerGadget::close(unsigned long flags)
{
    GDEBUG_CONDITION_STREAM(true, "GtPlusAccumulatorWorkOrderTriggerGadget - close(flags) : " << flags);

    if ( BaseClass::close(flags) != GADGET_OK ) return GADGET_FAIL;

    if ( flags!=0 && !triggered_in_close_ )
    // if ( !triggered_in_close_ )
    {
        triggered_in_close_ = true;

        GDEBUG_CONDITION_STREAM(true, "GtPlusAccumulatorWorkOrderTriggerGadget - trigger in close(flags) ... ");

        if ( needTriggerWorkOrderAllInClose() )
        {
            // never been triggered, so need to trigger with all data buffered
            if ( !triggerWorkOrderAllInClose() )
            {
                GDEBUG("triggerWorkOrderAllInClose() failed ... \n");
                return GADGET_FAIL;
            }
        }
        else
        {
            // need to trigger the last portion of kspace
            //if ( !triggerWorkOrder(NULL, true, true) )
            //{
            //    GDEBUG("Failed triggerWorkOrder(inClose)\n");
            //    return GADGET_FAIL;
            //}
        }
    }

    // return BaseClass::close(flags);
    return GADGET_OK;
}

GADGET_FACTORY_DECLARE(GtPlusAccumulatorWorkOrderTriggerGadget)

}
