#include "GtPlusAccumulatorWorkOrderTriggerGadget.h"

using namespace Gadgetron::gtPlus;

namespace Gadgetron
{

GtPlusAccumulatorWorkOrderTriggerGadget::GtPlusAccumulatorWorkOrderTriggerGadget()
: image_counter_(0), image_series_(100), first_kspace_scan_(true), triggered_in_close_(false), triggered_in_process_(false), triggered_in_process_last_acq_(false), 
    prev_dim1_(-1), curr_dim1_(-1), prev_dim2_(-1), curr_dim2_(-1), count_dim1_(0), verboseMode_(false), other_kspace_matching_Dim_(DIM_NONE)
{
    space_matrix_offset_E1_ = 0;
    space_matrix_offset_E2_ = 0;

    gtPlusISMRMRDReconUtil<ValueType>().clearAcquisitionHeaderISMRMRD(prev_acq_header_);
    memset(&meas_max_idx_ref_, 0, sizeof(ISMRMRD::EncodingCounters));
}

GtPlusAccumulatorWorkOrderTriggerGadget::~GtPlusAccumulatorWorkOrderTriggerGadget()
{

}

// extract necessary configuration information from the xml
int GtPlusAccumulatorWorkOrderTriggerGadget::process_config(ACE_Message_Block* mb)
{
    // gadget parameters
    image_series_ = this->get_int_value("image_series");

    noacceleration_triggerDim1_ = gtPlus_util_.getISMRMRDDimFromName(*(this->get_string_value("noacceleration_triggerDim1")));
    noacceleration_triggerDim2_ = gtPlus_util_.getISMRMRDDimFromName(*(this->get_string_value("noacceleration_triggerDim2")));
    noacceleration_numOfKSpace_triggerDim1_ = this->get_int_value("noacceleration_numOfKSpace_triggerDim1"); 

    interleaved_triggerDim1_ = gtPlus_util_.getISMRMRDDimFromName(*(this->get_string_value("interleaved_triggerDim1")));
    interleaved_triggerDim2_ = gtPlus_util_.getISMRMRDDimFromName(*(this->get_string_value("interleaved_triggerDim2")));
    interleaved_numOfKSpace_triggerDim1_ = this->get_int_value("interleaved_numOfKSpace_triggerDim1"); 

    embedded_triggerDim1_ = gtPlus_util_.getISMRMRDDimFromName(*(this->get_string_value("embedded_triggerDim1")));
    embedded_triggerDim2_ = gtPlus_util_.getISMRMRDDimFromName(*(this->get_string_value("embedded_triggerDim2")));
    embedded_numOfKSpace_triggerDim1_ = this->get_int_value("embedded_numOfKSpace_triggerDim1");

    separate_triggerDim1_ = gtPlus_util_.getISMRMRDDimFromName(*(this->get_string_value("separate_triggerDim1")));
    separate_triggerDim2_ = gtPlus_util_.getISMRMRDDimFromName(*(this->get_string_value("separate_triggerDim2")));
    separate_numOfKSpace_triggerDim1_ = this->get_int_value("separate_numOfKSpace_triggerDim1");

    other_kspace_matching_Dim_ = gtPlus_util_.getISMRMRDDimFromName(*(this->get_string_value("other_kspace_matching_Dim")));

    verboseMode_ = this->get_bool_value("verboseMode");

    // ---------------------------------------------------------------------------------------------------------
    // pass the xml file
    boost::shared_ptr<ISMRMRD::ismrmrdHeader> cfg = parseIsmrmrdXMLHeader(std::string(mb->rd_ptr()));

    // seq object
    ISMRMRD::ismrmrdHeader::encoding_sequence e_seq = cfg->encoding();
    // This only supports two encoding spaces where the recon_space is the same size
    // e.g. Parallel imaging reference scan collected with GRE and data with EPI
    if (e_seq.size() > 2)
    {
      GADGET_DEBUG2("Number of encoding spaces: %d\n", e_seq.size());
      GADGET_DEBUG1("This simple GtPlusAccumulatorWorkOrderTriggerGadget only supports one encoding space\n");
      return GADGET_FAIL;
    } 
    else if (e_seq.size() == 2)
    {
      if (! ((e_seq[0].reconSpace().matrixSize().x() == e_seq[1].reconSpace().matrixSize().x()) &
             (e_seq[0].reconSpace().matrixSize().y() == e_seq[1].reconSpace().matrixSize().y()) &
             (e_seq[0].reconSpace().matrixSize().z() == e_seq[1].reconSpace().matrixSize().z()) &
             (e_seq[0].reconSpace().fieldOfView_mm().x() == e_seq[1].reconSpace().fieldOfView_mm().x()) &
             (e_seq[0].reconSpace().fieldOfView_mm().y() == e_seq[1].reconSpace().fieldOfView_mm().y()) &
             (e_seq[0].reconSpace().fieldOfView_mm().z() == e_seq[1].reconSpace().fieldOfView_mm().z())) )
      {
	GADGET_DEBUG2("Number of encoding spaces: %d\n", e_seq.size());
	GADGET_DEBUG1("This simple GtPlusAccumulatorWorkOrderTriggerGadget only supports two encoding spaces with identical recon spaces.\n");
	return GADGET_FAIL;
      }
    }

    // find out the PAT mode
    ISMRMRD::ismrmrdHeader::parallelImaging_optional p_imaging_type = cfg->parallelImaging();
    ISMRMRD::parallelImagingType p_imaging = *p_imaging_type;

    workOrder_.acceFactorE1_ = (size_t)(p_imaging.accelerationFactor().kspace_encoding_step_1());
    workOrder_.acceFactorE2_ = (size_t)(p_imaging.accelerationFactor().kspace_encoding_step_2());
    GADGET_CONDITION_MSG(verboseMode_, "acceFactorE1_ is " << workOrder_.acceFactorE1_);
    GADGET_CONDITION_MSG(verboseMode_, "acceFactorE2_ is " << workOrder_.acceFactorE2_);

    ISMRMRD::calibrationModeType calib = *(p_imaging.calibrationMode());
    if ( calib == ISMRMRD::calibrationModeType::interleaved )
    {
        workOrder_.CalibMode_ = Gadgetron::gtPlus::ISMRMRD_interleaved;
        GADGET_CONDITION_MSG(verboseMode_, "Calibration mode is interleaved");

        if ( p_imaging.interleavingDimension().present() )
        {
            if ( *(p_imaging.interleavingDimension()) == ISMRMRD::interleavingDimensionType::phase )
            {
                workOrder_.InterleaveDim_ = Gadgetron::gtPlus::DIM_Phase;
            }

            if ( *(p_imaging.interleavingDimension()) == ISMRMRD::interleavingDimensionType::repetition )
            {
                workOrder_.InterleaveDim_ = Gadgetron::gtPlus::DIM_Repetition;
            }

            if ( *(p_imaging.interleavingDimension()) == ISMRMRD::interleavingDimensionType::average )
            {
                workOrder_.InterleaveDim_ = Gadgetron::gtPlus::DIM_Average;
            }

            if ( *(p_imaging.interleavingDimension()) == ISMRMRD::interleavingDimensionType::contrast )
            {
                workOrder_.InterleaveDim_ = Gadgetron::gtPlus::DIM_Contrast;
            }

            if ( *(p_imaging.interleavingDimension()) == ISMRMRD::interleavingDimensionType::other )
            {
                workOrder_.InterleaveDim_ = Gadgetron::gtPlus::DIM_other1;
            }

            GADGET_CONDITION_MSG(verboseMode_, "InterleaveDim is " << gtPlus_util_.getISMRMRDDimName(workOrder_.InterleaveDim_));
        }
    }

    if ( calib == ISMRMRD::calibrationModeType::embedded )
    {
        workOrder_.CalibMode_ = Gadgetron::gtPlus::ISMRMRD_embedded;
        GADGET_CONDITION_MSG(verboseMode_, "Calibration mode is embedded");
    }

    if ( calib == ISMRMRD::calibrationModeType::separate )
    {
        workOrder_.CalibMode_ = Gadgetron::gtPlus::ISMRMRD_separate;
        GADGET_CONDITION_MSG(verboseMode_, "Calibration mode is separate");
    }

    if ( calib == ISMRMRD::calibrationModeType::external )
    {
        workOrder_.CalibMode_ = Gadgetron::gtPlus::ISMRMRD_external;
    }

    if ( calib == ISMRMRD::calibrationModeType::other && workOrder_.acceFactorE1_==1 && workOrder_.acceFactorE2_==1 )
    {
        workOrder_.CalibMode_ = Gadgetron::gtPlus::ISMRMRD_noacceleration;
        // workOrder_.CalibMode_ = Gadgetron::gtPlus::ISMRMRD_interleaved;
        workOrder_.acceFactorE1_=1;
        // workOrder_.InterleaveDim_ = Gadgetron::gtPlus::DIM_Phase;
    }

    if ( calib == ISMRMRD::calibrationModeType::other && (workOrder_.acceFactorE1_>1 || workOrder_.acceFactorE2_>1) )
    {
        //workOrder_.CalibMode_ = Gadgetron::gtPlus::ISMRMRD_other;
        workOrder_.CalibMode_ = Gadgetron::gtPlus::ISMRMRD_interleaved;
        workOrder_.acceFactorE1_=2;
        workOrder_.InterleaveDim_ = Gadgetron::gtPlus::DIM_Phase;
    }

    // ---------------------------------------------------------------------------------------------------------

    // find out the encoding space 
    ISMRMRD::encodingSpaceType e_space = (*e_seq.begin()).encodedSpace();
    ISMRMRD::encodingSpaceType r_space = (*e_seq.begin()).reconSpace();
    ISMRMRD::encodingLimitsType e_limits = (*e_seq.begin()).encodingLimits();

    matrix_size_encoding_[0] = e_space.matrixSize().x();
    matrix_size_encoding_[1] = e_space.matrixSize().y();
    matrix_size_encoding_[2] = e_space.matrixSize().z();
    GADGET_CONDITION_MSG(verboseMode_, "Encoding matrix size: " << matrix_size_encoding_[0] << " " << matrix_size_encoding_[1] << " " << matrix_size_encoding_[2]);

    field_of_view_encoding_[0] = e_space.fieldOfView_mm().x();
    field_of_view_encoding_[1] = e_space.fieldOfView_mm().y();
    field_of_view_encoding_[2] = e_space.fieldOfView_mm().z();
    GADGET_CONDITION_MSG(verboseMode_, "Encoding field_of_view : " << field_of_view_encoding_[0] << " " << field_of_view_encoding_[1] << " " << field_of_view_encoding_[2]);

    // find the recon space
    matrix_size_recon_[0] = r_space.matrixSize().x();
    matrix_size_recon_[1] = r_space.matrixSize().y();
    matrix_size_recon_[2] = r_space.matrixSize().z();
    GADGET_CONDITION_MSG(verboseMode_, "Recon matrix size : " << matrix_size_recon_[0] << " " << matrix_size_recon_[1] << " " << matrix_size_recon_[2]);

    field_of_view_recon_[0] = r_space.fieldOfView_mm().x();
    field_of_view_recon_[1] = r_space.fieldOfView_mm().y();
    field_of_view_recon_[2] = r_space.fieldOfView_mm().z();
    GADGET_CONDITION_MSG(verboseMode_, "Recon field_of_view :  " << field_of_view_recon_[0] << " " << field_of_view_recon_[1] << " " << field_of_view_recon_[2]);

    // ---------------------------------------------------------------------------------------------------------
    // handle partial fourier
    //workOrder_.kSpaceCenterEncode1_ = e_limits.kspace_encoding_step_1().get().center();
    //GADGET_CONDITION_MSG(verboseMode_, "kSpaceCenterEncode1_ is " << workOrder_.kSpaceCenterEncode1_);

    //workOrder_.kSpaceCenterEncode2_ = e_limits.kspace_encoding_step_2().get().center();
    //GADGET_CONDITION_MSG(verboseMode_, "kSpaceCenterEncode2_ is " << workOrder_.kSpaceCenterEncode2_);

    workOrder_.kSpaceMaxEncode1_ = matrix_size_encoding_[1]-1; // e_limits.kspace_encoding_step_1().get().maximum();
    GADGET_CONDITION_MSG(verboseMode_, "matrix size kSpaceMaxEncode1_ is " << workOrder_.kSpaceMaxEncode1_);

    workOrder_.kSpaceMaxEncode2_ = matrix_size_encoding_[2]-1; // e_limits.kspace_encoding_step_2().get().maximum();
    GADGET_CONDITION_MSG(verboseMode_, "matrix size kSpaceMaxEncode2_ is " << workOrder_.kSpaceMaxEncode2_);

    space_size_[1] = workOrder_.kSpaceMaxEncode1_+1;
    space_size_[2] = workOrder_.kSpaceMaxEncode2_+1;

    max_sampled_E1_ = e_limits.kspace_encoding_step_1().get().maximum();
    max_sampled_E2_ = e_limits.kspace_encoding_step_2().get().maximum();

    GADGET_CONDITION_MSG(verboseMode_, "max_sampled_E1_ is " << max_sampled_E1_);
    GADGET_CONDITION_MSG(verboseMode_, "max_sampled_E2_ is " << max_sampled_E2_);

    center_line_E1_ = e_limits.kspace_encoding_step_1().get().center();
    center_line_E2_ = e_limits.kspace_encoding_step_2().get().center();

    GADGET_CONDITION_MSG(verboseMode_, "center_line_E1_ is " << center_line_E1_);
    GADGET_CONDITION_MSG(verboseMode_, "center_line_E2_ is " << center_line_E2_);

    workOrder_.kSpaceCenterEncode1_ = center_line_E1_;
    GADGET_CONDITION_MSG(verboseMode_, "kSpaceCenterEncode1_ is " << workOrder_.kSpaceCenterEncode1_);

    workOrder_.kSpaceCenterEncode2_ = center_line_E2_;
    GADGET_CONDITION_MSG(verboseMode_, "kSpaceCenterEncode2_ is " << workOrder_.kSpaceCenterEncode2_);

    // if partial fourier or asymmetric echo is used, correct the kSpaceCenter
    //if ( space_size_[1]-matrix_size_encoding_[1] > workOrder_.acceFactorE1_ )
    //{
    //    GADGET_CONDITION_MSG(verboseMode_, "Partial fourier along E1 ... ");
    //    //if ( GT_ABS(matrix_size_encoding_[1]/workOrder_.acceFactorE1_ - std::floor(matrix_size_encoding_[1]/workOrder_.acceFactorE1_)) > FLT_EPSILON )
    //    //{
    //    //    GADGET_WARN_MSG("matrix_size_[1] is not multiplied by acceFactorE1_ ... ");
    //    //    matrix_size_encoding_[1] = (std::floor(matrix_size_encoding_[1]/workOrder_.acceFactorE1_)+1)*workOrder_.acceFactorE1_;
    //    //}

    //    if ( 2*workOrder_.kSpaceCenterEncode1_ > (matrix_size_encoding_[1]+1) )
    //    {
    //        space_matrix_offset_E1_ = 0;

    //        workOrder_.start_E2_ = 0;
    //        workOrder_.end_E2_ = matrix_size_encoding_[1];
    //    }
    //    else
    //    {
    //        space_matrix_offset_E1_ = space_size_[1] - matrix_size_encoding_[1];

    //        workOrder_.start_E1_ = space_matrix_offset_E1_;
    //        workOrder_.end_E1_ = workOrder_.kSpaceMaxEncode1_;
    //    }
    //}
    //else
    //{
    //    space_matrix_offset_E1_ = 0;
    //}

    //if ( space_size_[2]-matrix_size_encoding_[2] > workOrder_.acceFactorE2_ )
    //{
    //    GADGET_CONDITION_MSG(verboseMode_, "Partial fourier along E2 ... ");
    //    //if ( GT_ABS(matrix_size_encoding_[2]/workOrder_.acceFactorE2_ - std::floor(matrix_size_encoding_[2]/workOrder_.acceFactorE2_)) > FLT_EPSILON )
    //    //{
    //    //    GADGET_WARN_MSG("matrix_size_[2] is not multiplied by acceFactorE2_ ... ");
    //    //    matrix_size_[2] = (std::floor(matrix_size_[2]/workOrder_.acceFactorE2_)+1)*workOrder_.acceFactorE2_;
    //    //}

    //    if ( 2*workOrder_.kSpaceCenterEncode2_ > (matrix_size_encoding_[2]+1) )
    //    {
    //        space_matrix_offset_E2_ = 0;

    //        workOrder_.start_E2_ = 0;
    //        workOrder_.end_E2_ = matrix_size_encoding_[2];
    //    }
    //    else
    //    {
    //        space_matrix_offset_E2_ = space_size_[2] - matrix_size_encoding_[2];

    //        workOrder_.start_E2_ = space_matrix_offset_E2_;
    //        workOrder_.end_E2_ = workOrder_.kSpaceMaxEncode2_;
    //    }
    //}
    //else
    //{
    //    space_matrix_offset_E2_ = 0;
    //}

    // ---------------------------------------------------------------------------------------------------------
    // encoding limits

    if ( GT_ABS(2*field_of_view_recon_[0]-field_of_view_encoding_[0]) < 1.0 )
    {
        meas_max_ro_ = e_space.matrixSize().x()/2;
    }
    else
    {
        meas_max_ro_ = r_space.matrixSize().x();
    }

    if (e_limits.kspace_encoding_step_1().present()) 
    {
        meas_max_idx_.kspace_encode_step_1 = matrix_size_encoding_[1]-1; // e_limits.kspace_encoding_step_1().get().maximum();
    }
    else
    {
        meas_max_idx_.kspace_encode_step_1 = 0;
        std::cout << "Setting number of kspace_encode_step_1 to 0" << std::endl;
        return GADGET_FAIL;
    }

    space_size_[0] = meas_max_ro_;

    if (e_limits.set().present())
    {
        if ( e_limits.set().get().maximum() > 0 )
            meas_max_idx_.set = e_limits.set().get().maximum() - 1;
        else
            meas_max_idx_.set = 0;

        if ( meas_max_idx_.set < 0 ) meas_max_idx_.set = 0;
    }
    else
    {
        meas_max_idx_.set = 0;
    }

    if (e_limits.phase().present())
    {
        if ( e_limits.phase().get().maximum() > 0 )
            meas_max_idx_.phase = e_limits.phase().get().maximum()-1;
        else
            meas_max_idx_.phase = 0;

        if ( meas_max_idx_.phase < 0 ) meas_max_idx_.phase = 0;
    }
    else
    {
        meas_max_idx_.phase = 0;
    }

    if (e_limits.kspace_encoding_step_2().present())
    {
        meas_max_idx_.kspace_encode_step_2 = matrix_size_encoding_[2]-1; // e_limits.kspace_encoding_step_2().get().maximum();
    }
    else
    {
        meas_max_idx_.kspace_encode_step_2 = 0;
    }

    if (e_limits.contrast().present())
    {
        if ( e_limits.contrast().get().maximum() > 0 )
            meas_max_idx_.contrast = e_limits.contrast().get().maximum()-1;
        else
            meas_max_idx_.contrast = 0;

        if ( meas_max_idx_.contrast < 0 ) meas_max_idx_.contrast = 0;
    }
    else
    {
        meas_max_idx_.contrast = 0;
    }

    if (e_limits.slice().present())
    {
        meas_max_idx_.slice = e_limits.slice().get().maximum();
    }
    else
    {
        meas_max_idx_.slice = 0;
    }

    if (e_limits.repetition().present())
    {
        meas_max_idx_.repetition = e_limits.repetition().get().maximum();
    }
    else
    {
        meas_max_idx_.repetition = 0;
    }

    if (e_limits.average().present())
    {
        meas_max_idx_.average = e_limits.average().get().maximum()-1;
    }
    else
    {
        meas_max_idx_.average = 0;
    }

    if (e_limits.segment().present())
    {
        // meas_max_idx_.segment = e_limits.segment().get().maximum()-1;
        meas_max_idx_.segment = 0;
    }
    else
    {
        meas_max_idx_.segment = 0;
    }

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
        GADGET_DEBUG1("Failed check readout status\n");
        return GADGET_FAIL;
    }

    size_t scan_counter = m1->getObjectPtr()->scan_counter;

    if ( scan_counter%1000 == 0 )
    {
        GADGET_CONDITION_MSG(verboseMode_, "--> receive scan : " << scan_counter);
    }

    // combine the segmentes
    m1->getObjectPtr()->idx.segment = 0;

    if ( (bIsNavigator || bIsRTFeedback || bIsHPFeedback || bIsDummyScan) && !bIsKSpace && !bIsRef )
    {
        m1->release();
        return GADGET_OK;
    }

    if ( !bIsRTFeedback && bIsKSpace && first_kspace_scan_ && m1->getObjectPtr()->center_sample>0 )
    {
        if ( (workOrder_.start_RO_<0) && (workOrder_.end_RO_<0) )
        {
            gtPlus_util_.findStartEndROAfterZeroFilling(m1->getObjectPtr()->center_sample, m1->getObjectPtr()->number_of_samples, workOrder_.start_RO_, workOrder_.end_RO_);

            GADGET_CONDITION_MSG(verboseMode_, "start_RO : " << workOrder_.start_RO_);
            GADGET_CONDITION_MSG(verboseMode_, "end_RO : " << workOrder_.end_RO_);

            workOrder_.kSpaceCenterRO_ = m1->getObjectPtr()->center_sample;
            workOrder_.kSpaceMaxRO_ = m1->getObjectPtr()->number_of_samples;
        }

        // if partial fourier or asymmetric echo is used, correct the kSpaceCenter
        if ( space_size_[1]-max_sampled_E1_ > workOrder_.acceFactorE1_ )
        {
            GADGET_CONDITION_MSG(verboseMode_, "Partial fourier along E1 ... ");

            if ( m1->getObjectPtr()->idx.user[5]>0 && GT_ABS(m1->getObjectPtr()->idx.user[5] - space_size_[1]/2 )<2 )
            {
                workOrder_.kSpaceCenterEncode1_ = m1->getObjectPtr()->idx.user[5];
            }

            if ( 2*workOrder_.kSpaceCenterEncode1_ > (max_sampled_E1_+1) )
            {
                space_matrix_offset_E1_ = 0;

                workOrder_.start_E1_ = 0;
                workOrder_.end_E1_ = max_sampled_E1_;
            }
            else
            {
                space_matrix_offset_E1_ = space_size_[1] - max_sampled_E1_ -1;

                workOrder_.start_E1_ = space_matrix_offset_E1_;
                workOrder_.end_E1_ = workOrder_.kSpaceMaxEncode1_;
            }
        }
        else
        {
            space_matrix_offset_E1_ = 0;
        }

        if ( space_size_[2]-max_sampled_E2_ > workOrder_.acceFactorE2_ )
        {
            GADGET_CONDITION_MSG(verboseMode_, "Partial fourier along E2 ... ");

            if ( m1->getObjectPtr()->idx.user[6]>0 && GT_ABS(m1->getObjectPtr()->idx.user[6] - space_size_[2]/2 )<2 )
            {
                workOrder_.kSpaceCenterEncode2_ = m1->getObjectPtr()->idx.user[6];
            }

            if ( 2*workOrder_.kSpaceCenterEncode2_ > (max_sampled_E2_+1) )
            {
                space_matrix_offset_E2_ = 0;

                workOrder_.start_E2_ = 0;
                workOrder_.end_E2_ = max_sampled_E2_;
            }
            else
            {
                space_matrix_offset_E2_ = space_size_[2] - max_sampled_E2_-1;

                workOrder_.start_E2_ = space_matrix_offset_E2_;
                workOrder_.end_E2_ = workOrder_.kSpaceMaxEncode2_;
            }
        }
        else
        {
            space_matrix_offset_E2_ = 0;
        }

        first_kspace_scan_ = false;
    }

    // hack for UCL data
    //if ( bIsKSpace && bIsRef )
    //{
    //    if ( m1->getObjectPtr()->idx.kspace_encode_step_1%2 == 1 )
    //    {
    //        bIsKSpace = false;
    //    }
    //}

    // store kspace read out
    if ( bIsKSpace )
    {
        if ( !storeImageData(m1, m2, bIsReflect) )
        {
            GADGET_DEBUG1("Failed check readout status\n");
            return GADGET_FAIL;
        }
    }

    // store ref read out
    if ( bIsRef )
    {
        if ( !storeRefData(m1, m2, bIsReflect) )
        {
            GADGET_DEBUG1("Failed check readout status\n");
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
        GADGET_DEBUG1("Failed triggerWorkOrder(m1)\n");
        return GADGET_FAIL;
    }

    m1->release();
    return GADGET_OK;
}

bool GtPlusAccumulatorWorkOrderTriggerGadget::needTriggerWorkOrderAllInClose()
{
    // already triggered for last acquisition
    if ( triggered_in_process_last_acq_ ) return false;

    // if never triggered in process(...)
    if ( !triggered_in_process_ && !triggered_in_process_last_acq_ ) return true;

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
    else if ( (workOrder_.CalibMode_ == ISMRMRD_noacceleration) )
    {
        return ((noacceleration_triggerDim1_==DIM_NONE)&&(noacceleration_triggerDim2_==DIM_NONE));
    }
    else
    {
        GADGET_ERROR_MSG("Unsupported calibration mode : " << workOrder_.CalibMode_);
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
        GADGET_ERROR_MSG("Unsupported calibration mode : " << workOrder_.CalibMode_);
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
            Gadgetron::gtPlus::ISMRMRDDIM& triggerDim1_, 
            Gadgetron::gtPlus::ISMRMRDDIM& triggerDim2_,
            int numOfKSpace_triggerDim1_)
{
    //bool is_first_acq_in_slice = ISMRMRD::FlagBit(ISMRMRD::ACQ_FIRST_IN_SLICE).isSet(m1->getObjectPtr()->flags);
    //if ( !is_first_acq_in_slice ) return true;

    bool is_last_acq = ((ISMRMRD::FlagBit(ISMRMRD::ACQ_LAST_IN_REPETITION).isSet(m1->getObjectPtr()->flags)) 
                    || (ISMRMRD::FlagBit(ISMRMRD::ACQ_LAST_IN_SLICE).isSet(m1->getObjectPtr()->flags)) ) 
                    && (m1->getObjectPtr()->idx.repetition==meas_max_idx_.repetition)
                    && (m1->getObjectPtr()->idx.slice==meas_max_idx_.slice)
                    && (m1->getObjectPtr()->idx.set==meas_max_idx_.set)
                    && (m1->getObjectPtr()->idx.contrast==meas_max_idx_.contrast)
                    && (m1->getObjectPtr()->idx.phase==meas_max_idx_.phase);

    curr_dim1_ = getDimValue(*(m1->getObjectPtr()), triggerDim1_);
    curr_dim2_ = getDimValue(*(m1->getObjectPtr()), triggerDim2_);

    if ( is_last_acq 
            && ( (triggerDim1_!=DIM_NONE) || (triggerDim2_!=DIM_NONE) ) )
    {
        GADGET_CONDITION_MSG(true, "Last scan in measurement - " << gtPlusISMRMRDReconUtil<ValueType>().getISMRMRDDimName(triggerDim1_ ) << " = " << curr_dim1_ << " - " << gtPlusISMRMRDReconUtil<ValueType>().getISMRMRDDimName(triggerDim2_ ) << " = " << curr_dim2_);

        if ( curr_dim1_==0 && curr_dim2_== 0 )
        {
            GADGET_CONDITION_MSG(true, "Last scan in measurement - not trigger ... ");
            return true;
        }

        triggered_in_process_last_acq_ = true;
        GADGET_CONDITION_MSG(true, "Last scan in measurement - triggered_in_process_last_acq_ : " << triggered_in_process_last_acq_);

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
            GADGET_ERROR_MSG("Unsupported calibration mode : " << workOrder_.CalibMode_);
            return false;
        }

        return true;
    }

    if ( prev_dim1_ == -1 )
    {
        prev_dim1_ = curr_dim1_;
        count_dim1_ = 0;
        GADGET_CONDITION_MSG(verboseMode_, "Current Dim1 : " << gtPlusISMRMRDReconUtil<ValueType>().getISMRMRDDimName(triggerDim1_ ) << " = " << curr_dim1_);
    }

    if ( prev_dim2_ == -1 )
    {
        prev_dim2_ = curr_dim2_;
        count_dim1_ = 0;
        GADGET_CONDITION_MSG(verboseMode_, "Current Dim2 : " << gtPlusISMRMRDReconUtil<ValueType>().getISMRMRDDimName(triggerDim2_ ) << " = " << curr_dim2_);
    }

    if ( prev_acq_header_.measurement_uid == 0 ) prev_acq_header_ = *(m1->getObjectPtr());

    bool workFlow_BufferKernel_ = false;
    bool workFlow_use_BufferedKernel_ = false;

    if ( prev_dim1_ != curr_dim1_ )
    {
        count_dim1_++;
        GADGET_CONDITION_MSG(verboseMode_, "Current Dim1 : " << gtPlusISMRMRDReconUtil<ValueType>().getISMRMRDDimName(triggerDim1_ ) << " = " << curr_dim1_);
        GADGET_CONDITION_MSG(verboseMode_, "Current Dim2 : " << gtPlusISMRMRDReconUtil<ValueType>().getISMRMRDDimName(triggerDim2_ ) << " = " << curr_dim2_);
        GADGET_CONDITION_MSG(verboseMode_, "count_dim1_ : " << count_dim1_);
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
        numOfAcquiredKSpaceForTriggerDim1 = numOfKSpace_triggerDim1_ * workOrder_.acceFactorE1_ * workOrder_.acceFactorE2_;
    }

    // trigger whenever the Dim2 is changed
    if (  triggerDim1_==DIM_NONE && triggerDim2_!=DIM_NONE  )
    {
        prev_dim1_ = curr_dim1_;
        prev_acq_header_ = *(m1->getObjectPtr());

        int prev_dim2_local_ = prev_dim2_;
        prev_dim2_ = curr_dim2_;

        if ( curr_dim2_!= prev_dim2_local_ )
        {
            GADGET_CONDITION_MSG(verboseMode_, "Trigger Dim2 : " << gtPlusISMRMRDReconUtil<ValueType>().getISMRMRDDimName(triggerDim2_ ) << " = " << prev_dim2_local_);
            GADGET_CHECK_RETURN_FALSE(triggerByDimEqual(triggerDim2_, prev_dim2_local_, workFlow_BufferKernel_, workFlow_use_BufferedKernel_));
            triggered_in_process_ = true;
        }
    }

    // trigger whenever the Dim1 is changed
    if (  triggerDim1_!=DIM_NONE && triggerDim2_==DIM_NONE  )
    {
        prev_dim2_ = curr_dim2_;

        int prev_dim1_local_ = prev_dim1_;
        prev_dim1_ = curr_dim1_;

        if ( numOfKSpace_triggerDim1_ > 0 )
        {
            if ( curr_dim1_!= prev_dim1_local_ )
            {
                if ( resetTriggerStatus(m1) )
                {
                    count_dim1_ = 0;
                    GADGET_CONDITION_MSG(verboseMode_, "Trigger Dim1 : " << gtPlusISMRMRDReconUtil<ValueType>().getISMRMRDDimName(triggerDim1_ ) << " = " << prev_dim1_local_);

                    workFlow_BufferKernel_ = false;
                    workFlow_use_BufferedKernel_ = true;
                    GADGET_CHECK_RETURN_FALSE(triggerByDimEqual(triggerDim1_, prev_dim1_local_, workFlow_BufferKernel_, workFlow_use_BufferedKernel_));
                    triggered_in_process_ = true;
                }

                if ( count_dim1_ == numOfAcquiredKSpaceForTriggerDim1 )
                {
                    GADGET_CONDITION_MSG(verboseMode_, "Trigger Dim1 : " << gtPlusISMRMRDReconUtil<ValueType>().getISMRMRDDimName(triggerDim1_ ) << " = " << prev_dim1_local_);

                    workFlow_BufferKernel_ = true;
                    workFlow_use_BufferedKernel_ = false;
                    GADGET_CHECK_RETURN_FALSE(triggerByDimLessEqual(triggerDim1_, prev_dim1_local_, workFlow_BufferKernel_, workFlow_use_BufferedKernel_));
                    triggered_in_process_ = true;
                }
                else if ( count_dim1_ > numOfAcquiredKSpaceForTriggerDim1 )
                {
                    GADGET_CONDITION_MSG(verboseMode_, "Trigger Dim1 : " << gtPlusISMRMRDReconUtil<ValueType>().getISMRMRDDimName(triggerDim1_ ) << " = " << prev_dim1_local_);

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
                GADGET_CONDITION_MSG(verboseMode_, "Trigger Dim1 : " << gtPlusISMRMRDReconUtil<ValueType>().getISMRMRDDimName(triggerDim1_ ) << " = " << prev_dim1_local_);
                GADGET_CHECK_RETURN_FALSE(triggerByDimEqual(triggerDim1_, prev_dim1_local_, workFlow_BufferKernel_, workFlow_use_BufferedKernel_));
                triggered_in_process_ = true;
            }

            prev_acq_header_ = *(m1->getObjectPtr());
        }
    }

    if (  triggerDim1_!=DIM_NONE && triggerDim2_!=DIM_NONE  )
    {
        int prev_dim1_local_ = prev_dim1_;
        int prev_dim2_local_ = prev_dim2_;

        prev_dim1_ = curr_dim1_;
        prev_dim2_ = curr_dim2_;

        if ( numOfKSpace_triggerDim1_ > 0 )
        {
            if ( (curr_dim2_!=prev_dim2_local_) || resetTriggerStatus(m1) )
            {
                count_dim1_ = 0;
                GADGET_CONDITION_MSG(verboseMode_, "Trigger Dim1 : " << gtPlusISMRMRDReconUtil<ValueType>().getISMRMRDDimName(triggerDim1_ ) << " = " << prev_dim1_local_ 
                    << "; Dim2 : " << gtPlusISMRMRDReconUtil<ValueType>().getISMRMRDDimName(triggerDim2_ ) << " = " << prev_dim2_local_);

                workFlow_BufferKernel_ = false;
                workFlow_use_BufferedKernel_ = true;

                GADGET_CHECK_RETURN_FALSE(triggerByDimEqual(triggerDim1_, prev_dim1_local_, triggerDim2_, prev_dim2_local_, workFlow_BufferKernel_, workFlow_use_BufferedKernel_));

                triggered_in_process_ = true;
            }

            if (curr_dim1_!=prev_dim1_local_)
            {
                if ( count_dim1_ == numOfAcquiredKSpaceForTriggerDim1 )
                {
                    GADGET_CONDITION_MSG(verboseMode_, "Trigger Dim1 : " << gtPlusISMRMRDReconUtil<ValueType>().getISMRMRDDimName(triggerDim1_ ) << " = " << prev_dim1_local_ 
                        << "; Dim2 : " << gtPlusISMRMRDReconUtil<ValueType>().getISMRMRDDimName(triggerDim2_ ) << " = " << prev_dim2_local_);

                    workFlow_BufferKernel_ = true;
                    workFlow_use_BufferedKernel_ = false;
                    GADGET_CHECK_RETURN_FALSE(triggerByDim1LessEqualDim2Equal(triggerDim1_, prev_dim1_local_, triggerDim2_, prev_dim2_local_, workFlow_BufferKernel_, workFlow_use_BufferedKernel_));
                    triggered_in_process_ = true;
                }
                else if ( count_dim1_ > numOfAcquiredKSpaceForTriggerDim1 )
                {
                    GADGET_CONDITION_MSG(verboseMode_, "Trigger Dim1 : " << gtPlusISMRMRDReconUtil<ValueType>().getISMRMRDDimName(triggerDim1_ ) << " = " << prev_dim1_local_ 
                        << "; Dim2 : " << gtPlusISMRMRDReconUtil<ValueType>().getISMRMRDDimName(triggerDim2_ ) << " = " << prev_dim2_local_);

                    workFlow_BufferKernel_ = false;
                    workFlow_use_BufferedKernel_ = true;
                    GADGET_CHECK_RETURN_FALSE(triggerByDimEqual(triggerDim1_, prev_dim1_local_, triggerDim2_, prev_dim2_local_, workFlow_BufferKernel_, workFlow_use_BufferedKernel_));
                    triggered_in_process_ = true;
                }
            }

            prev_acq_header_ = *(m1->getObjectPtr());
        }
        else
        {
            // trigger whenever the Dim2 is changed
            if ( curr_dim2_!= prev_dim2_local_ )
            {
                GADGET_CONDITION_MSG(verboseMode_, "Trigger Dim2 : " << gtPlusISMRMRDReconUtil<ValueType>().getISMRMRDDimName(triggerDim2_ ) << " = " << prev_dim2_local_);
                GADGET_CHECK_RETURN_FALSE(triggerByDimEqual(triggerDim2_, prev_dim2_local_, workFlow_BufferKernel_, workFlow_use_BufferedKernel_));
                triggered_in_process_ = true;
            }

            prev_acq_header_ = *(m1->getObjectPtr());
        }
    }

    return true;
}

bool GtPlusAccumulatorWorkOrderTriggerGadget::
triggerWorkOrderLastCountInClose(Gadgetron::gtPlus::ISMRMRDDIM& triggerDim1_, Gadgetron::gtPlus::ISMRMRDDIM& triggerDim2_, int numOfKSpace_triggerDim1_)
{
    GADGET_CONDITION_MSG(verboseMode_, "Current Dim1 InClose : " << gtPlusISMRMRDReconUtil<ValueType>().getISMRMRDDimName(triggerDim1_ ) << " = " << curr_dim1_);
    GADGET_CONDITION_MSG(verboseMode_, "Current Dim2 InClose : " << gtPlusISMRMRDReconUtil<ValueType>().getISMRMRDDimName(triggerDim2_ ) << " = " << curr_dim2_);

    if ( prev_dim1_ != curr_dim1_ )
    {
        count_dim1_++;
    }

    bool workFlow_BufferKernel_ = false;
    bool workFlow_use_BufferedKernel_ = false;

    int numOfAcquiredKSpaceForTriggerDim1 = numOfKSpace_triggerDim1_;
    if ( workOrder_.CalibMode_ == ISMRMRD_interleaved )
    {
        numOfAcquiredKSpaceForTriggerDim1 = numOfKSpace_triggerDim1_ * workOrder_.acceFactorE1_ * workOrder_.acceFactorE2_;
    }

    int prev_dim1_local_ = prev_dim1_;
    int prev_dim2_local_ = prev_dim2_;

    prev_dim1_ = curr_dim1_;
    prev_dim2_ = curr_dim2_;

    // trigger whenever the Dim2 is changed
    if (  triggerDim1_==DIM_NONE && triggerDim2_!=DIM_NONE  )
    {
        GADGET_CONDITION_MSG(verboseMode_, "Trigger Dim2 : " << gtPlusISMRMRDReconUtil<ValueType>().getISMRMRDDimName(triggerDim2_ ) << " = " << prev_dim2_local_);
        GADGET_CHECK_RETURN_FALSE(triggerByDimEqual(triggerDim2_, prev_dim2_local_, workFlow_BufferKernel_, workFlow_use_BufferedKernel_));
    }

    // trigger whenever the Dim1 is changed
    if (  triggerDim1_!=DIM_NONE && triggerDim2_==DIM_NONE  )
    {
        if ( numOfKSpace_triggerDim1_ > 0 )
        {
            if ( count_dim1_ <= numOfAcquiredKSpaceForTriggerDim1 )
            {
                GADGET_CONDITION_MSG(verboseMode_, "Trigger Dim1 : " << gtPlusISMRMRDReconUtil<ValueType>().getISMRMRDDimName(triggerDim1_ ) << " <= " << prev_dim1_local_);
                workFlow_BufferKernel_ = true;
                workFlow_use_BufferedKernel_ = false;
                GADGET_CHECK_RETURN_FALSE(triggerByDimLessEqual(triggerDim1_, prev_dim1_local_, workFlow_BufferKernel_, workFlow_use_BufferedKernel_));
            }
            else if ( count_dim1_ > numOfAcquiredKSpaceForTriggerDim1 )
            {
                GADGET_CONDITION_MSG(verboseMode_, "Trigger Dim1 : " << gtPlusISMRMRDReconUtil<ValueType>().getISMRMRDDimName(triggerDim1_ ) << " = " << prev_dim1_local_);
                workFlow_BufferKernel_ = false;
                workFlow_use_BufferedKernel_ = true;
                GADGET_CHECK_RETURN_FALSE(triggerByDimEqual(triggerDim1_, prev_dim1_local_, workFlow_BufferKernel_, workFlow_use_BufferedKernel_));
            }
        }
        else
        {
            GADGET_CONDITION_MSG(verboseMode_, "Trigger Dim1 : " << gtPlusISMRMRDReconUtil<ValueType>().getISMRMRDDimName(triggerDim1_ ) << " = " << prev_dim1_local_);
            GADGET_CHECK_RETURN_FALSE(triggerByDimEqual(triggerDim1_, prev_dim1_local_, workFlow_BufferKernel_, workFlow_use_BufferedKernel_));
        }
    }

    if (  triggerDim1_!=DIM_NONE && triggerDim2_!=DIM_NONE  )
    {
        if ( numOfKSpace_triggerDim1_ > 0 )
        {
            if ( count_dim1_ <= numOfAcquiredKSpaceForTriggerDim1 ) // no more data will be available, so have to do the recon
            {
                GADGET_CONDITION_MSG(verboseMode_, "Trigger Dim1 : " << gtPlusISMRMRDReconUtil<ValueType>().getISMRMRDDimName(triggerDim1_ ) << " <= " << prev_dim1_local_);
                workFlow_BufferKernel_ = true;
                workFlow_use_BufferedKernel_ = false;
                GADGET_CHECK_RETURN_FALSE(triggerByDim1LessEqualDim2Equal(triggerDim1_, prev_dim1_local_, triggerDim2_, prev_dim2_local_, workFlow_BufferKernel_, workFlow_use_BufferedKernel_));
            }
            else if ( count_dim1_ > numOfAcquiredKSpaceForTriggerDim1 )
            {
                GADGET_CONDITION_MSG(verboseMode_, "Trigger Dim1 : " << gtPlusISMRMRDReconUtil<ValueType>().getISMRMRDDimName(triggerDim1_ ) << " = " << prev_dim1_local_);
                workFlow_BufferKernel_ = false;
                workFlow_use_BufferedKernel_ = true;
                GADGET_CHECK_RETURN_FALSE(triggerByDimEqual(triggerDim1_, prev_dim1_local_, triggerDim2_, prev_dim2_local_, workFlow_BufferKernel_, workFlow_use_BufferedKernel_));
            }
        }
        else
        {
            GADGET_CONDITION_MSG(verboseMode_, "Trigger Dim2 : " << gtPlusISMRMRDReconUtil<ValueType>().getISMRMRDDimName(triggerDim2_ ) << " = " << prev_dim2_local_);
            GADGET_CHECK_RETURN_FALSE(triggerByDimEqual(triggerDim2_, prev_dim2_local_, workFlow_BufferKernel_, workFlow_use_BufferedKernel_));
        }
    }

    return true;
}

bool GtPlusAccumulatorWorkOrderTriggerGadget::checkStatus(uint64_t flag, int samples, 
    bool& bIsKSpace, bool& bIsRef, bool& bIsNoise, bool& bIsPhaseCorr, bool& bIsReflect, bool& bIsOther,
    bool& bIsNavigator, bool& bIsRTFeedback, bool& bIsHPFeedback, bool& bIsDummyScan)
{
    bIsNoise = ISMRMRD::FlagBit(ISMRMRD::ACQ_IS_NOISE_MEASUREMENT).isSet(flag);
    bool is_ref = ISMRMRD::FlagBit(ISMRMRD::ACQ_IS_PARALLEL_CALIBRATION).isSet(flag);
    bool is_ref_kspace = ISMRMRD::FlagBit(ISMRMRD::ACQ_IS_PARALLEL_CALIBRATION_AND_IMAGING).isSet(flag);
    bIsReflect = ISMRMRD::FlagBit(ISMRMRD::ACQ_IS_REVERSE).isSet(flag);
    bIsPhaseCorr = ISMRMRD::FlagBit(ISMRMRD::ACQ_IS_PHASECORR_DATA).isSet(flag);
    bIsNavigator = ISMRMRD::FlagBit(ISMRMRD::ACQ_IS_NAVIGATION_DATA).isSet(flag);
    bIsRTFeedback = ISMRMRD::FlagBit(ISMRMRD::ACQ_IS_RTFEEDBACK_DATA).isSet(flag);
    bIsHPFeedback = ISMRMRD::FlagBit(ISMRMRD::ACQ_IS_HPFEEDBACK_DATA).isSet(flag);
    bIsDummyScan = ISMRMRD::FlagBit(ISMRMRD::ACQ_IS_DUMMYSCAN_DATA).isSet(flag);

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

        idx.segment = 0; // combine the segments

        if ( workOrder_.data_.get_number_of_elements() <= 0 )
        {
            meas_max_channel_ = m1->getObjectPtr()->active_channels;

            int E1 = workOrder_.kSpaceMaxEncode1_+1;
            int E2 = workOrder_.kSpaceMaxEncode2_+1;
            if ( E2 == 0 ) E2 = 1;

            if ( E1 < matrix_size_encoding_[1] ) E1 = matrix_size_encoding_[1];
            if ( E2 < matrix_size_encoding_[2] ) E2 = matrix_size_encoding_[2];

            // find the loop counter boundary and allocate the buffer
            GADGET_CONDITION_MSG(verboseMode_, "[RO E1 Cha Slice E2 Con Phase Rep Set Seg] = [" 
                               << meas_max_ro_ 
                               << " " << E1 
                               << " " << meas_max_channel_ 
                               << " " << meas_max_idx_.slice+1 
                               << " " << E2 
                               << " " << meas_max_idx_.contrast+1 
                               << " " << meas_max_idx_.phase+1 
                               << " " << meas_max_idx_.repetition+1 
                               << " " << meas_max_idx_.set+1 
                               << " " << meas_max_idx_.segment+1 << "]");

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

            size_t N = dimensions_.size();
            for ( ii=0; ii<N; ii++ )
            {
                GADGET_CONDITION_MSG(verboseMode_, "dimensions_[" << ii << "] = " << dimensions_[ii]);
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
            }
            catch(...)
            {
                GADGET_DEBUG1("Failed create buffer\n");
                return false;
            }

            // allocate message buffer
            int matrix_size[10];
            for ( ii=0; ii<10; ii++ )
            {
                matrix_size[ii] = dimensions_[ii];
            }

            if (!(messageImage_ = new GtPlusGadgetImageArray(matrix_size))) 
            {
                GADGET_DEBUG1("Failed create buffer\n");
                return false;
            }
        }

        // if necessary, shift the E1/E2 indexes
        //if ( workOrder_.start_E1_ > 0 )
        //{
        //    idx.kspace_encode_step_1 += workOrder_.start_E1_;
        //}

        //if ( workOrder_.start_E2_ > 0 )
        //{
        //    idx.kspace_encode_step_2 += workOrder_.start_E2_;
        //}

        std::complex<float>* b = workOrder_.data_.begin();
        std::complex<float>* d = m2->getObjectPtr()->get_data_ptr();
        if (samples != static_cast<int>(dimensions_[0])) 
        {
            GADGET_DEBUG1("Wrong number of samples received\n");
            return false;
        }

        //Copy the data for all the channels
        hoNDArray<std::complex<float> > reflectBuf;
        if ( isReflect )
        {
            reflectBuf.create(samples);
        }

        std::vector<size_t> pos(10);
        for (int c = 0; c < m1->getObjectPtr()->active_channels; c++) 
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
            long long offsetBuffer = workOrder_.data_.calculate_offset(pos);

            if ( isReflect )
            {
                for ( int s=0; s<samples; s++ )
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
            GADGET_DEBUG1("Failed in fillImageInfo(m1, messageImage_, idx)\n");
            return false;
        }
    }
    catch(...)
    {
        GADGET_DEBUG1("Errors in GtPlusAccumulatorWorkOrderTriggerGadget::storeImageData(...) ... \n");
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

        idx.segment = 0; // combine the segments

        if ( workOrder_.ref_.get_number_of_elements() <= 0 )
        {
            meas_max_channel_ = m1->getObjectPtr()->active_channels;

            int E1 = workOrder_.kSpaceMaxEncode1_+1;
            int E2 = workOrder_.kSpaceMaxEncode2_+1;
            if ( E2 == 0 ) E2 = 1;

            if ( E1 < matrix_size_encoding_[1] ) E1 = matrix_size_encoding_[1];
            if ( E2 < matrix_size_encoding_[2] ) E2 = matrix_size_encoding_[2];

            size_t RO = meas_max_ro_;

            if ( (samples < meas_max_ro_) 
                && (( workOrder_.CalibMode_==ISMRMRD_separate || workOrder_.CalibMode_==ISMRMRD_external )) )
            {
                RO = samples;
            }

            // find the loop counter boundary and allocate the buffer
            GADGET_CONDITION_MSG(verboseMode_, "[RO E1 Cha Slice E2 Con Phase Rep Set Seg] = [" 
                               << RO 
                               << " " << E1 
                               << " " << meas_max_channel_ 
                               << " " << meas_max_idx_.slice+1 
                               << " " << E2 
                               << " " << meas_max_idx_.contrast+1 
                               << " " << meas_max_idx_.phase+1 
                               << " " << meas_max_idx_.repetition+1 
                               << " " << meas_max_idx_.set+1 
                               << " " << meas_max_idx_.segment+1 << "]");

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

            size_t N = dimensions_.size();
            for ( ii=0; ii<N; ii++ )
            {
                GADGET_CONDITION_MSG(verboseMode_, "ref dimensions_[" << ii << "] = " << dimensions_[ii]);
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
                GADGET_DEBUG1("Failed create ref buffer\n");
                return false;
            }
        }

        // if necessary, shift the E1/E2 indexes
        //if ( workOrder_.CalibMode_ == ISMRMRD_embedded )
        //{
        //    if ( workOrder_.start_E1_ > 0 )
        //    {
        //        idx.kspace_encode_step_1 += workOrder_.start_E1_;
        //    }

        //    if ( workOrder_.start_E2_ > 0 )
        //    {
        //        idx.kspace_encode_step_2 += workOrder_.start_E2_;
        //    }
        //}

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

            size_t ii;
            for ( ii=0; ii<ISMRMRD_USER_INTS; ii++ )
            {
                if ( idx.user[ii] > meas_max_idx_ref_.user[ii] ) meas_max_idx_ref_.user[ii] = idx.user[ii];
            }
        }

        std::complex<float>* b = workOrder_.ref_.begin();
        std::complex<float>* d = m2->getObjectPtr()->get_data_ptr();
        if (samples != static_cast<int>(dimensions_[0])) 
        {
            GADGET_DEBUG1("Wrong number of samples received\n");
            return false;
        }

        //Copy the data for all the channels
        hoNDArray<std::complex<float> > reflectBuf;
        if ( isReflect )
        {
            reflectBuf.create(samples);
        }

        std::vector<size_t> pos(10);
        for (int c = 0; c < m1->getObjectPtr()->active_channels; c++) 
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
            long long offsetBuffer = workOrder_.ref_.calculate_offset(pos);

            if ( isReflect )
            {
                for ( int s=0; s<samples; s++ )
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
    }
    catch(...)
    {
        GADGET_DEBUG1("Errors in GtPlusAccumulatorWorkOrderTriggerGadget::storeRefData(...) ... \n");
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
        int max_channel = 0;
        int max_col = 0;

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

            if ( readOutBuffer[a].acqHead_.active_channels > max_channel ) 
                max_channel = readOutBuffer[a].acqHead_.active_channels;
        }

        GADGET_CONDITION_MSG(verboseMode_, "[RO E1 Cha Slice E2 Contrast Phase Rep Set Seg] = [" 
                               << max_col 
                               << " " << max_idx.kspace_encode_step_1+1 
                               << " " << max_channel 
                               << " " << max_idx.slice+1 
                               << " " << max_idx.kspace_encode_step_2+1 
                               << " " << max_idx.contrast+1 
                               << " " << max_idx.phase+1 
                               << " " << max_idx.repetition+1 
                               << " " << max_idx.set+1 
                               << " " << max_idx.segment+1 << "]");

        // alloate buffer for data
        std::vector<size_t> dims(10);
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
            GADGET_DEBUG1("Failed create buffer\n");
            return false;
        }

        std::complex<float>* b = buf.begin();

        // copy the data
        int c;
        std::vector<size_t> pos(10);

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
        GADGET_DEBUG1("Errors in GtPlusAccumulatorWorkOrderTriggerGadget::fillBuffer(...) ... \n");
        return false;
    }

    return true;
}

bool GtPlusAccumulatorWorkOrderTriggerGadget::fillImageInfo(GadgetContainerMessage<ISMRMRD::AcquisitionHeader>* m1, GtPlusGadgetImageArray* messageImage, const ISMRMRD::EncodingCounters& idx)
{
    try
    {
        // fill the message info
        int offset = messageImage->get_offset(idx.slice, idx.kspace_encode_step_2, idx.contrast, idx.phase, idx.repetition, idx.set, idx.segment);

        // if it is the first acq in a slice, fill in all information
        bool is_first_acq_in_slice = ISMRMRD::FlagBit(ISMRMRD::ACQ_FIRST_IN_SLICE).isSet(m1->getObjectPtr()->flags);

        /*if ( is_first_acq_in_slice 
            || ( messageImage->imageArray_[offset].version==0 
                    && messageImage->imageArray_[offset].flags==0 
                    && messageImage->imageArray_[offset].measurement_uid==0 ) )*/
        if ( messageImage->imageArray_[offset].version==0 
                    && messageImage->imageArray_[offset].flags==0 
                    && messageImage->imageArray_[offset].measurement_uid==0 )
        {
            GADGET_CONDITION_MSG(verboseMode_, "--> buffer image header - offset = " << offset << " - [SLC E2 CON PHS REP SET] = [" 
                                                                      << idx.slice << " " 
                                                                      << idx.kspace_encode_step_2 << " " 
                                                                      << idx.contrast << " " 
                                                                      << idx.phase << " " 
                                                                      << idx.repetition << " " 
                                                                      << idx.set << "]");

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

            messageImage->imageArray_[offset].acquisition_time_stamp = m1->getObjectPtr()->acquisition_time_stamp;

            messageImage->imageArray_[offset].physiology_time_stamp[0] = m1->getObjectPtr()->physiology_time_stamp[0];
            messageImage->imageArray_[offset].physiology_time_stamp[1] = m1->getObjectPtr()->physiology_time_stamp[1];
            messageImage->imageArray_[offset].physiology_time_stamp[2] = m1->getObjectPtr()->physiology_time_stamp[2];

            messageImage->imageArray_[offset].image_data_type = ISMRMRD::DATA_COMPLEX_FLOAT;

            messageImage->imageArray_[offset].image_type = ISMRMRD::TYPE_MAGNITUDE;

            messageImage->imageArray_[offset].image_index = ++image_counter_;
            messageImage->imageArray_[offset].image_series_index = image_series_;

            // need to store the free user parameters
            memcpy(messageImage->imageArray_[offset].user_int, m1->getObjectPtr()->user_int, sizeof(int32_t)*8);
            memcpy(messageImage->imageArray_[offset].user_float, m1->getObjectPtr()->user_float, sizeof(float)*8);
        }

        // whether or not this acq is the first in a slice, we need to fill the TimeStamps and PMUTimeStamps
        messageImage->imageArray_[offset].time_stamps[idx.kspace_encode_step_1] = m1->getObjectPtr()->acquisition_time_stamp;
        messageImage->imageArray_[offset].pmu_time_stamps[idx.kspace_encode_step_1] = m1->getObjectPtr()->physiology_time_stamp[0];
    }
    catch(...)
    {
        GADGET_DEBUG1("Errors in GtPlusAccumulatorWorkOrderTriggerGadget::fillImageInfo(...) ... \n");
        return false;
    }

    return true;
}

size_t GtPlusAccumulatorWorkOrderTriggerGadget::
computeEncodedSizeE1(size_t centerE1, size_t maxE1)
{
    int E1;
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
    int E2;
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
triggerByDimEqual(Gadgetron::gtPlus::ISMRMRDDIM& triggerDim, size_t value, bool workFlow_BufferKernel_, bool workFlow_use_BufferedKernel_)
{
    try
    {
        GADGET_CONDITION_MSG(verboseMode_, "GtPlusAccumulatorWorkOrderTriggerGadget - triggerByDimEqual(triggerDim, value) ... ");

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
            GADGET_CONDITION_MSG(verboseMode_, "GtPlusAccumulatorWorkOrderTriggerGadget - phase correction signal found : " << phaseCorrBuffer_.size());

            if ( !fillBuffer(phaseCorrBuffer_, workOrder_.phaseCorr_, workOrder_.reflect_phaseCorr_) )
            {
                GADGET_DEBUG1("fillBuffer(phaseCorrBuffer_) failed ... \n");
                cm1->release();
                return false;
            }

            cm2->getObjectPtr()->phaseCorr_ = workOrder_.phaseCorr_;
            cm2->getObjectPtr()->reflect_phaseCorr_ = workOrder_.reflect_phaseCorr_;
        }

        if (!noiseBuffer_.empty())
        {
            GADGET_CONDITION_MSG(verboseMode_, "GtPlusAccumulatorGadget - noise signal found : " << noiseBuffer_.size());

            ReflectBufferType tmpBuf;
            if ( !fillBuffer(noiseBuffer_, workOrder_.noise_, tmpBuf) )
            {
                GADGET_DEBUG1("fillBuffer(noiseBuffer_) failed ... \n");
                cm1->release();
                return false;
            }

            cm2->getObjectPtr()->noise_ = workOrder_.noise_;
        }

        if (!otherBuffer_.empty())
        {
            GADGET_CONDITION_MSG(verboseMode_, "GtPlusAccumulatorGadget - other signal found : " << otherBuffer_.size());

            if ( !fillBuffer(otherBuffer_, workOrder_.other_, workOrder_.reflect_other_) )
            {
                GADGET_DEBUG1("fillBuffer(otherBuffer_) failed ... \n");
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
        GADGET_ERROR_MSG("Errors in GtPlusAccumulatorWorkOrderTriggerGadget::triggerByDimEqual(triggerDim, value) ... ");
        return false;
    }

    return true;
}

bool GtPlusAccumulatorWorkOrderTriggerGadget::
triggerByDimLessEqual(Gadgetron::gtPlus::ISMRMRDDIM& triggerDim, size_t value, bool workFlow_BufferKernel_, bool workFlow_use_BufferedKernel_)
{
    try
    {
        GADGET_CONDITION_MSG(verboseMode_, "GtPlusAccumulatorWorkOrderTriggerGadget - triggerByDimEqual(triggerDim, value) ... ");

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
            GADGET_CONDITION_MSG(verboseMode_, "GtPlusAccumulatorWorkOrderTriggerGadget - phase correction signal found : " << phaseCorrBuffer_.size());

            if ( !fillBuffer(phaseCorrBuffer_, workOrder_.phaseCorr_, workOrder_.reflect_phaseCorr_) )
            {
                GADGET_DEBUG1("fillBuffer(phaseCorrBuffer_) failed ... \n");
                cm1->release();
                return false;
            }

            cm2->getObjectPtr()->phaseCorr_ = workOrder_.phaseCorr_;
            cm2->getObjectPtr()->reflect_phaseCorr_ = workOrder_.reflect_phaseCorr_;
        }

        if (!noiseBuffer_.empty())
        {
            GADGET_CONDITION_MSG(verboseMode_, "GtPlusAccumulatorGadget - noise signal found : " << noiseBuffer_.size());

            ReflectBufferType tmpBuf;
            if ( !fillBuffer(noiseBuffer_, workOrder_.noise_, tmpBuf) )
            {
                GADGET_DEBUG1("fillBuffer(noiseBuffer_) failed ... \n");
                cm1->release();
                return false;
            }

            cm2->getObjectPtr()->noise_ = workOrder_.noise_;
        }

        if (!otherBuffer_.empty())
        {
            GADGET_CONDITION_MSG(verboseMode_, "GtPlusAccumulatorGadget - other signal found : " << otherBuffer_.size());

            if ( !fillBuffer(otherBuffer_, workOrder_.other_, workOrder_.reflect_other_) )
            {
                GADGET_DEBUG1("fillBuffer(otherBuffer_) failed ... \n");
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
        GADGET_ERROR_MSG("Errors in GtPlusAccumulatorWorkOrderTriggerGadget::triggerByDimLessEqual(triggerDim, value) ... ");
        return false;
    }

    return true;
}

bool GtPlusAccumulatorWorkOrderTriggerGadget::
triggerByDimEqual(Gadgetron::gtPlus::ISMRMRDDIM& triggerDim1, size_t value1, Gadgetron::gtPlus::ISMRMRDDIM& triggerDim2, size_t value2, bool workFlow_BufferKernel_, bool workFlow_use_BufferedKernel_)
{
    try
    {
        GADGET_CONDITION_MSG(verboseMode_, "GtPlusAccumulatorWorkOrderTriggerGadget - triggerByDimEqual(triggerDim1, value1, triggerDim2, value2) ... ");

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
            GADGET_CONDITION_MSG(verboseMode_, "GtPlusAccumulatorWorkOrderTriggerGadget - phase correction signal found : " << phaseCorrBuffer_.size());

            if ( !fillBuffer(phaseCorrBuffer_, workOrder_.phaseCorr_, workOrder_.reflect_phaseCorr_) )
            {
                GADGET_DEBUG1("fillBuffer(phaseCorrBuffer_) failed ... \n");
                cm1->release();
                return false;
            }

            cm2->getObjectPtr()->phaseCorr_ = workOrder_.phaseCorr_;
            cm2->getObjectPtr()->reflect_phaseCorr_ = workOrder_.reflect_phaseCorr_;
        }

        if (!noiseBuffer_.empty())
        {
            GADGET_CONDITION_MSG(verboseMode_, "GtPlusAccumulatorGadget - noise signal found : " << noiseBuffer_.size());

            ReflectBufferType tmpBuf;
            if ( !fillBuffer(noiseBuffer_, workOrder_.noise_, tmpBuf) )
            {
                GADGET_DEBUG1("fillBuffer(noiseBuffer_) failed ... \n");
                cm1->release();
                return false;
            }

            cm2->getObjectPtr()->noise_ = workOrder_.noise_;
        }

        if (!otherBuffer_.empty())
        {
            GADGET_CONDITION_MSG(verboseMode_, "GtPlusAccumulatorGadget - other signal found : " << otherBuffer_.size());

            if ( !fillBuffer(otherBuffer_, workOrder_.other_, workOrder_.reflect_other_) )
            {
                GADGET_DEBUG1("fillBuffer(otherBuffer_) failed ... \n");
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
        GADGET_ERROR_MSG("Errors in GtPlusAccumulatorWorkOrderTriggerGadget::triggerByDimEqual(triggerDim1, value1, triggerDim2, value2) ... ");
        return false;
    }

    return true;
}


bool GtPlusAccumulatorWorkOrderTriggerGadget::
triggerByDim1LessEqualDim2Equal(Gadgetron::gtPlus::ISMRMRDDIM& triggerDim1, size_t value1, Gadgetron::gtPlus::ISMRMRDDIM& triggerDim2, size_t value2, bool workFlow_BufferKernel_, bool workFlow_use_BufferedKernel_)
{
    try
    {
        GADGET_CONDITION_MSG(verboseMode_, "GtPlusAccumulatorWorkOrderTriggerGadget - triggerByDim1LessEqualDim2Equal(triggerDim1, value1, triggerDim2, value2) ... ");

        GadgetContainerMessage<GtPlusGadgetImageArray>* cm1 = new GadgetContainerMessage<GtPlusGadgetImageArray>();
        GadgetContainerMessage< WorkOrderType >* cm2 = new GadgetContainerMessage< WorkOrderType >();

        workOrder_.duplicate(*cm2->getObjectPtr());
        cm2->getObjectPtr()->workFlow_BufferKernel_ = workFlow_BufferKernel_;
        cm2->getObjectPtr()->workFlow_use_BufferedKernel_ = workFlow_use_BufferedKernel_;

        cm1->cont(cm2);

        // copy the image content
        GADGET_CHECK_RETURN_FALSE(gtPlusISMRMRDReconUtil<ValueType>().extractSubArrayForDim1LessEqualDim2Equal(workOrder_.data_, cm2->getObjectPtr()->data_, triggerDim1, value1, triggerDim2, value2));
        GADGET_CHECK_RETURN_FALSE(gtPlusISMRMRDReconUtil<unsigned short>().extractSubArrayForDim1LessEqualDim2Equal(workOrder_.reflect_, cm2->getObjectPtr()->reflect_, triggerDim1, value1, triggerDim2, value2));

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
            GADGET_CONDITION_MSG(verboseMode_, "GtPlusAccumulatorWorkOrderTriggerGadget - phase correction signal found : " << phaseCorrBuffer_.size());

            if ( !fillBuffer(phaseCorrBuffer_, workOrder_.phaseCorr_, workOrder_.reflect_phaseCorr_) )
            {
                GADGET_DEBUG1("fillBuffer(phaseCorrBuffer_) failed ... \n");
                cm1->release();
                return false;
            }

            cm2->getObjectPtr()->phaseCorr_ = workOrder_.phaseCorr_;
            cm2->getObjectPtr()->reflect_phaseCorr_ = workOrder_.reflect_phaseCorr_;
        }

        if (!noiseBuffer_.empty())
        {
            GADGET_CONDITION_MSG(verboseMode_, "GtPlusAccumulatorGadget - noise signal found : " << noiseBuffer_.size());

            ReflectBufferType tmpBuf;
            if ( !fillBuffer(noiseBuffer_, workOrder_.noise_, tmpBuf) )
            {
                GADGET_DEBUG1("fillBuffer(noiseBuffer_) failed ... \n");
                cm1->release();
                return false;
            }

            cm2->getObjectPtr()->noise_ = workOrder_.noise_;
        }

        if (!otherBuffer_.empty())
        {
            GADGET_CONDITION_MSG(verboseMode_, "GtPlusAccumulatorGadget - other signal found : " << otherBuffer_.size());

            if ( !fillBuffer(otherBuffer_, workOrder_.other_, workOrder_.reflect_other_) )
            {
                GADGET_DEBUG1("fillBuffer(otherBuffer_) failed ... \n");
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
        GADGET_ERROR_MSG("Errors in GtPlusAccumulatorWorkOrderTriggerGadget::triggerByDim1LessEqualDim2Equal(triggerDim1, value1, triggerDim2, value2) ... ");
        return false;
    }

    return true;
}

bool GtPlusAccumulatorWorkOrderTriggerGadget::triggerWorkOrderAllInClose()
{
    try
    {
        GADGET_CONDITION_MSG(verboseMode_, "GtPlusAccumulatorWorkOrderTriggerGadget - triggerWorkOrderAllInClose ... ");

        GadgetContainerMessage<GtPlusGadgetImageArray>* cm1 = new GadgetContainerMessage<GtPlusGadgetImageArray>();
        GadgetContainerMessage< WorkOrderType >* cm2 = new GadgetContainerMessage< WorkOrderType >();

        workOrder_.duplicate(*cm2->getObjectPtr());
        cm2->getObjectPtr()->workFlow_BufferKernel_ = false;
        cm2->getObjectPtr()->workFlow_use_BufferedKernel_ = false;

        cm1->cont(cm2);

        // copy the image content
        cm2->getObjectPtr()->data_ = workOrder_.data_;
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
            GADGET_CONDITION_MSG(verboseMode_, "GtPlusAccumulatorWorkOrderTriggerGadget - phase correction signal found : " << phaseCorrBuffer_.size());

            if ( !fillBuffer(phaseCorrBuffer_, workOrder_.phaseCorr_, workOrder_.reflect_phaseCorr_) )
            {
                GADGET_DEBUG1("fillBuffer(phaseCorrBuffer_) failed ... \n");
                cm1->release();
                return GADGET_FAIL;
            }

            cm2->getObjectPtr()->phaseCorr_ = workOrder_.phaseCorr_;
            cm2->getObjectPtr()->reflect_phaseCorr_ = workOrder_.reflect_phaseCorr_;
        }

        if (!noiseBuffer_.empty())
        {
            GADGET_CONDITION_MSG(verboseMode_, "GtPlusAccumulatorGadget - noise signal found : " << noiseBuffer_.size());

            ReflectBufferType tmpBuf;
            if ( !fillBuffer(noiseBuffer_, workOrder_.noise_, tmpBuf) )
            {
                GADGET_DEBUG1("fillBuffer(noiseBuffer_) failed ... \n");
                cm1->release();
                return GADGET_FAIL;
            }

            cm2->getObjectPtr()->noise_ = workOrder_.noise_;
        }

        if (!otherBuffer_.empty())
        {
            GADGET_CONDITION_MSG(verboseMode_, "GtPlusAccumulatorGadget - other signal found : " << otherBuffer_.size());

            if ( !fillBuffer(otherBuffer_, workOrder_.other_, workOrder_.reflect_other_) )
            {
                GADGET_DEBUG1("fillBuffer(otherBuffer_) failed ... \n");
                cm1->release();
                return GADGET_FAIL;
            }

            cm2->getObjectPtr()->other_ = workOrder_.other_;
            cm2->getObjectPtr()->reflect_other_ = workOrder_.reflect_other_;
        }

        // send to next gadget
        if (this->next()->putq(cm1) < 0) 
        {
            return GADGET_FAIL;
        }
    }
    catch(...)
    {
        GADGET_ERROR_MSG("Errors in GtPlusAccumulatorWorkOrderTriggerGadget::triggerWorkOrderAllInClose() ... ");
        return false;
    }

    return true;
}

size_t GtPlusAccumulatorWorkOrderTriggerGadget::
getDimValue(const ISMRMRD::AcquisitionHeader& acqHeader, Gadgetron::gtPlus::ISMRMRDDIM& dim)
{
    if ( dim == DIM_Encoding1 ) return acqHeader.idx.kspace_encode_step_1;
    if ( dim == DIM_Slice ) return acqHeader.idx.slice;
    if ( dim == DIM_Encoding2 ) return acqHeader.idx.kspace_encode_step_2;
    if ( dim == DIM_Contrast ) return acqHeader.idx.contrast;
    if ( dim == DIM_Phase ) return acqHeader.idx.phase;
    if ( dim == DIM_Repetition ) return acqHeader.idx.repetition;
    if ( dim == DIM_Set ) return acqHeader.idx.set;
    if ( dim == DIM_Segment ) return acqHeader.idx.segment;
    if ( dim == DIM_Average ) return acqHeader.idx.average;

    return 0;
}

void GtPlusAccumulatorWorkOrderTriggerGadget::
setDimValue(ISMRMRD::AcquisitionHeader& acqHeader, Gadgetron::gtPlus::ISMRMRDDIM& dim, size_t value)
{
    if ( dim == DIM_Encoding1 ) acqHeader.idx.kspace_encode_step_1 = value;
    if ( dim == DIM_Slice ) acqHeader.idx.slice = value;
    if ( dim == DIM_Encoding2 ) acqHeader.idx.kspace_encode_step_2 = value;
    if ( dim == DIM_Contrast ) acqHeader.idx.contrast = value;
    if ( dim == DIM_Phase ) acqHeader.idx.phase = value;
    if ( dim == DIM_Repetition ) acqHeader.idx.repetition = value;
    if ( dim == DIM_Set ) acqHeader.idx.set = value;
    if ( dim == DIM_Segment ) acqHeader.idx.segment = value;
    if ( dim == DIM_Average ) acqHeader.idx.average = value;

    return;
}

int GtPlusAccumulatorWorkOrderTriggerGadget::close(unsigned long flags)
{
    GADGET_CONDITION_MSG(true, "GtPlusAccumulatorWorkOrderTriggerGadget - close(flags) : " << flags);

    if ( BaseClass::close(flags) != GADGET_OK ) return GADGET_FAIL;

    // if ( flags!=0 && !triggered_in_close_ )
    if ( !triggered_in_close_ )
    {
        triggered_in_close_ = true;

        GADGET_CONDITION_MSG(true, "GtPlusAccumulatorWorkOrderTriggerGadget - trigger in close(flags) ... ");

        if ( needTriggerWorkOrderAllInClose() )
        {
            // never been triggered, so need to trigger with all data buffered
            if ( !triggerWorkOrderAllInClose() )
            {
                GADGET_DEBUG1("triggerWorkOrderAllInClose() failed ... \n");
                return GADGET_FAIL;
            }
        }
        else
        {
            // need to trigger the last portion of kspace
            //if ( !triggerWorkOrder(NULL, true, true) )
            //{
            //    GADGET_DEBUG1("Failed triggerWorkOrder(inClose)\n");
            //    return GADGET_FAIL;
            //}
        }
    }

    // return BaseClass::close(flags);
    return GADGET_OK;
}

GADGET_FACTORY_DECLARE(GtPlusAccumulatorWorkOrderTriggerGadget)

}
