#include "GtPlusAccumulatorGadget.h"

namespace Gadgetron
{

// --------------------------------------------------------------------

GadgetMessageImageExt::GadgetMessageImageExt() : ISMRMRD::ImageHeader()
{
    time_stamps.clear();
    pmu_time_stamps.clear();
}

GadgetMessageImageExt::~GadgetMessageImageExt() { }

void GadgetMessageImageExt::set_matrix_size(unsigned int index, ACE_UINT16 size)
{
    if (index < 3) 
    {
        matrix_size[index] = size;
    }

    if ( index == 1 )
    {
        time_stamps.clear();
        time_stamps.resize(matrix_size[1], -1);
        pmu_time_stamps.clear();
        pmu_time_stamps.resize(matrix_size[1], -1);
    }
}

void GadgetMessageImageExt::copy(GadgetMessageImageExt& aMessageImage)
{
    flags = aMessageImage.flags;

    matrix_size[0] = aMessageImage.matrix_size[0];
    matrix_size[1] = aMessageImage.matrix_size[1];
    matrix_size[2] = aMessageImage.matrix_size[2];

    channels = aMessageImage.channels;

    position[0] = aMessageImage.position[0];
    position[1] = aMessageImage.position[1];
    position[2] = aMessageImage.position[2];

    read_dir[0] = aMessageImage.read_dir[0];
    read_dir[1] = aMessageImage.read_dir[1];
    read_dir[2] = aMessageImage.read_dir[2];

    phase_dir[0] = aMessageImage.phase_dir[0];
    phase_dir[1] = aMessageImage.phase_dir[1];
    phase_dir[2] = aMessageImage.phase_dir[2];

    slice_dir[0] = aMessageImage.slice_dir[0];
    slice_dir[1] = aMessageImage.slice_dir[1];
    slice_dir[2] = aMessageImage.slice_dir[2];

    patient_table_position[0] = aMessageImage.patient_table_position[0];
    patient_table_position[1] = aMessageImage.patient_table_position[1];
    patient_table_position[2] = aMessageImage.patient_table_position[2];

    acquisition_time_stamp = aMessageImage.acquisition_time_stamp;

    physiology_time_stamp[0] = aMessageImage.physiology_time_stamp[0];
    physiology_time_stamp[1] = aMessageImage.physiology_time_stamp[1];
    physiology_time_stamp[2] = aMessageImage.physiology_time_stamp[2];

    image_data_type = aMessageImage.image_data_type;
    image_type = aMessageImage.image_type;
    image_index = aMessageImage.image_index;
    image_series_index = aMessageImage.image_series_index;

    memcpy(user_int, aMessageImage.user_int, sizeof(int32_t)*ISMRMRD_USER_INTS);
    memcpy(user_float, aMessageImage.user_float, sizeof(float)*ISMRMRD_USER_FLOATS);

    time_stamps = aMessageImage.time_stamps;
    pmu_time_stamps = aMessageImage.pmu_time_stamps;
}

void GadgetMessageImageExt::dump()
{
    std::cout << "GadgetMessageImageExt" << std::endl;
    std::cout << "----------------------------------------------------------" << std::endl;
    //dumpInfo();
    std::cout << "----------------------------------------------------------" << std::endl;
}

// --------------------------------------------------------------------

GadgetMessageImageArray::GadgetMessageImageArray() 
:   imageArray_(0)
{

}

GadgetMessageImageArray::GadgetMessageImageArray(int aSize[10])
{
    try
    {
        unsigned int ii;
        for ( ii=0; ii<10; ii++ )
        {
            matrix_size[ii] = aSize[ii];
        }

        unsigned int len = 1;
        for ( ii=3; ii<10; ii++ )
        {
            len *= matrix_size[ii];
        }

        if ( len > 0 )
        {
            imageArray_ = new GadgetMessageImageExt[len];
        }
    }
    catch(...)
    {
        std::cout << "Failed in allocate imageArray_" << std::endl;
    }
}

GadgetMessageImageArray::~GadgetMessageImageArray()
{
    if (imageArray_)
    {
        delete [] imageArray_;
    }
}

void GadgetMessageImageArray::resize(int aSize[10])
{
    try
    {
        unsigned int ii;
        for ( ii=0; ii<10; ii++ )
        {
            matrix_size[ii] = aSize[ii];
        }

        unsigned int len = 1;
        for ( ii=3; ii<10; ii++ )
        {
            len *= matrix_size[ii];
        }

        if ( imageArray_ ) 
        {
            delete [] imageArray_;
            imageArray_ = NULL;
        }

        if ( len > 0 )
        {
            imageArray_ = new GadgetMessageImageExt[len];
        }
    }
    catch(...)
    {
        std::cout << "Failed in resize GadgetMessageImageArray " << std::endl;
    }
}

void GadgetMessageImageArray::copy(GadgetMessageImageArray& imageArray)
{
    if (imageArray_) delete [] imageArray_;

    unsigned int ii;
    for ( ii=0; ii<10; ii++ )
    {
        matrix_size[ii] = imageArray.matrix_size[ii];
    }

    unsigned int len = 1;
    for ( ii=3; ii<10; ii++ )
    {
        len *= matrix_size[ii];
    }

    if ( len > 0 )
    {
        imageArray_ = new GadgetMessageImageExt[len];
    }

    for ( unsigned int i=0; i<len; i++ )
    {
        imageArray_[i] = imageArray.imageArray_[i];
    }
}

int GadgetMessageImageArray::get_offset(int slc, int e2, int con, int phs, int rep, int set, int seg)
{
    int offset = seg*matrix_size[8]*matrix_size[7]*matrix_size[6]*matrix_size[5]*matrix_size[4]*matrix_size[3]
                    + set*matrix_size[7]*matrix_size[6]*matrix_size[5]*matrix_size[4]*matrix_size[3]
                    + rep*matrix_size[6]*matrix_size[5]*matrix_size[4]*matrix_size[3]
                    + phs*matrix_size[5]*matrix_size[4]*matrix_size[3]
                    + con*matrix_size[4]*matrix_size[3]
                    + e2*matrix_size[3]
                    + slc;
    return offset;
}

void GadgetMessageImageArray::extractMessageImageArrayForSLC(int slc, GadgetMessageImageArray& imageArray)
{
    if ( slc >= matrix_size[3] )
    {
        std::cout << "extractMessageImageArrayForSLC error - slc >= matrix_size[3] " << std::endl;
        return;
    }

    int aSize[10];

    unsigned int ii;
    for ( ii=0; ii<10; ii++ )
    {
        aSize[ii] = matrix_size[ii];
    }

    aSize[3] = 1;

    imageArray.resize(aSize);

    int e2, con, phs, rep, set, seg;

    int E2 = matrix_size[4];
    int CON = matrix_size[5];
    int PHS = matrix_size[6];
    int REP = matrix_size[7];
    int SET = matrix_size[8];
    int SEG = matrix_size[9];

    for ( seg=0; seg<SEG; seg++ )
    {
        for ( set=0; set<SET; set++ )
        {
            for ( rep=0; rep<REP; rep++ )
            {
                for ( phs=0; phs<PHS; phs++ )
                {
                    for ( con=0; con<CON; con++ )
                    {
                        for ( e2=0; e2<E2; e2++ )
                        {
                            int offset = this->get_offset(slc, e2, con, phs, rep, set, seg);
                            int offsetSLC = imageArray.get_offset(0, e2, con, phs, rep, set, seg);

                            imageArray.imageArray_[offsetSLC] = imageArray_[offset];
                        }
                    }
                }
            }
        }
    }
}

void GadgetMessageImageArray::extractMessageImageArrayForREP(int rep, GadgetMessageImageArray& imageArray)
{
    if ( rep >= matrix_size[7] )
    {
        std::cout << "extractMessageImageArrayForSLC error - rep >= matrix_size[7] " << std::endl;
        return;
    }

    int aSize[10];

    unsigned int ii;
    for ( ii=0; ii<10; ii++ )
    {
        aSize[ii] = matrix_size[ii];
    }

    aSize[7] = 1;

    imageArray.resize(aSize);

    int e2, con, phs, slc, set, seg;

    int SLC = matrix_size[3];
    int E2 = matrix_size[4];
    int CON = matrix_size[5];
    int PHS = matrix_size[6];
    int SET = matrix_size[8];
    int SEG = matrix_size[9];

    for ( seg=0; seg<SEG; seg++ )
    {
        for ( set=0; set<SET; set++ )
        {
            for ( slc=0; slc<SLC; slc++ )
            {
                for ( phs=0; phs<PHS; phs++ )
                {
                    for ( con=0; con<CON; con++ )
                    {
                        for ( e2=0; e2<E2; e2++ )
                        {
                            int offset = this->get_offset(slc, e2, con, phs, rep, set, seg);
                            int offsetREP = imageArray.get_offset(slc, e2, con, phs, 0, set, seg);

                            imageArray.imageArray_[offsetREP] = imageArray_[offset];
                        }
                    }
                }
            }
        }
    }
}

void GadgetMessageImageArray::dump()
{
    unsigned int ii;
    std::cout << "GadgetMessageImageArray" << std::endl;
    std::cout << "==========================================================" << std::endl;
    std::cout << "matrix_size           : ";
    for ( ii=0; ii<10; ii++ )
    {
        std::cout << matrix_size[ii] << " ";
    }
    std::cout << std::endl;
    std::cout << "----------------------------------------------------------" << std::endl;
    if ( imageArray_ )
    {
        int slc, e2, con, phs, rep, set, seg;
        for ( seg=0; seg<matrix_size[9]; seg++ )
        {
            for ( set=0; set<matrix_size[8]; set++ )
            {
                for ( rep=0; rep<matrix_size[7]; rep++ )
                {
                    for ( phs=0; phs<matrix_size[6]; phs++ )
                    {
                        for ( con=0; con<matrix_size[5]; con++ )
                        {
                            for ( e2=0; e2<matrix_size[4]; e2++ )
                            {
                                for ( slc=0; slc<matrix_size[3]; slc++ )
                                {
                                    int offset = get_offset(slc, e2, con, phs, rep, set, seg);
                                    std::cout << "[Slice E2 Contrast Phase Rep Set Seg] = [" 
                                                << " " << slc 
                                                << " " << e2 
                                                << " " << con 
                                                << " " << phs 
                                                << " " << rep 
                                                << " " << set 
                                                << " " << seg << "]" << std::endl;

                                    imageArray_[offset].dump();
                                }
                            }
                        }
                    }
                }
            }
        }

    }
    std::cout << "==========================================================" << std::endl;
}

// --------------------------------------------------------------------

KSpaceBuffer::KSpaceBuffer() 
{

}

KSpaceBuffer::~KSpaceBuffer()
{

}

// --------------------------------------------------------------------

GtPlusAccumulatorGadget::GtPlusAccumulatorGadget()
    : messageImage_(0)
    , kspaceBuffer_(0)
    , image_counter_(0)
    , image_series_(0)
    , triggered_(false)
{

}

GtPlusAccumulatorGadget::~GtPlusAccumulatorGadget()
{
    if (messageImage_) delete messageImage_;
    if (kspaceBuffer_) delete kspaceBuffer_;
}

// extract necessary configuration information from the xml
int GtPlusAccumulatorGadget::process_config(ACE_Message_Block* mb)
{

    // allocate the kspace buffer
    if ( kspaceBuffer_ == NULL )
    {
        if (!(kspaceBuffer_ = new KSpaceBuffer)) 
        {
            GADGET_DEBUG1("Failed create buffer\n");
            return GADGET_FAIL;
        }
    }

    // image series
    image_series_ = this->get_int_value("image_series");

    // pass the xml file
    boost::shared_ptr<ISMRMRD::ismrmrdHeader> cfg = parseIsmrmrdXMLHeader(std::string(mb->rd_ptr()));

    // seq object
    ISMRMRD::ismrmrdHeader::encoding_sequence e_seq = cfg->encoding();
    if (e_seq.size() != 1)
    {
        GADGET_DEBUG2("Number of encoding spaces: %d\n", e_seq.size());
        GADGET_DEBUG1("This simple GtPlusAccumulatorGadget only supports one encoding space\n");
        return GADGET_FAIL;
    }

    ISMRMRD::encodingSpaceType e_space = (*e_seq.begin()).encodedSpace();
    ISMRMRD::encodingSpaceType r_space = (*e_seq.begin()).reconSpace();
    ISMRMRD::encodingLimitsType e_limits = (*e_seq.begin()).encodingLimits();

    GADGET_MSG("Matrix size: " << e_space.matrixSize().x() << " " << e_space.matrixSize().y() << " " << e_space.matrixSize().z());
    GADGET_MSG("Recon size: " << r_space.matrixSize().x() << " " << r_space.matrixSize().y() << " " << r_space.matrixSize().z());

    meas_max_ro_ = e_space.matrixSize().x()/2;

    field_of_view_[0] = r_space.fieldOfView_mm().x();
    field_of_view_[1] = r_space.fieldOfView_mm().y();
    field_of_view_[2] = r_space.fieldOfView_mm().z();
    GADGET_MSG("field_of_view_ is " << field_of_view_[0] << " " << field_of_view_[1] << " " << field_of_view_[2]);

    int newE1_ = field_of_view_[1]/(field_of_view_[0]/meas_max_ro_);

    if (e_limits.kspace_encoding_step_1().present()) 
    {
        meas_max_idx_.kspace_encode_step_1 = e_limits.kspace_encoding_step_1().get().maximum();
    }
    else
    {
        meas_max_idx_.kspace_encode_step_1 = 0;
        std::cout << "Setting number of kspace_encode_step_1 to 0" << std::endl;
        return GADGET_FAIL;
    }

    kspaceBuffer_->kSpaceCentreEncode1_ = e_limits.kspace_encoding_step_1().get().center();
    GADGET_MSG("kSpaceCentreEncode1_ is " << kspaceBuffer_->kSpaceCentreEncode1_);

    kspaceBuffer_->kSpaceCentreEncode2_ = e_limits.kspace_encoding_step_2().get().center();
    GADGET_MSG("kSpaceCentreEncode2_ is " << kspaceBuffer_->kSpaceCentreEncode2_);

    kspaceBuffer_->kSpaceMaxEncode1_ = e_limits.kspace_encoding_step_1().get().maximum()+1;
    GADGET_MSG("kSpaceMaxEncode1_ is " << kspaceBuffer_->kSpaceMaxEncode1_);

    kspaceBuffer_->kSpaceMaxEncode2_ = e_limits.kspace_encoding_step_2().get().maximum()+1;
    GADGET_MSG("kSpaceMaxEncode2_ is " << kspaceBuffer_->kSpaceMaxEncode2_);

    if (e_limits.set().present())
    {
        meas_max_idx_.set = e_limits.set().get().maximum() - 1;
        if ( meas_max_idx_.set < 0 ) meas_max_idx_.set = 0;
    }
    else
    {
        meas_max_idx_.set = 0;
    }

    if (e_limits.phase().present())
    {
        meas_max_idx_.phase = e_limits.phase().get().maximum()-1;
        if ( meas_max_idx_.phase < 0 ) meas_max_idx_.phase = 0;
    }
    else
    {
        meas_max_idx_.phase = 0;
    }

    if (e_limits.kspace_encoding_step_2().present())
    {
        meas_max_idx_.kspace_encode_step_2 = e_limits.kspace_encoding_step_2().get().maximum();
    }
    else
    {
        meas_max_idx_.kspace_encode_step_2 = 0;
    }

    if (e_limits.contrast().present())
    {
        meas_max_idx_.contrast = e_limits.contrast().get().maximum()-1;
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

    if (e_limits.segment().present())
    {
        // meas_max_idx_.segment = e_limits.segment().get().maximum()-1;
        meas_max_idx_.segment = 0;
    }
    else
    {
        meas_max_idx_.segment = 0;
    }

    // find out the PAT mode
    ISMRMRD::ismrmrdHeader::parallelImaging_optional p_imaging_type = cfg->parallelImaging();
    ISMRMRD::parallelImagingType p_imaging = *p_imaging_type;

    kspaceBuffer_->AccelFactE1_ = (unsigned int)(p_imaging.accelerationFactor().kspace_encoding_step_1());
    kspaceBuffer_->AccelFactE2_ = (unsigned int)(p_imaging.accelerationFactor().kspace_encoding_step_2());
    GADGET_MSG("AccelFactE1 is " << kspaceBuffer_->AccelFactE1_);
    GADGET_MSG("AccelFactE2 is " << kspaceBuffer_->AccelFactE2_);

    ISMRMRD::calibrationModeType calib = *(p_imaging.calibrationMode());
    kspaceBuffer_->CalibMode_ = calib;

    // find out the calibration mode
    if ( kspaceBuffer_->CalibMode_ == ISMRMRD::calibrationModeType::separate )
    {
        GADGET_MSG("Calibration mode is separate");
    }

    if ( kspaceBuffer_->CalibMode_ == ISMRMRD::calibrationModeType::embedded )
    {
        GADGET_MSG("Calibration mode is embedded");
    }

    if ( kspaceBuffer_->CalibMode_ == ISMRMRD::calibrationModeType::interleaved )
    {
        GADGET_MSG("Calibration mode is interleaved");

        if ( p_imaging.interleavingDimension().present() )
        {
            kspaceBuffer_->InterleaveDim_ = *(p_imaging.interleavingDimension());
            GADGET_MSG("InterleaveDim is " << kspaceBuffer_->InterleaveDim_);
        }
    }

    return GADGET_OK;
}

int GtPlusAccumulatorGadget::process(GadgetContainerMessage<ISMRMRD::AcquisitionHeader>* m1, GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2)
{
    // logic to control whether to store kspace and ref data
    bool bIsKSpace, bIsRef, bIsNoise, bIsPhaseCorr, bIsReflect, bIsOther;
    if ( !checkStatus(m1->getObjectPtr()->flags, m1->getObjectPtr()->number_of_samples, bIsKSpace, bIsRef, bIsNoise, bIsPhaseCorr, bIsReflect, bIsOther) )
    {
        GADGET_DEBUG1("Failed check readout status\n");
        return GADGET_FAIL;
    }

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
        ISMRMRD::AcquisitionHeader* pMDH = m1->getObjectPtr();
        hoNDArray< ValueType >* pRefLine = m2->getObjectPtr();

        ReadOutBuffer item;
        item.acqHead_ = *pMDH;
        item.data_ = *pRefLine;
        item.isReflect_ = bIsReflect;
        refBuffer_.push_back(item);
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

        ReadOutBuffer item;
        item.acqHead_ = *pMDH;
        item.data_ = *pRefLine;
        item.isReflect_ = bIsReflect;
        otherBuffer_.push_back(item);
    }

    m1->release();
    return GADGET_OK;
}

bool GtPlusAccumulatorGadget::checkStatus(uint64_t flag, int samples, bool& bIsKSpace, bool& bIsRef, bool& bIsNoise, bool& bIsPhaseCorr, bool& bIsReflect, bool& bIsOther)
{
    bIsNoise = ISMRMRD::FlagBit(ISMRMRD::ACQ_IS_NOISE_MEASUREMENT).isSet(flag);
    bool is_ref = ISMRMRD::FlagBit(ISMRMRD::ACQ_IS_PARALLEL_CALIBRATION).isSet(flag);
    bool is_ref_kspace = ISMRMRD::FlagBit(ISMRMRD::ACQ_IS_PARALLEL_CALIBRATION_AND_IMAGING).isSet(flag);
    bIsReflect = ISMRMRD::FlagBit(ISMRMRD::ACQ_IS_REVERSE).isSet(flag);
    bIsPhaseCorr = ISMRMRD::FlagBit(ISMRMRD::ACQ_IS_PHASECORR_DATA).isSet(flag);

    bIsKSpace = false;
    bIsRef = false;
    bIsOther = false;

    if ( bIsNoise || bIsPhaseCorr )
    {
        return true;
    }

    // in interleaved mode, only store the image data
    if ( kspaceBuffer_->CalibMode_==ISMRMRD::calibrationModeType::interleaved )
    {
        bIsKSpace = true;
        bIsRef = false;
    }

    // in embedded, kspace stores only the undersampled lines
    // ref stores all lines used for references
    if ( kspaceBuffer_->CalibMode_==ISMRMRD::calibrationModeType::embedded )
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
    if ( kspaceBuffer_->CalibMode_==ISMRMRD::calibrationModeType::separate 
    || kspaceBuffer_->CalibMode_==ISMRMRD::calibrationModeType::external )
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

bool GtPlusAccumulatorGadget::storeImageData(GadgetContainerMessage<ISMRMRD::AcquisitionHeader>* m1, GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2, bool isReflect)
{
    try
    {
        unsigned int ii;
        int samples =  m1->getObjectPtr()->number_of_samples;
        ISMRMRD::EncodingCounters idx = m1->getObjectPtr()->idx;

        if ( kspaceBuffer_->buffer_.get_number_of_elements() <= 0 )
        {
            meas_max_channel_ = m1->getObjectPtr()->active_channels;

            int E1 = 2*kspaceBuffer_->kSpaceCentreEncode1_;
            int E2 = 2*kspaceBuffer_->kSpaceCentreEncode2_;

            // find the loop counter boundary and allocate the buffer
            GADGET_MSG("[RO E1 Cha Slice E2 Con Phase Rep Set Seg] = [" 
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

            unsigned int N = dimensions_.size();
            for ( ii=0; ii<N; ii++ )
            {
                GADGET_MSG("dimensions_[" << ii << "] = " << dimensions_[ii]);
            }

            // allocate data buffer
            try
            {
                kspaceBuffer_->buffer_.create(&dimensions_);

                std::vector<unsigned int> reflect_dimensions_(dimensions_);
                reflect_dimensions_[0] = 1;
                reflect_dimensions_[2] = 1;
                kspaceBuffer_->reflect_.create(&reflect_dimensions_);
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

            if (!(messageImage_ = new GadgetMessageImageArray(matrix_size))) 
            {
                GADGET_DEBUG1("Failed create buffer\n");
                return false;
            }
        }

        std::complex<float>* b = kspaceBuffer_->buffer_.begin();
        std::complex<float>* d = m2->getObjectPtr()->get_data_ptr();
        if (samples != static_cast<int>(dimensions_[0])) 
        {
            GADGET_DEBUG1("Wrong number of samples received\n");
            return false;
        }

        //Copy the data for all the channels
        std::vector<unsigned int> pos(10);
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
            int offsetBuffer = kspaceBuffer_->buffer_.calculate_offset(pos);

            memcpy(b+offsetBuffer, d+c*samples, sizeof(std::complex<float>)*samples);

            pos[2] = 0;
            offsetBuffer = kspaceBuffer_->reflect_.calculate_offset(pos);
            kspaceBuffer_->reflect_.at(offsetBuffer) = isReflect;
        }

        if ( !fillImageInfo(m1, messageImage_, m1->getObjectPtr()->idx) )
        {
            GADGET_DEBUG1("Failed in fillImageInfo(m1, messageImage_, m1->getObjectPtr()->idx)\n");
            return false;
        }
    }
    catch(...)
    {
        GADGET_DEBUG1("Errors in GtPlusAccumulatorGadget::storeImageData(...) ... \n");
        return false;
    }

    return true;
}

bool GtPlusAccumulatorGadget::
fillBuffer(ReadOutBufferType& readOutBuffer, BufferType& buf, ReflectBufferType& reflectBuf)
{
    try
    {
        // find the maximal dimension of all buffered ICE readout
        unsigned int numOfReadOuts = readOutBuffer.size();
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

        unsigned int a;
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

        GADGET_MSG("[RO E1 Cha Slice E2 Contrast Phase Rep Set Seg] = [" 
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
        std::vector<unsigned int> dims(10);
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

            std::vector<unsigned int> reflect_dims(dims);
            reflect_dims[0] = 1;
            reflect_dims[2] = 1;
            reflectBuf.create(&reflect_dims);
        }
        catch(...)
        {
            GADGET_DEBUG1("Failed create buffer\n");
            return false;
        }

        std::complex<float>* b = buf.begin();

        // copy the data
        int c;
        std::vector<unsigned int> pos(10);

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
                int offsetBuffer = buf.calculate_offset(pos);

                memcpy(b+offsetBuffer, d+c*readOutBuffer[a].acqHead_.number_of_samples, sizeof(std::complex<float>)*readOutBuffer[a].acqHead_.number_of_samples);

                pos[2] = 0;
                offsetBuffer = reflectBuf.calculate_offset(pos);
                reflectBuf.at(offsetBuffer) = readOutBuffer[a].isReflect_;
            }
        }
    }
    catch(...)
    {
        GADGET_DEBUG1("Errors in GtPlusAccumulatorGadget::fillBuffer(...) ... \n");
        return false;
    }

    return true;
}

bool GtPlusAccumulatorGadget::fillImageInfo(GadgetContainerMessage<ISMRMRD::AcquisitionHeader>* m1, GadgetMessageImageArray* messageImage, const ISMRMRD::EncodingCounters& idx)
{
    try
    {
        // fill the message info
        int offset = messageImage->get_offset(idx.slice, idx.kspace_encode_step_2, idx.contrast, idx.phase, idx.repetition, idx.set, idx.segment);

        // if it is the first acq in a slice, fill in all information
        bool is_first_acq_in_slice = ISMRMRD::FlagBit(ISMRMRD::ACQ_FIRST_IN_SLICE).isSet(m1->getObjectPtr()->flags);

        if ( is_first_acq_in_slice )
        {
            messageImage->imageArray_[offset].version = m1->getObjectPtr()->version;
            messageImage->imageArray_[offset].flags = m1->getObjectPtr()->flags;
            messageImage->imageArray_[offset].measurement_uid = m1->getObjectPtr()->measurement_uid;

            //messageImage->imageArray_[offset].matrix_size[0] = dimensions_[0];
            //messageImage->imageArray_[offset].matrix_size[1] = dimensions_[1];
            //messageImage->imageArray_[offset].matrix_size[2] = dimensions_[2];

            messageImage->imageArray_[offset].set_matrix_size(0, dimensions_[0]);
            messageImage->imageArray_[offset].set_matrix_size(1, dimensions_[1]);
            messageImage->imageArray_[offset].set_matrix_size(2, dimensions_[2]);

            messageImage->imageArray_[offset].field_of_view[0] = field_of_view_[0];
            messageImage->imageArray_[offset].field_of_view[1] = field_of_view_[1];
            messageImage->imageArray_[offset].field_of_view[2] = field_of_view_[2];

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
        GADGET_DEBUG1("Errors in GtPlusAccumulatorGadget::fillImageInfo(...) ... \n");
        return false;
    }

    return true;
}

int GtPlusAccumulatorGadget::close(unsigned long flags)
{
    if ( !triggered_ )
    {
        triggered_ = true;

        GADGET_MSG("GtPlusAccumulatorGadget - trigger next gadget ... ");

        GadgetContainerMessage<GadgetMessageImageArray>* cm1 = 
            new GadgetContainerMessage<GadgetMessageImageArray>();

        GadgetContainerMessage< KSpaceBuffer >* cm2 = 
            new GadgetContainerMessage< KSpaceBuffer >();

        cm1->cont(cm2);

        // copy the image content
        cm2->getObjectPtr()->buffer_ = kspaceBuffer_->buffer_;
        cm2->getObjectPtr()->reflect_ = kspaceBuffer_->reflect_;

        // copy the message image array
        cm1->getObjectPtr()->copy(*messageImage_);

        if (!refBuffer_.empty())
        {
            GADGET_MSG("GtPlusAccumulatorGadget - ref signal found : " << refBuffer_.size());

            if ( !fillBuffer(refBuffer_, kspaceBuffer_->ref_, kspaceBuffer_->refReflect_) )
            {
                GADGET_DEBUG1("fillBuffer(refBuffer_) failed ... \n");
                cm1->release();
                return GADGET_FAIL;
            }

            cm2->getObjectPtr()->ref_ = kspaceBuffer_->ref_;
            cm2->getObjectPtr()->refReflect_ = kspaceBuffer_->refReflect_;
        }

        if (!phaseCorrBuffer_.empty())
        {
            GADGET_MSG("GtPlusAccumulatorGadget - phase correction signal found : " << phaseCorrBuffer_.size());

            if ( !fillBuffer(phaseCorrBuffer_, kspaceBuffer_->phaseCorr_, kspaceBuffer_->phaseCorrReflect_) )
            {
                GADGET_DEBUG1("fillBuffer(phaseCorrBuffer_) failed ... \n");
                cm1->release();
                return GADGET_FAIL;
            }

            cm2->getObjectPtr()->phaseCorr_ = kspaceBuffer_->phaseCorr_;
            cm2->getObjectPtr()->phaseCorrReflect_ = kspaceBuffer_->phaseCorrReflect_;
        }

        if (!noiseBuffer_.empty())
        {
            GADGET_MSG("GtPlusAccumulatorGadget - noise signal found : " << noiseBuffer_.size());

            ReflectBufferType tmpBuf;
            if ( !fillBuffer(noiseBuffer_, kspaceBuffer_->noise_, tmpBuf) )
            {
                GADGET_DEBUG1("fillBuffer(noiseBuffer_) failed ... \n");
                cm1->release();
                return GADGET_FAIL;
            }

            cm2->getObjectPtr()->noise_ = kspaceBuffer_->noise_;
        }

        if (!otherBuffer_.empty())
        {
            GADGET_MSG("GtPlusAccumulatorGadget - other signal found : " << otherBuffer_.size());

            ReflectBufferType tmpBuf;
            if ( !fillBuffer(otherBuffer_, kspaceBuffer_->other_, tmpBuf) )
            {
                GADGET_DEBUG1("fillBuffer(otherBuffer_) failed ... \n");
                cm1->release();
                return GADGET_FAIL;
            }

            cm2->getObjectPtr()->other_ = kspaceBuffer_->other_;
        }

        // send to next gadget
        if (this->next()->putq(cm1) < 0) 
        {
            return GADGET_FAIL;
        }
    }

    return BaseClass::close(flags);
}

GADGET_FACTORY_DECLARE(GtPlusAccumulatorGadget)

}
