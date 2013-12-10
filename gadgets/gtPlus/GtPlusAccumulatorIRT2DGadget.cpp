
#include "GtPlusAccumulatorIRT2DGadget.h"
#include "GadgetIsmrmrdReadWrite.h"

namespace Gadgetron
{

GtPlusAccumulatorIRT2DGadget::GtPlusAccumulatorIRT2DGadget() : prev_rep_(-1), cur_rep_(-1), num_scan_buffered_(0)
{

}

GtPlusAccumulatorIRT2DGadget::~GtPlusAccumulatorIRT2DGadget()
{

}

int GtPlusAccumulatorIRT2DGadget::process_config(ACE_Message_Block* mb)
{
    return BaseClass::process_config(mb);
}

bool GtPlusAccumulatorIRT2DGadget::
copyBufferForREP(BufferType& buf, int rep, BufferType& bufREP)
{
    try
    {
        boost::shared_ptr< std::vector<unsigned int> > dims = buf.get_dimensions();

        boost::shared_ptr< std::vector<unsigned int> > dimsREP = dims;
        (*dimsREP)[7] = 1;

        try
        {
            bufREP.create(dimsREP);
        }
        catch(...)
        {
            GADGET_DEBUG1("Failed create buffer for REP \n");
            return false;
        }

        // copy the memory over
        int RO = (*dims)[0];
        int E1 = (*dims)[1];
        int CHA = (*dims)[2];
        int SLC = (*dims)[3];
        int E2 = (*dims)[4];
        int CON = (*dims)[5];
        int PHS = (*dims)[6];
        int REP = (*dims)[7];
        int SET = (*dims)[8];
        int SEG = (*dims)[9];

        int e2, con, phs, slc, set, seg;

        std::vector<unsigned int> pos(10);

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
                                pos[0] = 0;
                                pos[1] = 0;
                                pos[2] = 0;
                                pos[3] = slc;
                                pos[4] = e2;
                                pos[5] = con;
                                pos[6] = phs;
                                pos[7] = rep;
                                pos[8] = set;
                                pos[9] = seg;
                                int offsetBuffer = buf.calculate_offset(pos);

                                // buffer slc
                                pos[0] = 0;
                                pos[1] = 0;
                                pos[2] = 0;
                                pos[3] = slc;
                                pos[4] = e2;
                                pos[5] = con;
                                pos[6] = phs;
                                pos[7] = 0;
                                pos[8] = set;
                                pos[9] = seg;
                                int offsetBufferREP = bufREP.calculate_offset(pos);

                                // copy the image content
                                memcpy(bufREP.begin()+offsetBufferREP, buf.begin()+offsetBuffer, sizeof(std::complex<float>)*RO*E1*CHA);
                            }
                        }
                    }
                }
            }
        }
    }
    catch(...)
    {
        GADGET_DEBUG1("Errors in GtPlusAccumulatorIRT2DGadget::copyBufferForREP(...) ... \n");
        return false;
    }

    return true;
}

bool GtPlusAccumulatorIRT2DGadget::
copyReflectBufferForREP(ReflectBufferType& buf, int rep, ReflectBufferType& bufREP)
{
    try
    {
        boost::shared_ptr< std::vector<unsigned int> > dims = buf.get_dimensions();

        boost::shared_ptr< std::vector<unsigned int> > dimsREP = dims;
        (*dimsREP)[7] = 1;

        try
        {
            bufREP.create(dimsREP);
        }
        catch(...)
        {
            GADGET_DEBUG1("Failed create buffer for REP \n");
            return false;
        }

        // copy the memory over
        int RO = (*dims)[0];
        int E1 = (*dims)[1];
        int CHA = (*dims)[2];
        int SLC = (*dims)[3];
        int E2 = (*dims)[4];
        int CON = (*dims)[5];
        int PHS = (*dims)[6];
        int REP = (*dims)[7];
        int SET = (*dims)[8];
        int SEG = (*dims)[9];

        int e2, con, phs, slc, set, seg;

        std::vector<unsigned int> pos(10);

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
                                pos[0] = 0;
                                pos[1] = 0;
                                pos[2] = 0;
                                pos[3] = slc;
                                pos[4] = e2;
                                pos[5] = con;
                                pos[6] = phs;
                                pos[7] = rep;
                                pos[8] = set;
                                pos[9] = seg;
                                int offsetBuffer = buf.calculate_offset(pos);

                                // buffer slc
                                pos[0] = 0;
                                pos[1] = 0;
                                pos[2] = 0;
                                pos[3] = slc;
                                pos[4] = e2;
                                pos[5] = con;
                                pos[6] = phs;
                                pos[7] = 0;
                                pos[8] = set;
                                pos[9] = seg;
                                int offsetBufferREP = bufREP.calculate_offset(pos);

                                // copy the image content
                                memcpy(bufREP.begin()+offsetBufferREP, buf.begin()+offsetBuffer, sizeof(unsigned short)*RO*E1*CHA);
                            }
                        }
                    }
                }
            }
        }
    }
    catch(...)
    {
        GADGET_DEBUG1("Errors in GtPlusAccumulatorSLCGadget::copyReflectBufferForSLC(...) ... \n");
        return false;
    }

    return true;
}

bool GtPlusAccumulatorIRT2DGadget::triggerREP(int rep)
{
    try
    {
        GadgetContainerMessage<GadgetMessageImageArray>* cm1 = 
            new GadgetContainerMessage<GadgetMessageImageArray>();

        GadgetContainerMessage< KSpaceBuffer >* cm2 = 
            new GadgetContainerMessage< KSpaceBuffer >();

        cm1->cont(cm2);

        // copy the kspace data for this REP
        if ( !copyBufferForREP(kspaceBuffer_->buffer_, 0, cm2->getObjectPtr()->buffer_) ) 
        {
            GADGET_DEBUG1("Unable to copyBufferForREP\n");
            cm1->release();
            return false;
        }

        if ( !copyReflectBufferForREP(kspaceBuffer_->reflect_, 0, cm2->getObjectPtr()->reflect_) ) 
        {
            GADGET_DEBUG1("Unable to copyReflectBufferForREP\n");
            cm1->release();
            return false;
        }

        // fill buffer with zeros, ready for next REP
        kspaceBuffer_->buffer_.fill(0);
        kspaceBuffer_->reflect_.fill(0);

        // copy the message image array for this REP
        GadgetMessageImageArray aMessageArray;
        messageImage_->extractMessageImageArrayForREP(0, aMessageArray);
        cm1->getObjectPtr()->copy(aMessageArray);

        if (!refBuffer_.empty())
        {
            GADGET_MSG("GtPlusAccumulatorIRT2DGadget - ref signal found : " << refBuffer_.size());

            BufferType refCurr;
            ReflectBufferType refReflectCurr;
            if ( !fillBuffer(refBuffer_, refCurr, refReflectCurr) )
            {
                GADGET_DEBUG1("fillBuffer(refBuffer_, refCurr, refReflectCurr) failed ... \n");
                cm1->release();
                return GADGET_FAIL;
            }

            if ( !copyBufferForREP(refCurr, rep, cm2->getObjectPtr()->ref_) ) 
            {
                GADGET_DEBUG1("Unable to copyBufferForREP(refCurr, rep, cm2->getObjectPtr()->ref_)\n");
                cm1->release();
                return false;
            }

            if ( !copyReflectBufferForREP(refReflectCurr, rep, cm2->getObjectPtr()->refReflect_) ) 
            {
                GADGET_DEBUG1("Unable to copyReflectBufferForREP(refReflectCurr, rep, cm2->getObjectPtr()->refReflect_)\n");
                cm1->release();
                return false;
            }
        }

        if (!phaseCorrBuffer_.empty())
        {
            GADGET_MSG("GtPlusAccumulatorIRT2DGadget - phase correction signal found : " << phaseCorrBuffer_.size());

            BufferType phsCorrCurr;
            ReflectBufferType phsCorrReflectCurr;
            if ( !fillBuffer(phaseCorrBuffer_, phsCorrCurr, phsCorrReflectCurr) )
            {
                GADGET_DEBUG1("fillBuffer(phaseCorrBuffer_, phsCorrCurr, phsCorrReflectCurr) failed ... \n");
                cm1->release();
                return GADGET_FAIL;
            }

            if ( !copyBufferForREP(phsCorrCurr, rep, cm2->getObjectPtr()->phaseCorr_) ) 
            {
                GADGET_DEBUG1("Unable to copyBufferForREP(phsCorrCurr, rep, cm2->getObjectPtr()->phaseCorr_)\n");
                cm1->release();
                return false;
            }

            if ( !copyReflectBufferForREP(phsCorrReflectCurr, rep, cm2->getObjectPtr()->phaseCorrReflect_) ) 
            {
                GADGET_DEBUG1("Unable to copyReflectBufferForREP(phsCorrReflectCurr, rep, cm2->getObjectPtr()->phaseCorrReflect_)\n");
                cm1->release();
                return false;
            }
        }

        if (!noiseBuffer_.empty())
        {
            GADGET_MSG("GtPlusAccumulatorIRT2DGadget - noise signal found : " << noiseBuffer_.size());

            BufferType noiseCurr;
            ReflectBufferType tmpBuf;
            if ( !fillBuffer(noiseBuffer_, noiseCurr, tmpBuf) )
            {
                GADGET_DEBUG1("fillBuffer(noiseBuffer_, noiseCurr, tmpBuf) failed ... \n");
                cm1->release();
                return false;
            }

            if ( !copyBufferForREP(noiseCurr, rep, cm2->getObjectPtr()->noise_) ) 
            {
                GADGET_DEBUG1("Unable to copyBufferForREP(noiseCurr, rep, cm2->getObjectPtr()->noise_)\n");
                cm1->release();
                return false;
            }
        }

        if (!otherBuffer_.empty())
        {
            GADGET_MSG("GtPlusAccumulatorIRT2DGadget - other signal found : " << otherBuffer_.size());

            BufferType otherCurr;
            ReflectBufferType tmpBuf;
            if ( !fillBuffer(otherBuffer_, otherCurr, tmpBuf) )
            {
                GADGET_DEBUG1("fillBuffer(otherBuffer_, otherCurr, tmpBuf) failed ... \n");
                cm1->release();
                return false;
            }

            if ( !copyBufferForREP(otherCurr, rep, cm2->getObjectPtr()->other_) ) 
            {
                GADGET_DEBUG1("Unable to copyBufferForSLC(otherCurr, rep, cm2->getObjectPtr()->other_)\n");
                cm1->release();
                return false;
            }
        }

        // send to next gadget
        if (this->next()->putq(cm1) < 0) 
        {
            return false;
        }
    }
    catch(...)
    {
        GADGET_DEBUG1("Errors in GtPlusAccumulatorIRT2DGadget::triggerREP(rep) ... \n");
        return false;
    }

    return true;
}

bool GtPlusAccumulatorIRT2DGadget::
storeImageData(GadgetContainerMessage<ISMRMRD::AcquisitionHeader>* m1, GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2, bool isReflect)
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
                               << " " << 1 
                               << " " << meas_max_idx_.set+1 
                               << " " << meas_max_idx_.segment+1 << "]");

            dimensions_.push_back(meas_max_ro_);
            dimensions_.push_back(E1);
            dimensions_.push_back(meas_max_channel_);
            dimensions_.push_back(meas_max_idx_.slice+1);
            dimensions_.push_back(E2);
            dimensions_.push_back(meas_max_idx_.contrast+1);
            dimensions_.push_back(meas_max_idx_.phase+1);
            dimensions_.push_back(1);
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
            pos[7] = 0;
            pos[8] = idx.set;
            pos[9] = idx.segment;
            int offsetBuffer = kspaceBuffer_->buffer_.calculate_offset(pos);

            memcpy(b+offsetBuffer, d+c*samples, sizeof(std::complex<float>)*samples);

            pos[2] = 0;
            offsetBuffer = kspaceBuffer_->reflect_.calculate_offset(pos);
            kspaceBuffer_->reflect_.at(offsetBuffer) = isReflect;
        }

        idx.repetition = 0;
        if ( !fillImageInfo(m1, messageImage_, idx) )
        {
            GADGET_DEBUG1("Failed in fillImageInfo(m1, messageImage_, idx)\n");
            return false;
        }
    }
    catch(...)
    {
        GADGET_DEBUG1("Errors in GtPlusAccumulatorIRT2DGadget::storeImageData(...) ... \n");
        return false;
    }

    return true;
}

int GtPlusAccumulatorIRT2DGadget::process(GadgetContainerMessage<ISMRMRD::AcquisitionHeader>* m1, 
        GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2)
{
    // check whether a new REP starts
    bool isLastScanInSlice = ISMRMRD::FlagBit(ISMRMRD::ACQ_LAST_IN_SLICE).isSet(m1->getObjectPtr()->flags);

    bool isNewRep = false;
    cur_rep_ = m1->getObjectPtr()->idx.repetition;

    if ( prev_rep_==-1 )
    {
        prev_rep_ = cur_rep_;
    }
    else
    {
        if ( cur_rep_!=prev_rep_ )
        {
            isNewRep = true;
        }
    }

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

        num_scan_buffered_++;
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

        if ( isNewRep )
        {
            refBuffer_.clear();
        }

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

        if ( isNewRep )
        {
            phaseCorrBuffer_.clear();
        }

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

        if ( isNewRep )
        {
            noiseBuffer_.clear();
        }

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

        if ( isNewRep )
        {
            otherBuffer_.clear();
        }

        otherBuffer_.push_back(item);
    }

   // if a new rep comes, it indicates the previous one is complete and can be sent out
    if ( isLastScanInSlice )
    {
        // GADGET_MSG("Repetition " << prev_rep_ << " is complete ... ");

        if ( !triggerREP(prev_rep_) ) 
        {
            GADGET_DEBUG1("Unable to trigger this rep ... \n");
            return GADGET_FAIL;
        }

        prev_rep_ = cur_rep_;

        GADGET_ERROR_MSG("GtPlusAccumulatorIRT2DGadget - trigger next gadget for REP " << prev_rep_ << " - scan buffered - " << num_scan_buffered_ << " ... ");
        num_scan_buffered_ = 0;
    }

    m1->release();
    return GADGET_OK;
}

GADGET_FACTORY_DECLARE(GtPlusAccumulatorIRT2DGadget)
}
