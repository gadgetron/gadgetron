#include "GtPlusAccumulatorSLCGadget.h"
#include "GadgetIsmrmrdReadWrite.h"

namespace Gadgetron
{

GtPlusAccumulatorSLCGadget::GtPlusAccumulatorSLCGadget() : prev_slc_(-1), cur_slc_(-1)
{

}

GtPlusAccumulatorSLCGadget::~GtPlusAccumulatorSLCGadget()
{

}

int GtPlusAccumulatorSLCGadget::process_config(ACE_Message_Block* mb)
{
    return BaseClass::process_config(mb);
}

bool GtPlusAccumulatorSLCGadget::
copyBufferForSLC(BufferType& buf, int slc, BufferType& bufSLC)
{
    try
    {
        boost::shared_ptr< std::vector<unsigned int> > dims = buf.get_dimensions();

        boost::shared_ptr< std::vector<unsigned int> > dimsSLC = dims;
        (*dimsSLC)[3] = 1;

        try
        {
            bufSLC.create(dimsSLC);
        }
        catch(...)
        {
            GADGET_DEBUG1("Failed create buffer for SLC \n");
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

        int e2, con, phs, rep, set, seg;

        std::vector<unsigned int> pos(10);

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
                                pos[3] = 0;
                                pos[4] = e2;
                                pos[5] = con;
                                pos[6] = phs;
                                pos[7] = rep;
                                pos[8] = set;
                                pos[9] = seg;
                                int offsetBufferSLC = bufSLC.calculate_offset(pos);

                                // copy the image content
                                memcpy(bufSLC.begin()+offsetBufferSLC, buf.begin()+offsetBuffer, sizeof(std::complex<float>)*RO*E1*CHA);
                            }
                        }
                    }
                }
            }
        }
    }
    catch(...)
    {
        GADGET_DEBUG1("Errors in GtPlusAccumulatorSLCGadget::copyBufferForSLC(...) ... \n");
        return false;
    }

    return true;
}

bool GtPlusAccumulatorSLCGadget::
copyReflectBufferForSLC(ReflectBufferType& buf, int slc, ReflectBufferType& bufSLC)
{
    try
    {
        boost::shared_ptr< std::vector<unsigned int> > dims = buf.get_dimensions();

        boost::shared_ptr< std::vector<unsigned int> > dimsSLC = dims;
        (*dimsSLC)[3] = 1;

        try
        {
            bufSLC.create(dimsSLC);
        }
        catch(...)
        {
            GADGET_DEBUG1("Failed create buffer for SLC \n");
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

        int e2, con, phs, rep, set, seg;

        std::vector<unsigned int> pos(10);

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
                                pos[3] = 0;
                                pos[4] = e2;
                                pos[5] = con;
                                pos[6] = phs;
                                pos[7] = rep;
                                pos[8] = set;
                                pos[9] = seg;
                                int offsetBufferSLC = bufSLC.calculate_offset(pos);

                                // copy the image content
                                memcpy(bufSLC.begin()+offsetBufferSLC, buf.begin()+offsetBuffer, sizeof(unsigned short)*RO*E1*CHA);
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

bool GtPlusAccumulatorSLCGadget::triggerSLC(int slc)
{
    try
    {
        GadgetContainerMessage<GadgetMessageImageArray>* cm1 = 
            new GadgetContainerMessage<GadgetMessageImageArray>();

        GadgetContainerMessage< KSpaceBuffer >* cm2 = 
            new GadgetContainerMessage< KSpaceBuffer >();

        cm1->cont(cm2);

        // copy the kspace data for this SLC
        if ( !copyBufferForSLC(kspaceBuffer_->buffer_, slc, cm2->getObjectPtr()->buffer_) ) 
        {
            GADGET_DEBUG1("Unable to copyBufferForSLC\n");
            cm1->release();
            return false;
        }

        if ( !copyReflectBufferForSLC(kspaceBuffer_->reflect_, slc, cm2->getObjectPtr()->reflect_) ) 
        {
            GADGET_DEBUG1("Unable to copyReflectBufferForSLC\n");
            cm1->release();
            return false;
        }

        // copy the message image array for this SLC
        GadgetMessageImageArray aMessageArraySLC;
        messageImage_->extractMessageImageArrayForSLC(slc, aMessageArraySLC);
        cm1->getObjectPtr()->copy(aMessageArraySLC);

        if (!refBuffer_.empty())
        {
            GADGET_MSG("GtPlusAccumulatorSLCGadget - ref signal found : " << refBuffer_.size());

            BufferType refCurr;
            ReflectBufferType refReflectCurr;
            if ( !fillBuffer(refBuffer_, refCurr, refReflectCurr) )
            {
                GADGET_DEBUG1("fillBuffer(refBuffer_, refCurr, refReflectCurr) failed ... \n");
                cm1->release();
                return GADGET_FAIL;
            }

            if ( !copyBufferForSLC(refCurr, slc, cm2->getObjectPtr()->ref_) ) 
            {
                GADGET_DEBUG1("Unable to copyBufferForSLC(refCurr, slc, cm2->getObjectPtr()->ref_)\n");
                cm1->release();
                return false;
            }

            if ( !copyReflectBufferForSLC(refReflectCurr, slc, cm2->getObjectPtr()->refReflect_) ) 
            {
                GADGET_DEBUG1("Unable to copyBufferForSLC(refReflectCurr, slc, cm2->getObjectPtr()->refReflect_)\n");
                cm1->release();
                return false;
            }
        }

        if (!phaseCorrBuffer_.empty())
        {
            GADGET_MSG("GtPlusAccumulatorSLCGadget - phase correction signal found : " << phaseCorrBuffer_.size());

            BufferType phsCorrCurr;
            ReflectBufferType phsCorrReflectCurr;
            if ( !fillBuffer(phaseCorrBuffer_, phsCorrCurr, phsCorrReflectCurr) )
            {
                GADGET_DEBUG1("fillBuffer(phaseCorrBuffer_, phsCorrCurr, phsCorrReflectCurr) failed ... \n");
                cm1->release();
                return GADGET_FAIL;
            }

            if ( !copyBufferForSLC(phsCorrCurr, slc, cm2->getObjectPtr()->phaseCorr_) ) 
            {
                GADGET_DEBUG1("Unable to copyBufferForSLC(phsCorrCurr, slc, cm2->getObjectPtr()->phaseCorr_)\n");
                cm1->release();
                return false;
            }

            if ( !copyReflectBufferForSLC(phsCorrReflectCurr, slc, cm2->getObjectPtr()->phaseCorrReflect_) ) 
            {
                GADGET_DEBUG1("Unable to copyBufferForSLC(phsCorrReflectCurr, slc, cm2->getObjectPtr()->phaseCorrReflect_)\n");
                cm1->release();
                return false;
            }
        }

        if (!noiseBuffer_.empty())
        {
            GADGET_MSG("GtPlusAccumulatorSLCGadget - noise signal found : " << noiseBuffer_.size());

            BufferType noiseCurr;
            ReflectBufferType tmpBuf;
            if ( !fillBuffer(noiseBuffer_, noiseCurr, tmpBuf) )
            {
                GADGET_DEBUG1("fillBuffer(noiseBuffer_, noiseCurr, tmpBuf) failed ... \n");
                cm1->release();
                return false;
            }

            if ( !copyBufferForSLC(noiseCurr, slc, cm2->getObjectPtr()->noise_) ) 
            {
                GADGET_DEBUG1("Unable to copyBufferForSLC(noiseCurr, slc, cm2->getObjectPtr()->noise_)\n");
                cm1->release();
                return false;
            }
        }

        if (!otherBuffer_.empty())
        {
            GADGET_MSG("GtPlusAccumulatorSLCGadget - other signal found : " << otherBuffer_.size());

            BufferType otherCurr;
            ReflectBufferType tmpBuf;
            if ( !fillBuffer(otherBuffer_, otherCurr, tmpBuf) )
            {
                GADGET_DEBUG1("fillBuffer(otherBuffer_, otherCurr, tmpBuf) failed ... \n");
                cm1->release();
                return false;
            }

            if ( !copyBufferForSLC(otherCurr, slc, cm2->getObjectPtr()->other_) ) 
            {
                GADGET_DEBUG1("Unable to copyBufferForSLC(otherCurr, slc, cm2->getObjectPtr()->other_)\n");
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
        GADGET_DEBUG1("Errors in triggerSLC(slc) ... \n");
        return false;
    }

    return true;
}

int GtPlusAccumulatorSLCGadget::process(GadgetContainerMessage<ISMRMRD::AcquisitionHeader>* m1, 
        GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2)
{
    cur_slc_ = m1->getObjectPtr()->idx.slice;

    BaseClass::process(m1, m2);

    if ( prev_slc_==-1 )
    {
        prev_slc_ = cur_slc_;
    }

   // if a new slice comes, it indicates the previous one is complete and can be sent out
    if ( cur_slc_!=prev_slc_ )
    {
        GADGET_MSG("Slice " << prev_slc_ << " is complete ... ");

        // send out prev slice
        GADGET_MSG("GtPlusAccumulatorSLCGadget - trigger next gadget for SLC " << prev_slc_ << " ... ");
        if ( !triggerSLC(prev_slc_) ) 
        {
            GADGET_DEBUG1("Unable to trigger this slc ... \n");
            return GADGET_FAIL;
        }

        prev_slc_ = cur_slc_;
    }

    return GADGET_OK;
}

int GtPlusAccumulatorSLCGadget::close(unsigned long flags)
{
    // the last slice is still not sent out yet
    if ( !triggered_ )
    {
        GADGET_MSG("GtPlusAccumulatorSLCGadget - trigger next gadget for SLC " << cur_slc_ << " ... ");

        if ( !triggerSLC(cur_slc_) ) 
        {
            GADGET_DEBUG1("Unable to trigger this slc ... \n");
            return GADGET_FAIL;
        }

        triggered_ = true;
    }

    // the base class shall do nothing
    triggered_ = true;
    return BaseClass::close(flags);
}

GADGET_FACTORY_DECLARE(GtPlusAccumulatorSLCGadget)

}
