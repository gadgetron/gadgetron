#include "GtPlusAccumulatorPerfAIFGadget.h"
#include "GadgetIsmrmrdReadWrite.h"

namespace Gadgetron
{

GtPlusAccumulatorPerfAIFGadget::GtPlusAccumulatorPerfAIFGadget() : cur_rep_(0)
{

}

GtPlusAccumulatorPerfAIFGadget::~GtPlusAccumulatorPerfAIFGadget()
{

}

int GtPlusAccumulatorPerfAIFGadget::process_config(ACE_Message_Block* mb)
{
    return BaseClass::process_config(mb);
}

int GtPlusAccumulatorPerfAIFGadget::process(GadgetContainerMessage<ISMRMRD::AcquisitionHeader>* m1, 
        GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2)
{
    bool bIsKSpace, bIsRef, bIsNoise, bIsPhaseCorr, bIsReflect, bIsOther;
    if ( !checkStatus(m1->getObjectPtr()->flags, m1->getObjectPtr()->number_of_samples, bIsKSpace, bIsRef, bIsNoise, bIsPhaseCorr, bIsReflect, bIsOther) )
    {
        GADGET_DEBUG1("Failed check readout status\n");
        return GADGET_FAIL;
    }

    // Last scan for measurement of the first slice can indicate the number of repetition
    bool is_last_scan_in_slice = ISMRMRD::FlagBit(ISMRMRD::ACQ_LAST_IN_SLICE).isSet(m1->getObjectPtr()->flags);
    if ( is_last_scan_in_slice && m1->getObjectPtr()->idx.slice==0 && !bIsOther )
    {
        GADGET_MSG("Repetition " << cur_rep_ << " is complete ... ");
        cur_rep_++;
    }

    BaseClass::process(m1, m2);

    // if the other data is stored, need to correct the repetition
    if ( bIsOther )
    {
        if ( !otherBuffer_.empty() )
        {
            otherBuffer_[otherBuffer_.size()-1].acqHead_.idx.repetition = cur_rep_;
        }
    }

    return GADGET_OK;
}

GADGET_FACTORY_DECLARE(GtPlusAccumulatorPerfAIFGadget)
}
