#pragma once

#include "GtPlusAccumulatorGadget.h"

namespace Gadgetron
{

class EXPORTGTPLUS GtPlusAccumulatorPerfAIFGadget : public GtPlusAccumulatorGadget
{
public:
    GADGET_DECLARE(GtPlusAccumulatorPerfAIFGadget);

    typedef GtPlusAccumulatorGadget BaseClass;

    typedef BaseClass::ValueType ValueType;
    typedef BaseClass::ReadOutBufferType ReadOutBufferType;
    typedef BaseClass::BufferType BufferType;
    typedef BaseClass::ReflectBufferType ReflectBufferType;

    GtPlusAccumulatorPerfAIFGadget();
    ~GtPlusAccumulatorPerfAIFGadget();

protected:

    virtual int process_config(ACE_Message_Block* mb);

    virtual int process(Gadgetron::GadgetContainerMessage< ISMRMRD::AcquisitionHeader >* m1, Gadgetron::GadgetContainerMessage< Gadgetron::hoNDArray< std::complex<float> > > * m2);

    int cur_rep_;
};

}
