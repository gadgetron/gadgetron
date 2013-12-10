#pragma once

#include "GtPlusAccumulatorGadget.h"

namespace Gadgetron
{

class EXPORTGTPLUS GtPlusAccumulatorSLCGadget : public GtPlusAccumulatorGadget
{
public:
    GADGET_DECLARE(GtPlusAccumulatorSLCGadget);

    typedef GtPlusAccumulatorGadget BaseClass;

    typedef BaseClass::ValueType ValueType;
    typedef BaseClass::ReadOutBufferType ReadOutBufferType;
    typedef BaseClass::BufferType BufferType;
    typedef BaseClass::ReflectBufferType ReflectBufferType;

    GtPlusAccumulatorSLCGadget();
    ~GtPlusAccumulatorSLCGadget();

    virtual int close(unsigned long flags);

protected:

    virtual int process_config(ACE_Message_Block* mb);

    virtual bool copyBufferForSLC(BufferType& buf, int slc, BufferType& bufSLC);
    virtual bool copyReflectBufferForSLC(ReflectBufferType& buf, int slc, ReflectBufferType& bufSLC);

    virtual bool triggerSLC(int slc);

    virtual int process(Gadgetron::GadgetContainerMessage< ISMRMRD::AcquisitionHeader >* m1, Gadgetron::GadgetContainerMessage< Gadgetron::hoNDArray< std::complex<float> > > * m2);

    int prev_slc_;
    int cur_slc_;
};

}
