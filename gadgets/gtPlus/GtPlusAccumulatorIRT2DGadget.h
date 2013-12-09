#pragma once

#include "GtPlusAccumulatorGadget.h"

namespace Gadgetron
{

class EXPORTGTPLUS GtPlusAccumulatorIRT2DGadget : public GtPlusAccumulatorGadget
{
public:
    GADGET_DECLARE(GtPlusAccumulatorIRT2DGadget);

    typedef GtPlusAccumulatorGadget BaseClass;

    typedef BaseClass::ValueType ValueType;
    typedef BaseClass::ReadOutBufferType ReadOutBufferType;
    typedef BaseClass::BufferType BufferType;
    typedef BaseClass::ReflectBufferType ReflectBufferType;

    GtPlusAccumulatorIRT2DGadget();
    ~GtPlusAccumulatorIRT2DGadget();

protected:

    virtual int process_config(ACE_Message_Block* mb);

    // here, every 2D kspace is stored and send out for every new repetition
    virtual bool storeImageData(GadgetContainerMessage<ISMRMRD::AcquisitionHeader>* m1, GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2, bool isReflect);

    virtual bool triggerREP(int rep);
    virtual int process(Gadgetron::GadgetContainerMessage< ISMRMRD::AcquisitionHeader >* m1, Gadgetron::GadgetContainerMessage< Gadgetron::hoNDArray< std::complex<float> > > * m2);

    virtual bool copyBufferForREP(BufferType& buf, int rep, BufferType& bufREP);
    virtual bool copyReflectBufferForREP(ReflectBufferType& buf, int rep, ReflectBufferType& bufREP);

    int prev_rep_;
    int cur_rep_;

    int num_scan_buffered_;
};

}
