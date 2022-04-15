
#include "GenericReconImageToImageArrayGadget.h"

namespace Gadgetron { 

    GenericReconImageToImageArrayGadget::GenericReconImageToImageArrayGadget() : BaseClass(), process_called_times_(0)
    {
    }

    GenericReconImageToImageArrayGadget::~GenericReconImageToImageArrayGadget()
    {
    }

    int GenericReconImageToImageArrayGadget::process_config(ACE_Message_Block* mb)
    {
        // pass the xml file
        ISMRMRD::IsmrmrdHeader h;
        try
        {
          deserialize(mb->rd_ptr(),h);
        }
        catch (...)
        {
          GDEBUG("Error parsing ISMRMRD Header");
          throw;
          return GADGET_FAIL;
        }

        size_t NE = h.encoding.size();
        GDEBUG_CONDITION_STREAM(verbose.value(), "Number of encoding spaces: " << NE);

        return GADGET_OK;
    }

    int GenericReconImageToImageArrayGadget::process(GadgetContainerMessage<ISMRMRD::ImageHeader>* m1, GadgetContainerMessage< hoNDArray<ValueType> >* m2, GadgetContainerMessage<ISMRMRD::MetaContainer>* m3)
    {
        GDEBUG_CONDITION_STREAM(verbose.value(), "GenericReconImageToImageArrayGadget<T, D>::process(...) starts ... ");

        process_called_times_++;

        Gadgetron::GadgetContainerMessage<Gadgetron::IsmrmrdImageArray>* cm1 = new Gadgetron::GadgetContainerMessage<Gadgetron::IsmrmrdImageArray>();

        cm1->getObjectPtr()->headers_.create(1, 1, 1); // must allocate for [N S LOC]
        memcpy(&cm1->getObjectPtr()->headers_(0, 0, 0), m1->getObjectPtr(), sizeof(ISMRMRD::ImageHeader));

        size_t RO = m2->getObjectPtr()->get_size(0);
        size_t E1 = m2->getObjectPtr()->get_size(1);
        size_t E2 = m2->getObjectPtr()->get_size(2);
        size_t CHA = m2->getObjectPtr()->get_size(3);

        cm1->getObjectPtr()->data_.create(RO, E1, E2, CHA, 1, 1, 1);
        memcpy(cm1->getObjectPtr()->data_ .begin(), m2->getObjectPtr()->begin(), m2->getObjectPtr()->get_number_of_bytes());

        cm1->getObjectPtr()->meta_.resize(1);
        cm1->getObjectPtr()->meta_[0] = *m3->getObjectPtr();

        if(verbose.value())
        {
            GDEBUG_STREAM("--> GenericReconImageToImageArrayGadget, receive " << process_called_times_ << " images : [Cha Slice Con Phase Rep Set Ave] : [" 
                << m1->getObjectPtr()->channels << " " << m1->getObjectPtr()->slice << " " << m1->getObjectPtr()->contrast 
                << " " << m1->getObjectPtr()->phase << " " << m1->getObjectPtr()->repetition << " " << m1->getObjectPtr()->set << " " << m1->getObjectPtr()->average << "]");
        }

        if (this->next()->putq(cm1) < 0)
        {
            m1->release();
            cm1->release();
            return GADGET_FAIL;
        }

        m1->release();
        return GADGET_OK;
    }

    int GenericReconImageToImageArrayGadget::close(unsigned long flags)
    {
        GDEBUG_CONDITION_STREAM(true, "GenericReconImageToImageArrayGadget - close(flags) : " << flags);
        if ( BaseClass::close(flags) != GADGET_OK ) return GADGET_FAIL;
        return GADGET_OK;
    }

    GADGET_FACTORY_DECLARE(GenericReconImageToImageArrayGadget)
}
