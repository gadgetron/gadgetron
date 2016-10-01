#include "DependencyDataSenderGadget.h"
#include "ismrmrd/xml.h"

#ifdef USE_OMP
    #include "omp.h"
#endif // USE_OMP

namespace Gadgetron {

    DependencyDataSenderGadget::DependencyDataSenderGadget() : coil_sen_head_gadget_(NULL), scc_head_gadget_(NULL)
    {
    }

    DependencyDataSenderGadget::~DependencyDataSenderGadget()
    {
    }

    int DependencyDataSenderGadget::process_config( ACE_Message_Block* mb )
    {
        coil_sen_head_gadget_ = this->get_controller()->find_gadget(coil_sen_head_gadget.value());
        scc_head_gadget_ = this->get_controller()->find_gadget(scc_head_gadget.value());

        GADGET_CHECK_RETURN( (coil_sen_head_gadget_!=NULL || scc_head_gadget_!=NULL), GADGET_FAIL);

        return GADGET_OK;
    }

    int DependencyDataSenderGadget::process( GadgetContainerMessage<ISMRMRD::AcquisitionHeader>* m1, GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2 )
    {
        bool is_scc = m1->getObjectPtr()->isFlagSet(ISMRMRD::ISMRMRD_ACQ_IS_SURFACECOILCORRECTIONSCAN_DATA);

        if(is_scc)
        {
            if (scc_head_gadget_->putq( m1 ) == -1)
            {
                GERROR( "DependencyDataSenderGadget::process, passing data on to scc gadget" );
                return GADGET_FAIL;
            }
        }
        else
        {
            if (coil_sen_head_gadget_->putq( m1 ) == -1)
            {
                GERROR( "DependencyDataSenderGadget::process, passing data on to coil sen gadget" );
                return GADGET_FAIL;
            }
        }

        return GADGET_OK;
    }

    GADGET_FACTORY_DECLARE( DependencyDataSenderGadget )
}
