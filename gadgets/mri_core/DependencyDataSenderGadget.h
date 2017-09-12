
/** \file   DependencyDataSenderGadget.h
    \brief  This gadget is a part of the depedency handling chain.

            This gadget will send coil sen map data and scc data to corresponding chains.

    \author Hui Xue
*/

#pragma once

#include "Gadget.h"
#include "hoNDArray.h"
#include "gadgetron_mricore_export.h"

#include "GadgetStreamController.h"
#include "ismrmrd/ismrmrd.h"
#include "ismrmrd/xml.h"

#include <complex>

namespace Gadgetron {

    class EXPORTGADGETSMRICORE DependencyDataSenderGadget : public Gadget2<ISMRMRD::AcquisitionHeader, hoNDArray< std::complex<float> > >
    {
    public:
        GADGET_DECLARE( DependencyDataSenderGadget );

        DependencyDataSenderGadget();
        virtual ~DependencyDataSenderGadget();

    protected:

        virtual int process_config( ACE_Message_Block* mb );
        virtual int process( GadgetContainerMessage<ISMRMRD::AcquisitionHeader>* m1, GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2 );

        GADGET_PROPERTY(verbose, bool, "Whether to print more information", false);
        GADGET_PROPERTY( coil_sen_head_gadget, std::string, "Head gadget of coil sen branch", "AccCoilSen" );
        GADGET_PROPERTY( scc_head_gadget, std::string, "Head gadget of scc branch", "AccSCC" );

        Gadget* coil_sen_head_gadget_;
        Gadget* scc_head_gadget_;
    };
}
