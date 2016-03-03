/** \file   OneEncodingGadget.h
    \brief  This is the class gadget to make sure EPI Flash Ref lines are in the same encoding space as the imaging lines.
    \author Hui Xue
*/

#ifndef ONEENCODINGGADGET_H
#define ONEENCODINGGADGET_H

#include "Gadget.h"
#include "hoNDArray.h"
#include "gadgetron_epi_export.h"

#include <ismrmrd/ismrmrd.h>
#include <complex>

namespace Gadgetron {

    class EXPORTGADGETS_EPI OneEncodingGadget :
        public Gadget2<ISMRMRD::AcquisitionHeader, hoNDArray< std::complex<float> > >
    {
    public:
        OneEncodingGadget();
        virtual ~OneEncodingGadget();

    protected:
        virtual int process(GadgetContainerMessage< ISMRMRD::AcquisitionHeader>* m1, GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2);
    };
}

#endif // ONEENCODINGGADGET_H
