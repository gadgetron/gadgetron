#ifndef IMAGEFINISHGADGET_H
#define IMAGEFINISHGADGET_H

#include "Gadget.h"
#include "hoNDArray.h"
#include "GadgetMRIHeaders.h"
#include "GadgetStreamController.h"
#include "gadgetron_mricore_export.h"

#include <ismrmrd/ismrmrd.h>
#include <complex>

namespace Gadgetron{

    class EXPORTGADGETSMRICORE ImageFinishGadget : public Gadget1 < ISMRMRD::ImageHeader >
    {
    protected:
        virtual int process(GadgetContainerMessage<ISMRMRD::ImageHeader>* m1);
    };
}

#endif //IMAGEFINISHGADGET_H
