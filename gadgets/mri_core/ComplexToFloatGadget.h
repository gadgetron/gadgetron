/** \file   ComplexToFloatGadget.h
    \brief  This Gadget converts complex float values to float format.
    \author Hui Xue
*/

#ifndef ComplexToFloatGadget_H_
#define ComplexToFloatGadget_H_

#include "Gadget.h"
#include "hoNDArray.h"
#include "ismrmrd/meta.h"
#include "gadgetron_mricore_export.h"

#include <ismrmrd/ismrmrd.h>

namespace Gadgetron
{
    class EXPORTGADGETSMRICORE ComplexToFloatGadget:public Gadget2<ISMRMRD::ImageHeader, hoNDArray< std::complex<float> > >
    {
    public:

        GADGET_DECLARE(ComplexToFloatGadget);

        typedef std::complex<float> ValueType;

        ComplexToFloatGadget();
        virtual ~ComplexToFloatGadget();

    protected:
        virtual int process(GadgetContainerMessage<ISMRMRD::ImageHeader>* m1, GadgetContainerMessage< hoNDArray< ValueType > >* m2);
    };
}

#endif // ComplexToFloatGadget
