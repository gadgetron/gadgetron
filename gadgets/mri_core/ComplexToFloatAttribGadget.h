/** \file   ComplexToFloatAttribGadget.h
    \brief  This Gadget converts complex float values to float format.
    \author Hui Xue
*/

#ifndef ComplexToFloatAttribGadget_H_
#define ComplexToFloatAttribGadget_H_

#include "Gadget.h"
#include "hoNDArray.h"
#include "ismrmrd_meta.h"
#include "gadgetron_mricore_export.h"

#include <ismrmrd.h>

namespace Gadgetron
{
    class EXPORTGADGETSMRICORE ComplexToFloatAttribGadget:public Gadget3<ISMRMRD::ImageHeader, hoNDArray< std::complex<float> >, ISMRMRD::MetaContainer >
    {
    public:

        GADGET_DECLARE(ComplexToFloatAttribGadget);

        typedef std::complex<float> ValueType;

        ComplexToFloatAttribGadget();
        virtual ~ComplexToFloatAttribGadget();

    protected:
        virtual int process(GadgetContainerMessage<ISMRMRD::ImageHeader>* m1, GadgetContainerMessage< hoNDArray< ValueType > >* m2, GadgetContainerMessage<ISMRMRD::MetaContainer>* m3);
    };
}

#endif // ComplexToFloatAttribGadget
