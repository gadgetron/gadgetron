/** \file   ImageFinishAttribGadget.h
    \brief  Image finish gadget with meta attributes support.
    \author Hui Xue
*/

#ifndef ImageFinishAttribGadget_H
#define ImageFinishAttribGadget_H

#include "Gadget.h"
#include "hoNDArray.h"
#include "GadgetMRIHeaders.h"
#include "GadgetStreamController.h"
#include "ismrmrd/meta.h"
#include "gadgetron_mricore_export.h"

#include <ismrmrd/ismrmrd.h>
#include <complex>

namespace Gadgetron{

    template <typename T> class EXPORTGADGETSMRICORE ImageFinishAttribGadget : public Gadget3<ISMRMRD::ImageHeader, hoNDArray< T >, ISMRMRD::MetaContainer >
    {
    protected:
        virtual int process(GadgetContainerMessage<ISMRMRD::ImageHeader>* m1, GadgetContainerMessage< hoNDArray< T > >* m2, GadgetContainerMessage<ISMRMRD::MetaContainer>* m3);
    };

    class EXPORTGADGETSMRICORE ImageFinishAttribGadgetUSHORT :
        public ImageFinishAttribGadget<ACE_UINT16>
    {
    public:
        GADGET_DECLARE(ImageFinishAttribGadgetUSHORT);
    };

    class EXPORTGADGETSMRICORE ImageFinishAttribGadgetSHORT :
        public ImageFinishAttribGadget<ACE_INT16>
    {
    public:
        GADGET_DECLARE(ImageFinishAttribGadgetSHORT);
    };

    class EXPORTGADGETSMRICORE ImageFinishAttribGadgetFLOAT :
        public ImageFinishAttribGadget<float>
    {
    public:
        GADGET_DECLARE(ImageFinishAttribGadgetFLOAT);
    };

    class EXPORTGADGETSMRICORE ImageFinishAttribGadgetCPLX :
        public ImageFinishAttribGadget< std::complex<float> >
    {
    public:
        GADGET_DECLARE(ImageFinishAttribGadgetCPLX);
    };
}

#endif //ImageFinishAttribGadget_H
