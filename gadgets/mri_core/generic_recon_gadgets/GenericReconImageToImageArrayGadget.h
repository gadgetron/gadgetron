/** \file   GenericReconImageToImageArrayGadget.h
    \brief  The conversion gadget, from incoming images to IsmrmrdImageArray
    \author Hui Xue
*/

#pragma once

#include <complex>
#include "Gadget.h"
#include "ismrmrd/ismrmrd.h"
#include "ismrmrd/meta.h"
#include "ismrmrd/xml.h"
#include "hoNDArray.h"

#include "mri_core_def.h"

namespace Gadgetron {

class GenericReconImageToImageArrayGadget : public Gadget3<ISMRMRD::ImageHeader, hoNDArray<std::complex<float>>, ISMRMRD::MetaContainer>
{
public:
    typedef std::complex<float> ValueType;
    typedef Gadget3< ISMRMRD::ImageHeader, hoNDArray< ValueType >, ISMRMRD::MetaContainer > BaseClass;

    GenericReconImageToImageArrayGadget();
    ~GenericReconImageToImageArrayGadget();

    virtual int close(unsigned long flags);

    GADGET_PROPERTY(verbose, bool, "Whether to print more information", false);

protected:

    virtual int process_config(ACE_Message_Block* mb);
    virtual int process(GadgetContainerMessage<ISMRMRD::ImageHeader>* m1, GadgetContainerMessage< hoNDArray<ValueType> >* m2, GadgetContainerMessage<ISMRMRD::MetaContainer>* m3);

    int process_called_times_;
};

}
