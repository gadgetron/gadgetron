/** \file   ImageFFTGadget.h
    \brief  This Gadget image to its fft.
    \author Hui Xue
*/

#pragma once

#include "Gadget.h"
#include "hoNDArray.h"
#include "ismrmrd/meta.h"

#include <ismrmrd/ismrmrd.h>
#include <ismrmrd/meta.h>

namespace Gadgetron
{
    class ImageFFTGadget :public Gadget3<ISMRMRD::ImageHeader, hoNDArray< std::complex<float> >, ISMRMRD::MetaContainer >
    {
    public:
        typedef std::complex<float> ValueType;
        typedef hoNDArray< ValueType > ArrayType;

        ImageFFTGadget();
        virtual ~ImageFFTGadget();

    protected:
        virtual int process(GadgetContainerMessage<ISMRMRD::ImageHeader>* m1, GadgetContainerMessage< hoNDArray< ValueType > >* m2, GadgetContainerMessage <ISMRMRD::MetaContainer> * m3);
    };
}
