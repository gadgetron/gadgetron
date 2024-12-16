/** \file   ComplexToFloatGadget.h
    \brief  This Gadget converts complex float values to float format.
    \author Hui Xue
*/

#pragma once

#include "Gadget.h"
#include "hoNDArray.h"
#include "ismrmrd/meta.h"

#include <ismrmrd/ismrmrd.h>

namespace Gadgetron
{
    class ImageResizingGadget :public Gadget2<ISMRMRD::ImageHeader, hoNDArray< std::complex<float> > >
    {
    public:
        typedef std::complex<float> ValueType;
        typedef hoNDArray< ValueType > ArrayType;

        ImageResizingGadget();
        virtual ~ImageResizingGadget();

        GADGET_PROPERTY(new_RO, size_t, "New image size along RO; if 0, no effect", 0);
        GADGET_PROPERTY(new_E1, size_t, "New image size along E1; if 0, no effect", 0);
        GADGET_PROPERTY(new_E2, size_t, "New image size along E2; if 0, no effect", 0);

        GADGET_PROPERTY(scale_factor_RO, double, "Scale factors; if 0, no effect", 1.0);
        GADGET_PROPERTY(scale_factor_E1, double, "Scale factors; if 0, no effect", 1.0);
        GADGET_PROPERTY(scale_factor_E2, double, "Scale factors; if 0, no effect", 1.0);

        // BSpline interpolation was used
        GADGET_PROPERTY(order_interpolator, size_t, "Order of interpolator; higher order may increase noise level", 5);

    protected:
        virtual int process(GadgetContainerMessage<ISMRMRD::ImageHeader>* m1, GadgetContainerMessage< hoNDArray< ValueType > >* m2);
    };
}
