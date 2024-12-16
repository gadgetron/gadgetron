/** \file   CoilComputationGadget.h
    \brief  This Gadget computes coil senstivities based on image data.
    \author Johannes Mayer
*/

#pragma once

#include "Gadget.h"
#include "hoNDArray.h"

#include "ismrmrd/meta.h"
#include <ismrmrd/ismrmrd.h>

namespace Gadgetron {
    class CoilComputationGadget :public Gadget2<ISMRMRD::ImageHeader, hoNDArray< std::complex<float> > >
    {
    public:
        GADGET_PROPERTY(ks_, size_t, "Correlation matrix size in plane.", 7);
        GADGET_PROPERTY(kz_, size_t, "Correlation matrix size in slice direction.", 5);
        GADGET_PROPERTY(power_, size_t, "Number of iterations to apply power method", 3);

        typedef std::complex<float> ValueType;
        typedef hoNDArray< ValueType > ArrayType;

        CoilComputationGadget();
        virtual ~CoilComputationGadget();


    protected:
        virtual int process(GadgetContainerMessage<ISMRMRD::ImageHeader>* m1, GadgetContainerMessage<ArrayType>* m2);
    };
}
