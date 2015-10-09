/** \file   GenericReconPartialFourierHandlingPOCSGadget.h
    \brief  This is the class gadget for both 2DT and 3DT cartesian partial fourier handling using the POCS.
    \author Hui Xue
*/

#pragma once

#include "GenericReconPartialFourierHandlingGadget.h"

namespace Gadgetron {

    class EXPORTGADGETSMRICORE GenericReconPartialFourierHandlingPOCSGadget : public GenericReconPartialFourierHandlingGadget
    {
    public:
        GADGET_DECLARE(GenericReconPartialFourierHandlingPOCSGadget);

        typedef float real_value_type;
        typedef std::complex<real_value_type> ValueType;
        typedef ValueType T;

        typedef GenericReconPartialFourierHandlingGadget BaseClass;

        GenericReconPartialFourierHandlingPOCSGadget();
        ~GenericReconPartialFourierHandlingPOCSGadget();

        /// ------------------------------------------------------------------------------------
        /// parameters to control the reconstruction
        /// ------------------------------------------------------------------------------------

        GADGET_PROPERTY(partial_fourier_POCS_iters, int, "Number of iterations for POCS PF handling", 6);
        GADGET_PROPERTY(partial_fourier_POCS_thres, double, "Threshold for POSC PF handling", 0.01);
        GADGET_PROPERTY(partial_fourier_POCS_transitBand, int, "Transition band width for POCS PF handling", 24);
        GADGET_PROPERTY(partial_fourier_POCS_transitBand_E2, int, "Transition band width for POCS PF handling for E2 dimension", 16);

        // ------------------------------------------------------------------------------------

    protected:

        virtual int perform_partial_fourier_handling();
    };
}
