/** \file   GenericReconPartialFourierHandlingFilterGadget.h
    \brief  This is the class gadget for both 2DT and 3DT cartesian partial fourier handling using the kspace filter.
    \author Hui Xue
*/

#pragma once

#include "GenericReconPartialFourierHandlingGadget.h"

namespace Gadgetron {

    class EXPORTGADGETSMRICORE GenericReconPartialFourierHandlingFilterGadget : public GenericReconPartialFourierHandlingGadget
    {
    public:
        GADGET_DECLARE(GenericReconPartialFourierHandlingFilterGadget);

        typedef float real_value_type;
        typedef std::complex<real_value_type> ValueType;
        typedef ValueType T;

        typedef GenericReconPartialFourierHandlingGadget BaseClass;

        GenericReconPartialFourierHandlingFilterGadget();
        ~GenericReconPartialFourierHandlingFilterGadget();

        /// ------------------------------------------------------------------------------------
        /// parameters to control the reconstruction
        /// ------------------------------------------------------------------------------------

        // ------------------------------------------------------------------------------------

        GADGET_PROPERTY(partial_fourier_filter_RO_width, double, "Partial fourier filter width for tapered hanning for RO dimension", 0.15);
        GADGET_PROPERTY(partial_fourier_filter_E1_width, double, "Partial fourier filter width for tapered hanning for E1 dimension", 0.15);
        GADGET_PROPERTY(partial_fourier_filter_E2_width, double, "Partial fourier filter width for tapered hanning for E2 dimension", 0.15);

        GADGET_PROPERTY(partial_fourier_filter_densityComp, bool, "Whether to apply density compensation for RO dimension", false);

        // ------------------------------------------------------------------------------------

    protected:

        // partial fourier filters, avoid recomputing
        hoNDArray<T> filter_pf_RO_;
        hoNDArray<T> filter_pf_E1_;
        hoNDArray<T> filter_pf_E2_;

        virtual int perform_partial_fourier_handling();
    };
}
