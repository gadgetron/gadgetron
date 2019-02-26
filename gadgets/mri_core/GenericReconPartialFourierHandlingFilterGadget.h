/** \file   GenericReconPartialFourierHandlingFilterGadget.h
    \brief  This is the class gadget for both 2DT and 3DT cartesian partial fourier handling using the kspace filter.
    \author Hui Xue
*/

#pragma once

#include "GenericReconPartialFourierHandlingGadget.h"

namespace Gadgetron {

    class GenericReconPartialFourierHandlingFilterGadget : public GenericReconPartialFourierHandlingGadget
    {
    public:

        typedef float real_value_type;
        typedef std::complex<real_value_type> ValueType;
        typedef ValueType T;

        typedef GenericReconPartialFourierHandlingGadget BaseClass;

        using GenericReconPartialFourierHandlingGadget::GenericReconPartialFourierHandlingGadget;

        ~GenericReconPartialFourierHandlingFilterGadget() override =default;

        /// ------------------------------------------------------------------------------------
        /// parameters to control the reconstruction
        /// ------------------------------------------------------------------------------------

        // ------------------------------------------------------------------------------------

        NODE_PROPERTY(partial_fourier_filter_RO_width, double, "Partial fourier filter width for tapered hanning for RO dimension", 0.15);
        NODE_PROPERTY(partial_fourier_filter_E1_width, double, "Partial fourier filter width for tapered hanning for E1 dimension", 0.15);
        NODE_PROPERTY(partial_fourier_filter_E2_width, double, "Partial fourier filter width for tapered hanning for E2 dimension", 0.15);

        NODE_PROPERTY(partial_fourier_filter_densityComp, bool, "Whether to apply density compensation for RO dimension", false);

        // ------------------------------------------------------------------------------------

    protected:

        // partial fourier filters, avoid recomputing. Must be mutable and locked to respect PureGadgets promise of being thread safe

        mutable std::mutex filter_mutex;
        mutable hoNDArray<T> filter_pf_RO_;
        mutable hoNDArray<T> filter_pf_E1_;
        mutable hoNDArray<T> filter_pf_E2_;

        hoNDArray <std::complex<float>> perform_partial_fourier_handling(const hoNDArray <std::complex<float>> &kspace, size_t start_RO, size_t end_RO,
                                         size_t start_E1, size_t end_E1, size_t start_E2, size_t end_E2) const override ;
    };
}
