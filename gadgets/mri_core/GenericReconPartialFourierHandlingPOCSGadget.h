/** \file   GenericReconPartialFourierHandlingPOCSGadget.h
    \brief  This is the class gadget for both 2DT and 3DT cartesian partial fourier handling using the POCS.
    \author Hui Xue
*/

#pragma once

#include "GenericReconPartialFourierHandlingGadget.h"

namespace Gadgetron {

    class GenericReconPartialFourierHandlingPOCSGadget : public GenericReconPartialFourierHandlingGadget {
    public:
        typedef float real_value_type;
        typedef std::complex<real_value_type> ValueType;
        typedef ValueType T;

        typedef GenericReconPartialFourierHandlingGadget BaseClass;

        using GenericReconPartialFourierHandlingGadget::GenericReconPartialFourierHandlingGadget;

        /// ------------------------------------------------------------------------------------
        /// parameters to control the reconstruction
        /// ------------------------------------------------------------------------------------

        NODE_PROPERTY(partial_fourier_POCS_iters, size_t, "Number of iterations for POCS PF handling", 6);
        NODE_PROPERTY(partial_fourier_POCS_thres, double, "Threshold for POSC PF handling", 0.01);
        NODE_PROPERTY(partial_fourier_POCS_transitBand, size_t, "Transition band width for POCS PF handling", 24);
        NODE_PROPERTY(partial_fourier_POCS_transitBand_E2, size_t,
            "Transition band width for POCS PF handling for E2 dimension", 16);

        // ------------------------------------------------------------------------------------

    protected:
        hoNDArray<std::complex<float>> perform_partial_fourier_handling(
            const hoNDArray<std::complex<float>>& kspace_buffer, size_t start_RO, size_t end_RO, size_t start_E1,
            size_t end_E1, size_t start_E2, size_t end_E2) const override;
    };
}
