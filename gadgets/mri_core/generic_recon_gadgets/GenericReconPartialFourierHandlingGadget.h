/** \file   GenericReconPartialFourierHandlingGadget.h
    \brief  This is the class gadget for both 2DT and 3DT cartesian partial fourier handling, working on the IsmrmrdImageArray.

            This class is a part of general cartesian recon chain.

            The sampled kspace region is found from image meta in fields:

            sampling_limits_RO
            sampling_limits_E1
            sampling_limits_E2

    \author Hui Xue
*/

#pragma once


#include "hoNDArray_utils.h"
#include "hoNDArray_elemwise.h"
#include "hoNDFFT.h"

#include "mri_core_partial_fourier.h"
#include "mri_core_data.h"
#include "PureGadget.h"

namespace Gadgetron {

class GenericReconPartialFourierHandlingGadget : public Core::PureGadget<IsmrmrdImageArray,IsmrmrdImageArray>
    {
    public:

        typedef float real_value_type;
        typedef std::complex<real_value_type> ValueType;
        typedef ValueType T;
        using BaseClass = Core::PureGadget<IsmrmrdImageArray,IsmrmrdImageArray>;

        GenericReconPartialFourierHandlingGadget(const Core::Context& context, const Core::GadgetProperties& props);

        virtual ~GenericReconPartialFourierHandlingGadget() = default;

        IsmrmrdImageArray process_function(IsmrmrdImageArray array) const override;
        /// ------------------------------------------------------------------------------------
        /// parameters to control the reconstruction
        /// ------------------------------------------------------------------------------------
        NODE_PROPERTY(skip_processing_meta_field, std::string, "If this meta field exists, pass the incoming image array to next gadget without processing", "Skip_processing_after_recon");
        NODE_PROPERTY(verbose, bool, "Verbose",false);
        NODE_PROPERTY(perform_timing, bool, "Perform timing",false);

    protected:
        size_t num_encoding_spaces;
        // --------------------------------------------------
        // variables for protocol
        // --------------------------------------------------

        // acceleration factor for E1 and E2
        std::vector<double> acceFactorE1_;
        std::vector<double> acceFactorE2_;


        // --------------------------------------------------
        // implementation functions
        // --------------------------------------------------
        virtual hoNDArray <std::complex<float>> perform_partial_fourier_handling(const hoNDArray <std::complex<float>> &kspace_buffer, size_t start_RO, size_t end_RO,
                                         size_t start_E1, size_t end_E1, size_t start_E2, size_t end_E2) const  = 0;

    };
}
