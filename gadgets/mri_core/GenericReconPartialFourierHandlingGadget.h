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

#include "GenericReconBase.h"

#include "hoNDArray_utils.h"
#include "hoNDArray_elemwise.h"
#include "hoNDFFT.h"

#include "mri_core_partial_fourier.h"

namespace Gadgetron {

    class EXPORTGADGETSMRICORE GenericReconPartialFourierHandlingGadget : public GenericReconImageBase
    {
    public:
        GADGET_DECLARE(GenericReconPartialFourierHandlingGadget);

        typedef float real_value_type;
        typedef std::complex<real_value_type> ValueType;
        typedef ValueType T;

        typedef GenericReconImageBase BaseClass;

        GenericReconPartialFourierHandlingGadget();
        ~GenericReconPartialFourierHandlingGadget();

        /// ------------------------------------------------------------------------------------
        /// parameters to control the reconstruction
        /// ------------------------------------------------------------------------------------
        GADGET_PROPERTY(skip_processing_meta_field, std::string, "If this meta field exists, pass the incoming image array to next gadget without processing", "Skip_processing_after_recon");

    protected:

        // --------------------------------------------------
        // variables for protocol
        // --------------------------------------------------

        // acceleration factor for E1 and E2
        std::vector<double> acceFactorE1_;
        std::vector<double> acceFactorE2_;

        // --------------------------------------------------
        // variable for recon
        // --------------------------------------------------

        // sampled range
        size_t startRO_;
        size_t endRO_;
        size_t startE1_;
        size_t endE1_;
        size_t startE2_;
        size_t endE2_;

        // kspace buffer
        hoNDArray<T> kspace_buf_;

        // results of pf handling
        hoNDArray<T> pf_res_;

        // --------------------------------------------------
        // functional functions
        // --------------------------------------------------

        // default interface function
        virtual int process_config(ACE_Message_Block* mb);
        virtual int process(Gadgetron::GadgetContainerMessage< IsmrmrdImageArray >* m1);


        // --------------------------------------------------
        // implementation functions
        // --------------------------------------------------
        virtual int perform_partial_fourier_handling() = 0;

    };
}
