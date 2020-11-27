/** \file   GenericReconNoiseStdMapComputingGadget.h
    \brief  This is the class gadget to compute standard deviation map, working on the IsmrmrdImageArray.

            This class is a part of general cartesian recon chain. It computes the std map on incoming SNR images.

\author     Hui Xue
*/

#pragma once

#include "GenericReconBase.h"

namespace Gadgetron {

    class EXPORTGADGETSMRICORE GenericReconNoiseStdMapComputingGadget : public GenericReconImageBase
    {
    public:
        GADGET_DECLARE(GenericReconNoiseStdMapComputingGadget);

        typedef float real_value_type;
        typedef std::complex<real_value_type> ValueType;
        typedef ValueType T;

        typedef GenericReconImageBase BaseClass;

        GenericReconNoiseStdMapComputingGadget();
        ~GenericReconNoiseStdMapComputingGadget();

        /// ------------------------------------------------------------------------------------
        /// parameters to control the reconstruction
        /// ------------------------------------------------------------------------------------

        /// ------------------------------------------------------------------------------------
        /// start N index to compute std map
        GADGET_PROPERTY(start_N_for_std_map, int, "Start N index to compute std map", 5);

        // ------------------------------------------------------------------------------------

    protected:

        // --------------------------------------------------
        // functional functions
        // --------------------------------------------------

        // default interface function
        virtual int process_config(ACE_Message_Block* mb);
        virtual int process(Gadgetron::GadgetContainerMessage< IsmrmrdImageArray >* m1);

    };
}
