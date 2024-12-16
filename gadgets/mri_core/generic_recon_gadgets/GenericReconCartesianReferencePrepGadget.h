/** \file   GenericReconCartesianReferencePrepGadget.h
    \brief  This is the class gadget for both 2DT and 3DT cartesian reconstruction to prepare the reference data, working on the ReconData.
    \author Hui Xue
*/

#pragma once

#include "GenericReconBase.h"

#include "hoNDArray_utils.h"
#include "hoNDArray_elemwise.h"

namespace Gadgetron {

    class GenericReconCartesianReferencePrepGadget : public GenericReconDataBase
    {
    public:
        typedef GenericReconDataBase BaseClass;

        GenericReconCartesianReferencePrepGadget();
        ~GenericReconCartesianReferencePrepGadget();

        /// ------------------------------------------------------------------------------------
        /// parameters to control the reconstruction
        /// ------------------------------------------------------------------------------------

        /// ref preparation
        /// whether to average all N for ref generation
        /// for the interleaved mode, the sampling times will be counted and used for averaging
        /// it is recommended to set N as the interleaved dimension
        GADGET_PROPERTY(average_all_ref_N, bool, "Whether to average all N for ref generation", true);
        /// whether to average all S for ref generation
        GADGET_PROPERTY(average_all_ref_S, bool, "Whether to average all S for ref generation", false);

        /// pick specific N or S for ref, these options overwrites average_all_ref_N and average_all_ref_S
        GADGET_PROPERTY(N_for_ref, int, "If N_for_ref >=0, this N will be used for ref preparation", -1);
        GADGET_PROPERTY(S_for_ref, int, "If S_for_ref >=0, this S will be used for ref preparation", -1);

        /// some reconstruction will benefit to fill back the ref data into data array for the embedded mode
        GADGET_PROPERTY(ref_fill_into_data_embedded, bool, "If true and calibration is in embedded mode, fill the full sampled data from ref array into the data array", false);

        /// whether to update ref for every incoming ReconData; for some applications, we may want to only compute ref data once
        /// if false, the ref will only be prepared for the first incoming ReconData
        GADGET_PROPERTY(prepare_ref_always, bool, "Whether to prepare ref for every incoming ReconData", true);

    protected:

        // --------------------------------------------------
        // variables for protocol
        // --------------------------------------------------

        /// indicate whether ref has been prepared for an encoding space
        std::vector<bool> ref_prepared_;

        // for every encoding space
        // calibration mode
        std::vector<mrd::CalibrationMode> calib_mode_;

        // --------------------------------------------------
        // variable for recon
        // --------------------------------------------------

        // --------------------------------------------------
        // gadget functions
        // --------------------------------------------------
        // default interface function
        virtual int process_config(const mrd::Header& header);
        virtual int process(Gadgetron::GadgetContainerMessage< mrd::ReconData >* m1);
    };
}
