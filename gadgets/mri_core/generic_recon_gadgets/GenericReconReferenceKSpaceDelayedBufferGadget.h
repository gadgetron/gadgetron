/** \file   GenericReconReferenceKSpaceDelayedBufferGadget.h
    \brief  Generic chain does requires that reference data is acquried right before the imaging data, depending on the triggering scheme
            Sometimes for the separate acceleration mode, ref data for all SLC may be acuquired all at once at the beginning of scan

            This gadget will buffer the ref data for every slice and only send them down stream when imaging data for a slice arrives.

    \author Hui Xue
*/

#pragma once

#include "GenericReconBase.h"

#include "hoNDArray_utils.h"
#include "hoNDArray_elemwise.h"

namespace Gadgetron {

    class EXPORTGADGETSMRICORE GenericReconReferenceKSpaceDelayedBufferGadget : public GenericReconKSpaceReadoutBase
    {
    public:
        GADGET_DECLARE(GenericReconReferenceKSpaceDelayedBufferGadget);

        typedef GenericReconKSpaceReadoutBase BaseClass;

        GenericReconReferenceKSpaceDelayedBufferGadget();
        ~GenericReconReferenceKSpaceDelayedBufferGadget();

        /// ------------------------------------------------------------------------------------
        /// parameters to control the reconstruction
        /// ------------------------------------------------------------------------------------

    protected:

        // --------------------------------------------------
        // variables for protocol
        // --------------------------------------------------

        // for every encoding space
        // calibration mode
        std::vector<Gadgetron::ismrmrdCALIBMODE> calib_mode_;

        // number of slices
        std::vector< size_t > SLC_;

        // reference data buffer for every slice, ecoding-slice-acq
        std::vector< std::vector< std::vector< GadgetContainerMessage<ISMRMRD::AcquisitionHeader>* > > > ref_buf_;

        // check has imaging data arrived yet
        std::vector< std::vector<bool> > imaging_data_arrived_;

        // --------------------------------------------------
        // variable for recon
        // --------------------------------------------------

        // ref buffer for every slice

        // --------------------------------------------------
        // gadget functions
        // --------------------------------------------------
        // default interface function
        virtual int process_config(ACE_Message_Block* mb);
        virtual int process(Gadgetron::GadgetContainerMessage< ISMRMRD::AcquisitionHeader >* m1);

        // close call
        int close(unsigned long flags);
    };
}
