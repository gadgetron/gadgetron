/** \file   GenericReconEigenChannelGadget.h
    \brief  This is the class gadget for both 2DT and 3DT cartesian reconstruction to convert the data into eigen channel, working on the mrd::ReconData.
            If incoming data has the ref, ref data will be used to compute KLT coefficients
    \author Hui Xue
*/

#pragma once

#include "GenericReconBase.h"

#include "hoNDArray_utils.h"
#include "hoNDArray_elemwise.h"
#include "hoNDKLT.h"

namespace Gadgetron {

    class GenericReconEigenChannelGadget : public GenericReconDataBase
    {
    public:
        typedef GenericReconDataBase BaseClass;
        typedef hoNDKLT< std::complex<float> > KLTType;

        GenericReconEigenChannelGadget();
        ~GenericReconEigenChannelGadget();

        /// ------------------------------------------------------------------------------------
        /// parameters to control the reconstruction
        /// ------------------------------------------------------------------------------------

        /// compute KLT coefficients
        /// whether to average all N for coefficient computation
        /// for the interleaved mode, the sampling times will be counted and used for averaging
        GADGET_PROPERTY(average_all_ref_N, bool, "Whether to average all N for ref generation", true);
        /// whether to average all S for coefficient computation
        GADGET_PROPERTY(average_all_ref_S, bool, "Whether to average all S for ref generation", false);

        /// if update_eigen_channel_coefficients==true, every incoming ReconData will be used to compute KLT coefficients
        /// and the older one will be replaced
        /// if update_eigen_channel_coefficients==false, the KLT coefficients will be computed only once for the first incoming ReconData
        GADGET_PROPERTY(update_eigen_channel_coefficients, bool, "Whether to update KLT coefficients for eigen channel computation", false);

        /// optionally, upstream coil compression can be applied
        /// if upstream_coil_compression==true, only kept channels will be sent out to next gadgets and other channels will be removed
        /// no matter whether upstream_coil_compression is true or false, all channels will be converted into eigen channel
        GADGET_PROPERTY(upstream_coil_compression, bool, "Whether to perform upstream coil compression", true);
        /// the logic here is that if upstream_coil_compression_num_modesKept>0, only upstream_coil_compression_num_modesKept channels will be kept
        /// if upstream_coil_compression_num_modesKept<=0 and upstream_coil_compression_thres>0, this threshold will be used to determine how many channels to keep
        /// the first N and first S will be used to compute number of channels to keep
        GADGET_PROPERTY(upstream_coil_compression_thres, double, "Threadhold for upstream coil compression", -1);
        GADGET_PROPERTY(upstream_coil_compression_num_modesKept, int, "Number of modes to keep for upstream coil compression", 0);

    protected:

        // --------------------------------------------------
        // variables for protocol
        // --------------------------------------------------

        // for every encoding space
        // calibration mode
        std::vector<mrd::CalibrationMode> calib_mode_;

        // --------------------------------------------------
        // variable for recon
        // --------------------------------------------------

        // store the KLT coefficients for N, S, SLC at every encoding space
        std::vector< std::vector< std::vector< std::vector< KLTType > > > > KLT_;

        // --------------------------------------------------
        // gadget functions
        // --------------------------------------------------
        // default interface function
        virtual int process_config(const mrd::Header& header);
        virtual int process(Gadgetron::GadgetContainerMessage< mrd::ReconData >* m1);
    };
}
