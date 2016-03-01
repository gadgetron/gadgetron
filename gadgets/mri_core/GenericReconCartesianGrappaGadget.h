/** \file   GenericReconCartesianGrappaGadget.h
    \brief  This is the class gadget for both 2DT and 3DT cartesian grappa and grappaone reconstruction, working on the IsmrmrdReconData.
    \author Hui Xue
*/

#pragma once

#include "GenericReconGadget.h"

namespace Gadgetron {

    /// define the recon status
    template <typename T>
    class EXPORTGADGETSMRICORE GenericReconCartesianGrappaObj
    {
    public:

        GenericReconCartesianGrappaObj() {}
        virtual ~GenericReconCartesianGrappaObj() {}

        // ------------------------------------
        /// recon outputs
        // ------------------------------------
        /// reconstructed images, headers and meta attributes
        IsmrmrdImageArray recon_res_;

        /// gfactor, [RO E1 E2 uncombinedCHA+1 N S SLC]
        hoNDArray<typename realType<T>::Type> gfactor_;

        // ------------------------------------
        /// buffers used in the recon
        // ------------------------------------
        /// [RO E1 E2 CHA Nor1 Sor1 SLC]
        hoNDArray<T> ref_calib_;

        /// reference data ready for coil map estimation
        /// [RO E1 E2 CHA Nor1 Sor1 SLC]
        hoNDArray<T> ref_coil_map_;

        /// for combined imgae channel
        /// convolution kernel, [RO E1 E2 srcCHA - uncombinedCHA dstCHA - uncombinedCHA Nor1 Sor1 SLC]
        hoNDArray<T> kernel_;
        /// image domain kernel, [RO E1 E2 srcCHA - uncombinedCHA dstCHA - uncombinedCHA Nor1 Sor1 SLC]
        hoNDArray<T> kernelIm_;
        /// image domain unmixing coefficients, [RO E1 E2 srcCHA - uncombinedCHA Nor1 Sor1 SLC]
        hoNDArray<T> unmixing_coeff_;

        /// coil sensitivity map, [RO E1 E2 dstCHA - uncombinedCHA Nor1 Sor1 SLC]
        hoNDArray<T> coil_map_;
    };
}

namespace Gadgetron {

    class EXPORTGADGETSMRICORE GenericReconCartesianGrappaGadget : public GenericReconGadget
    {
    public:
        GADGET_DECLARE(GenericReconCartesianGrappaGadget);

        typedef GenericReconGadget BaseClass;
        typedef Gadgetron::GenericReconCartesianGrappaObj< std::complex<float> > ReconObjType;

        GenericReconCartesianGrappaGadget();
        ~GenericReconCartesianGrappaGadget();

        /// ------------------------------------------------------------------------------------
        /// parameters to control the reconstruction
        /// ------------------------------------------------------------------------------------

        /// ------------------------------------------------------------------------------------
        /// image sending
        GADGET_PROPERTY(send_out_gfactor, bool, "Whether to send out gfactor map", false);
        GADGET_PROPERTY(send_out_snr_map, bool, "Whether to send out SNR map", false);

        /// ------------------------------------------------------------------------------------
        /// Grappa parameters
        GADGET_PROPERTY(grappa_kSize_RO, int, "Grappa kernel size RO", 5);
        GADGET_PROPERTY(grappa_kSize_E1, int, "Grappa kernel size E1", 4);
        GADGET_PROPERTY(grappa_kSize_E2, int, "Grappa kernel size E2", 4);
        GADGET_PROPERTY(grappa_reg_lamda, double, "Grappa regularization threshold", 0.0005);
        GADGET_PROPERTY(grappa_calib_over_determine_ratio, double, "Grappa calibration overdermination ratio", 45);

    protected:

        // --------------------------------------------------
        // variable for recon
        // --------------------------------------------------
        // record the recon kernel, coil maps etc. for every encoding space
        std::vector< ReconObjType > recon_obj_;

        // --------------------------------------------------
        // gadget functions
        // --------------------------------------------------
        // default interface function
        virtual int process_config(ACE_Message_Block* mb);
        virtual int process(Gadgetron::GadgetContainerMessage< IsmrmrdReconData >* m1);

        // --------------------------------------------------
        // recon step functions
        // --------------------------------------------------

        // calibration, if only one dst channel is prescribed, the GrappaOne is used
        virtual void perform_calib(IsmrmrdReconBit& recon_bit, ReconObjType& recon_obj, size_t encoding);

        // unwrapping or coil combination
        virtual void perform_unwrapping(IsmrmrdReconBit& recon_bit, ReconObjType& recon_obj, size_t encoding);
    };
}
