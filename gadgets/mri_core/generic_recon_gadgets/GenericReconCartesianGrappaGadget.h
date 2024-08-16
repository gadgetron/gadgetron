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
        /// [RO E1 E2 srcCHA Nor1 Sor1 SLC]
        hoNDArray<T> ref_calib_;
        /// [RO E1 E2 dstCHA Nor1 Sor1 SLC]
        hoNDArray<T> ref_calib_dst_;

        /// reference data ready for coil map estimation
        /// [RO E1 E2 dstCHA Nor1 Sor1 SLC]
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
        ~GenericReconCartesianGrappaGadget() override;

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

        /// ------------------------------------------------------------------------------------
        /// down stream coil compression
        /// if downstream_coil_compression==true, down stream coil compression is used
        /// if downstream_coil_compression_num_modesKept > 0, this number of channels will be used as the dst channels
        /// if downstream_coil_compression_num_modesKept==0 and downstream_coil_compression_thres>0, the number of dst channels will be determined  by this threshold
        GADGET_PROPERTY(downstream_coil_compression, bool, "Whether to perform downstream coil compression", true);
        GADGET_PROPERTY(downstream_coil_compression_thres, double, "Threadhold for downstream coil compression", 0.002);
        GADGET_PROPERTY(downstream_coil_compression_num_modesKept, size_t, "Number of modes to keep for downstream coil compression", 0);

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
        virtual int process_config(ACE_Message_Block* mb) override;
        virtual int process(Gadgetron::GadgetContainerMessage< IsmrmrdReconData >* m1) override;
        virtual int close(unsigned long flags) override;

        // --------------------------------------------------
        // recon step functions
        // --------------------------------------------------

        // if downstream coil compression is used, determine number of channels used and prepare the ref_calib_dst_
        virtual void prepare_down_stream_coil_compression_ref_data(const hoNDArray< std::complex<float> >& ref_src, hoNDArray< std::complex<float> >& ref_coil_map, hoNDArray< std::complex<float> >& ref_dst, size_t encoding);

        // calibration, if only one dst channel is prescribed, the GrappaOne is used
        virtual void perform_calib(IsmrmrdReconBit& recon_bit, ReconObjType& recon_obj, size_t encoding);

        // unwrapping or coil combination
        virtual void perform_unwrapping(IsmrmrdReconBit& recon_bit, ReconObjType& recon_obj, size_t encoding);

        // compute snr map
        virtual void compute_snr_map(ReconObjType& recon_obj, hoNDArray< std::complex<float> >& snr_map);

    };
}
