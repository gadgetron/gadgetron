/** \file   GenericReconCartesianSpiritGadget.h
    \brief  This is the class gadget for both 2DT and 3DT cartesian Spirit reconstruction, working on the ReconData.
    \author Hui Xue
*/

#pragma once

#include "GenericReconGadget.h"

namespace Gadgetron {

    /// define the recon status
    template <typename T>
    class GenericReconCartesianSpiritObj
    {
    public:

        GenericReconCartesianSpiritObj() {}
        virtual ~GenericReconCartesianSpiritObj() {}

        // ------------------------------------
        /// recon outputs
        // ------------------------------------
        /// reconstructed images, headers and meta attributes
        mrd::ImageArray recon_res_;

        /// full kspace reconstructed
        hoNDArray<T> full_kspace_;

        // ------------------------------------
        /// buffers used in the recon
        // ------------------------------------
        /// [RO E1 E2 dstCHA Nor1 Sor1 SLC]
        hoNDArray<T> ref_calib_;

        /// reference data ready for coil map estimation
        /// [RO E1 E2 dstCHA Nor1 Sor1 SLC]
        hoNDArray<T> ref_coil_map_;

        /// for combined imgae channel
        /// convolution kernel, [convRO convE1 convE2 dstCHA dstCHA Nor1 Sor1 SLC]
        hoNDArray<T> kernel_;
        /// image domain kernel for 2D, [RO E1 dstCHA dstCHA Nor1 Sor1 SLC]
        hoNDArray<T> kernelIm2D_;

        /// due to the iterative nature of SPIRIT method, the complete memory storage of 3D kernel is not feasible
        /// the RO decouplling is used for 3D spirit
        /// image domain kernel 3D, [convE1 convE2 dstCHA dstCHA RO Nor1 Sor1 SLC]
        hoNDArray<T> kernelIm3D_;

        /// coil sensitivity map, [RO E1 E2 dstCHA Nor1 Sor1 SLC]
        hoNDArray<T> coil_map_;

        /// an estimate of gfactor
        /// gfactor, [RO E1 E2 1 N S SLC]
        hoNDArray<typename realType<T>::Type> gfactor_;
    };
}

namespace Gadgetron {

    class GenericReconCartesianSpiritGadget : public GenericReconGadget
    {
    public:
        typedef GenericReconGadget BaseClass;
        typedef Gadgetron::GenericReconCartesianSpiritObj< std::complex<float> > ReconObjType;

        GenericReconCartesianSpiritGadget();
        ~GenericReconCartesianSpiritGadget();

        /// ------------------------------------------------------------------------------------
        /// Spirit parameters
        GADGET_PROPERTY(spirit_kSize_RO, int, "Spirit kernel size RO", 7);
        GADGET_PROPERTY(spirit_kSize_E1, int, "Spirit kernel size E1", 7);
        GADGET_PROPERTY(spirit_kSize_E2, int, "Spirit kernel size E2", 5);
        GADGET_PROPERTY(spirit_reg_lamda, double, "Spirit regularization threshold", 0);
        GADGET_PROPERTY(spirit_calib_over_determine_ratio, double, "Spirit calibration overdermination ratio", 45);
        GADGET_PROPERTY(spirit_iter_max, int, "Spirit maximal number of iterations", 0);
        GADGET_PROPERTY(spirit_iter_thres, double, "Spirit threshold to stop iteration", 0);
        GADGET_PROPERTY(spirit_print_iter, bool, "Spirit print out iterations", false);

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
        virtual int process_config(const mrd::Header& header);
        virtual int process(Gadgetron::GadgetContainerMessage< mrd::ReconData >* m1);

        // --------------------------------------------------
        // recon step functions
        // --------------------------------------------------

        // calibration, if only one dst channel is prescribed, the SpiritOne is used
        virtual void perform_calib(mrd::ReconAssembly& recon_bit, ReconObjType& recon_obj, size_t encoding);

        // unwrapping or coil combination
        virtual void perform_unwrapping(mrd::ReconAssembly& recon_bit, ReconObjType& recon_obj, size_t encoding);

        // perform spirit unwrapping
        // kspace, kerIm, full_kspace: [RO E1 CHA N S SLC]
        void perform_spirit_unwrapping(hoNDArray< std::complex<float> >& kspace, hoNDArray< std::complex<float> >& kerIm, hoNDArray< std::complex<float> >& full_kspace);

        // perform coil combination
        void perform_spirit_coil_combine(ReconObjType& recon_obj);
    };
}
