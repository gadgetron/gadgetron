/** \file   GenericReconCartesianFFTGadget.h
    \brief  This is the class gadget for both 2DT and 3DT cartesian FFT reconstruction, working on the IsmrmrdReconData.
    \author Stanislas Rapacchi
*/

#pragma once

#include "GenericReconGadget.h"

namespace Gadgetron {

    /// define the recon status
    template <typename T>
    class GenericReconCartesianFFTObj
    {
    public:

        GenericReconCartesianFFTObj()  = default;
        virtual ~GenericReconCartesianFFTObj() = default;

        // ------------------------------------
        /// recon outputs
        // ------------------------------------
        /// reconstructed images, headers and meta attributes
        IsmrmrdImageArray recon_res_;

        // ------------------------------------
        /// buffers used in the recon
        // ------------------------------------
        /// [RO E1 E2 srcCHA Nor1 Sor1 SLC]
        hoNDArray<T> ref_calib_;
        /// [RO E1 E2 dstCHA Nor1 Sor1 SLC]

        /// reference data ready for coil map estimation
        /// [RO E1 E2 dstCHA Nor1 Sor1 SLC]
        hoNDArray<T> ref_coil_map_;

        /// coil sensitivity map, [RO E1 E2 dstCHA - uncombinedCHA Nor1 Sor1 SLC]
        hoNDArray<T> coil_map_;
    };
}

namespace Gadgetron {

    class GenericReconCartesianFFTGadget : public GenericReconGadget
    {
    public:
        typedef GenericReconGadget BaseClass;
        typedef Gadgetron::GenericReconCartesianFFTObj< std::complex<float> > ReconObjType;

        GenericReconCartesianFFTGadget() = default;
        ~GenericReconCartesianFFTGadget() override = default;

        /// ------------------------------------------------------------------------------------
        /// parameters to control the reconstruction
        /// ------------------------------------------------------------------------------------


    protected:
        // --------------------------------------------------
        // gadget functions
        // --------------------------------------------------
        // default interface function
        int process_config(ACE_Message_Block* mb) override;
        int process(Gadgetron::GadgetContainerMessage< IsmrmrdReconData >* m1) override ;

	// fft and coil combination
        void perform_fft_combine(IsmrmrdReconBit& recon_bit, ReconObjType& recon_obj, size_t encoding);

    };
}
