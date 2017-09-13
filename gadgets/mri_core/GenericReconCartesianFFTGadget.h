/** \file   GenericReconCartesianFFTGadget.h
    \brief  This is the class gadget for both 2DT and 3DT cartesian FFT reconstruction, working on the IsmrmrdReconData.
    \author Stanislas Rapacchi
*/

#pragma once

#include "GenericReconGadget.h"

namespace Gadgetron {

    /// define the recon status
    template <typename T>
    class EXPORTGADGETSMRICORE GenericReconCartesianFFTObj
    {
    public:

        GenericReconCartesianFFTObj() {}
        virtual ~GenericReconCartesianFFTObj() {}

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
        hoNDArray<T> ref_calib_dst_;

        /// reference data ready for coil map estimation
        /// [RO E1 E2 dstCHA Nor1 Sor1 SLC]
        hoNDArray<T> ref_coil_map_;

        /// coil sensitivity map, [RO E1 E2 dstCHA - uncombinedCHA Nor1 Sor1 SLC]
        hoNDArray<T> coil_map_;
    };
}

namespace Gadgetron {

    class EXPORTGADGETSMRICORE GenericReconCartesianFFTGadget : public GenericReconGadget
    {
    public:
        GADGET_DECLARE(GenericReconCartesianFFTGadget);

        typedef GenericReconGadget BaseClass;
        typedef Gadgetron::GenericReconCartesianFFTObj< std::complex<float> > ReconObjType;

        GenericReconCartesianFFTGadget();
        ~GenericReconCartesianFFTGadget();

        /// ------------------------------------------------------------------------------------
        /// parameters to control the reconstruction
        /// ------------------------------------------------------------------------------------


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

	// fft and coil combination
        virtual void perform_fft_combine(IsmrmrdReconBit& recon_bit, ReconObjType& recon_obj, size_t encoding);
       
    };
}
