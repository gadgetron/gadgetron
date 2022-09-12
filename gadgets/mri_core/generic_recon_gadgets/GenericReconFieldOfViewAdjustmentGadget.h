/** \file   GenericReconFieldOfViewAdjustmentGadget.h
    \brief  This is the class gadget for both 2DT and 3DT reconstruction, working on the IsmrmrdImageArray.
            This gadget will adjust the image field-of-view and/or image size to the protocol prescribed values.

            This class is a part of general cartesian recon chain.

\author     Hui Xue
*/

#pragma once

#include "GenericReconBase.h"

namespace Gadgetron {

    class EXPORTGADGETSMRICORE GenericReconFieldOfViewAdjustmentGadget : public GenericReconImageBase
    {
    public:
        GADGET_DECLARE(GenericReconFieldOfViewAdjustmentGadget);

        typedef GenericReconImageBase BaseClass;

        GenericReconFieldOfViewAdjustmentGadget();
        ~GenericReconFieldOfViewAdjustmentGadget();

    protected:

        // --------------------------------------------------
        // variables for protocol
        // --------------------------------------------------

        // encoding FOV and recon FOV
        std::vector< std::vector<float> > encoding_FOV_;
        std::vector< std::vector<float> > recon_FOV_;
        // recon size
        std::vector< std::vector<size_t> > recon_size_;

        // --------------------------------------------------
        // variable for recon
        // --------------------------------------------------

        // kspace filter
        hoNDArray< std::complex<float> > filter_RO_;
        hoNDArray< std::complex<float> > filter_E1_;
        hoNDArray< std::complex<float> > filter_E2_;

        // kspace buffer
        hoNDArray< std::complex<float> > kspace_buf_;

        // results of filtering
        hoNDArray< std::complex<float> > res_;

        // number of times the process function is called
        size_t process_called_times_;

        // --------------------------------------------------
        // functional functions
        // --------------------------------------------------

        // default interface function
        virtual int process_config(ACE_Message_Block* mb);
        virtual int process(Gadgetron::GadgetContainerMessage< IsmrmrdImageArray >* m1);



        // adjust FOV
        int adjust_FOV(IsmrmrdImageArray& data);

        // perform fft or ifft
        void perform_fft(size_t E2, const hoNDArray< std::complex<float> >& input, hoNDArray< std::complex<float> >& output);
        void perform_ifft(size_t E2, const hoNDArray< std::complex<float> >& input, hoNDArray< std::complex<float> >& output);
    };
}
