/** \file   GenericReconOpenReconDemoGadget.h
    \brief  This is the class gadget to demo a set of functions for OpenRecon.
            a) RGB images
            b) Attach contours and landmarks to an image
            c) Report images

            The input data will be "thrown-away" and 
    \author Hui Xue
*/

#pragma once

#include "GenericReconGadget.h"

namespace Gadgetron {

    class EXPORTGADGETSMRICORE GenericReconOpenReconDemoGadget : public GenericReconGadget
    {
    public:
        GADGET_DECLARE(GenericReconOpenReconDemoGadget);

        typedef GenericReconGadget BaseClass;

        GenericReconOpenReconDemoGadget();
        ~GenericReconOpenReconDemoGadget() override;

    protected:

        virtual int process_config(ACE_Message_Block* mb) override;
        virtual int process(Gadgetron::GadgetContainerMessage< IsmrmrdReconData >* m1) override;
        virtual int close(unsigned long flags) override;

        ISMRMRD::AcquisitionHeader acq_header_;
        bool triggered_in_close_;

        void create_rgb_3D_image(IsmrmrdImageArray& im);
        void create_rgb_2D_image(IsmrmrdImageArray& im);
        void create_image_with_contours(IsmrmrdImageArray& im);
        void create_reports(IsmrmrdImageArray& reports);

        void fill_image_header(IsmrmrdImageArray& res, bool is_rgb);
    };
}
