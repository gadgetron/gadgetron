/**
\file   GenericImageReconArrayToImageGadget.h
\brief  Convert image array to single image
\author Hui Xue
*/

#pragma once

#include "GenericImageReconGadget.h"

namespace Gadgetron {

    class EXPORTGADGETSMRICORE GenericImageReconArrayToImageGadget : public GenericImageReconGadget
    {
    public:
        GADGET_DECLARE(GenericImageReconArrayToImageGadget);

        typedef GenericImageReconGadget BaseClass;

        typedef BaseClass::ValueType ValueType;
        typedef BaseClass::Image2DType ImageType;
        typedef BaseClass::Image2DBufferType ImageBufferType;
        typedef BaseClass::ImgArrayType ImgArrayType;

        GenericImageReconArrayToImageGadget();
        ~GenericImageReconArrayToImageGadget();

        virtual int close(unsigned long flags);

    protected:

        virtual int process_config(ACE_Message_Block* mb);
        virtual int processImageBuffer(ImageBufferType& ori);
    };
}
