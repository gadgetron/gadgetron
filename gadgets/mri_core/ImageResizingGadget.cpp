/*
*       ImageResizingGadget.cpp
*       Author: Hui Xue
*/

#include "ImageResizingGadget.h"
#include "hoNDArray_elemwise.h"
#include "hoNDImage_util.h"

namespace Gadgetron
{
    ImageResizingGadget::ImageResizingGadget()
    {
    }

    ImageResizingGadget::~ImageResizingGadget()
    {
    }

    int ImageResizingGadget::process(GadgetContainerMessage<ISMRMRD::ImageHeader>* m1, GadgetContainerMessage< hoNDArray< ValueType > >* m2)
    {
        ArrayType* input_array = m2->getObjectPtr();

        ArrayType output_array;

        std::vector<size_t> dims;
        input_array->get_dimensions(dims);

        size_t RO = dims[0];
        size_t E1 = dims[1];
        size_t E2 = 1;

        if (this->new_RO.value() > 0)
        {
            dims[0] = this->new_RO.value();
        }
        else if (this->scale_factor_RO.value()>0)
        {
            dims[0] = (size_t)(dims[0]*this->scale_factor_RO.value() + 0.5);
        }

        if (this->new_E1.value() > 0)
        {
            dims[1] = this->new_E1.value();
        }
        else if (this->scale_factor_E1.value() > 0)
        {
            dims[1] = (size_t)(dims[1] * this->scale_factor_E1.value() + 0.5);
        }

        if (input_array->get_number_of_dimensions() > 2)
        {
            E2 = dims[2];
            if (this->new_E2.value() > 0)
            {
                dims[2] = this->new_E2.value();
            }
            else if (this->scale_factor_E2.value() > 0)
            {
                dims[2] = (size_t)(dims[2] * this->scale_factor_E2.value() + 0.5);
            }
        }

        if (input_array->get_number_of_elements() == RO*E1)
        {
            typedef hoNDImage<ValueType, 2> ImageType;
            ImageType input_image(RO, E1, input_array->begin()), output_image;

            hoNDBoundaryHandlerBorderValue<ImageType> bhBorderValue(input_image);
            hoNDInterpolatorBSpline<ImageType, 2> interp(input_image, bhBorderValue, this->order_interpolator.value());

            Gadgetron::resampleImage(input_image, interp, dims, output_image);

            output_array = output_image;
        }
        else if (input_array->get_number_of_elements() == RO*E1*E2)
        {
            typedef hoNDImage<ValueType, 3> ImageType;
            ImageType input_image(RO, E1, E2, input_array->begin()), output_image;

            hoNDBoundaryHandlerBorderValue<ImageType> bhBorderValue(input_image);
            hoNDInterpolatorBSpline<ImageType, 3> interp(input_image, bhBorderValue, this->order_interpolator.value());

            Gadgetron::resampleImage(input_image, interp, dims, output_image);

            output_array = output_image;
        }
        else 
        {
            GERROR_STREAM("ImageResizingGadget, only support 2D or 3D input images ... ");
            std::ostringstream ostr;
            ostr << "[";
            for (auto i=0; i<dims.size(); i++) ostr << " " << dims[i];
            ostr << "]";
            GERROR_STREAM("ImageResizingGadget, image size is " << ostr.str());

            if (this->next()->putq(m1) < 0)
            {
                GERROR_STREAM("ImageResizingGadget, failed to pass images to next gadget ... ");
                return GADGET_FAIL;
            }

            return GADGET_OK;
        }

        *m2->getObjectPtr() = output_array;

        m1->getObjectPtr()->matrix_size[0] = output_array.get_size(0);
        m1->getObjectPtr()->matrix_size[1] = output_array.get_size(1);
        m1->getObjectPtr()->matrix_size[2] = output_array.get_size(2);

        if (this->next()->putq(m1) < 0)
        {
            GERROR_STREAM("ImageResizingGadget, failed to pass images to next gadget ... ");
            return GADGET_FAIL;
        }

        return GADGET_OK;
    }

    GADGET_FACTORY_DECLARE(ImageResizingGadget)
}
