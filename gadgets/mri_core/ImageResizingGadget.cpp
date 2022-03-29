/*
*       ImageResizingGadget.cpp
*       Author: Hui Xue & Andrew Dupuis
*/

#include "ImageResizingGadget.h"
#include "Node.h"
#include "Types.h"
#include "log.h"
#include "hoNDArray_elemwise.h"
#include "hoNDImage_util.h"

using namespace Gadgetron::Core;

namespace {

    template<class T>
    Image<T> resize(Image<T> &image, size_t new_RO, double scale_factor_RO, size_t new_E1, double scale_factor_E1, size_t new_E2, double scale_factor_E2 , size_t order_interpolator) {
        auto &header = std::get<ISMRMRD::ImageHeader>(image);
        auto &input_array = std::get<hoNDArray<T>>(image);
        hoNDArray<T> output_array;

        std::vector<size_t> dims;
        input_array.get_dimensions(dims);

        size_t RO = dims[0];
        size_t E1 = dims[1];
        size_t E2 = 1;

        if (new_RO > 0)
        {
            dims[0] = new_RO;
        }
        else if (scale_factor_RO > 0)
        {
            dims[0] = (size_t)(dims[0]*scale_factor_RO + 0.5);
        }

        if (new_E1 > 0)
        {
            dims[1] = new_E1;
        }
        else if (scale_factor_E1 > 0)
        {
            dims[1] = (size_t)(dims[1] * scale_factor_E1 + 0.5);
        }

        if (input_array.get_number_of_dimensions() > 2)
        {
            E2 = dims[2];
            if (new_E2 > 0)
            {
                dims[2] = new_E2;
            }
            else if (scale_factor_E2 > 0)
            {
                dims[2] = (size_t)(dims[2] * scale_factor_E2 + 0.5);
            }
        }

        if (input_array.get_number_of_dimensions() == 2)
        {
            typedef hoNDImage<T, 2> ImageType;
            ImageType input_image(RO, E1, input_array.begin());
            ImageType output_image;

            hoNDBoundaryHandlerBorderValue<ImageType> bhBorderValue(input_image);
            hoNDInterpolatorBSpline<ImageType, 2> interp(input_image, bhBorderValue, order_interpolator);

            Gadgetron::resampleImage(input_image, interp, dims, output_image);

            output_array = output_image;
        }
        else if (input_array.get_number_of_dimensions() == 3)
        {
            typedef hoNDImage<T, 3> ImageType;
            
            ImageType input_image(RO, E1, E2, input_array.begin());
            ImageType output_image;

            hoNDBoundaryHandlerBorderValue<ImageType> bhBorderValue(input_image);
            hoNDInterpolatorBSpline<ImageType, 3> interp(input_image, bhBorderValue, order_interpolator);

            Gadgetron::resampleImage(input_image, interp, dims, output_image);

            output_array = output_image;
        }
        else 
        {
            GERROR_STREAM("ImageResizingGadget, only support 2D or 3D input images ... ");
            return image;
        }

        header.matrix_size[0] = output_array.get_size(0);
        header.matrix_size[1] = output_array.get_size(1);
        header.matrix_size[2] = output_array.get_size(2);
        return image;
    }
}

namespace Gadgetron {

    ImageResizingGadget::ImageResizingGadget(const Context &context, const GadgetProperties &properties)
        : ChannelGadget(context, properties) {}

    void ImageResizingGadget::process(InputChannel<Core::AnyImage> &input, OutputChannel &output) {

        // Lambda, performs image resizing
        auto resizeAndPushImage = [&](auto image){ 
            output.push(resize(image, new_RO, new_E1, new_E2, scale_factor_RO, scale_factor_E1, scale_factor_E2, order_interpolator));   
        };

        // Add all the images from input channel to vector of ImageEntries
        for (auto image : input) {
          visit(resizeAndPushImage, image);
        }
    }
    GADGETRON_GADGET_EXPORT(ImageResizingGadget);
}