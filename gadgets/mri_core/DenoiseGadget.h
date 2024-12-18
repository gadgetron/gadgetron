
#ifndef GADGETRON_DENOISEGADGET_H
#define GADGETRON_DENOISEGADGET_H

#include "Gadget.h"
#include "hoNDArray.h"

#include "PureGadget.h"
#include <string>

namespace Gadgetron {

    using DenoiseSupportedTypes =
        std::variant<mrd::Image<float>, mrd::Image<std::complex<float>>, mrd::ImageArray>;

    class DenoiseGadget
        : public Core::PureGadget<DenoiseSupportedTypes, DenoiseSupportedTypes> {

    public:
        using Core::PureGadget<DenoiseSupportedTypes, DenoiseSupportedTypes>::PureGadget;

        DenoiseSupportedTypes process_function(DenoiseSupportedTypes input) const;
        NODE_PROPERTY(image_std, float, "Standard deviation of the noise in the produced image", 1);
        NODE_PROPERTY(search_radius, int, "Standard deviation of the noise in the produced image", 25);
        NODE_PROPERTY(denoiser, std::string, "Type of denoiser - non_local_means or non_local_bayes", "non_local_bayes");

    protected:
        template <class T>
        mrd::Image<T> denoise(mrd::Image<T> image) const;
        mrd::ImageArray denoise(mrd::ImageArray image_array) const;

        template <class T> hoNDArray<T> denoise_function(const hoNDArray<T>&) const;
    };

}

#endif // GADGETRON_DENOISEGADGET_H
