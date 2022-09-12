
#include "DenoiseGadget.h"
#include "GadgetronTimer.h"
#include "non_local_bayes.h"
#include "non_local_means.h"

namespace Gadgetron {
    template <class T>
    Gadgetron::hoNDArray<T> Gadgetron::DenoiseGadget::denoise_function(const Gadgetron::hoNDArray<T>& input) const {

        if (denoiser == "non_local_bayes") {
            return Denoise::non_local_bayes(input, image_std, search_radius);
        } else if (denoiser == "non_local_means") {
            return Denoise::non_local_means(input, image_std, search_radius);
        } else {
            throw std::invalid_argument(std::string("DenoiseGadget: Unknown denoiser type: ") + std::string(denoiser));
        }
    }

    DenoiseSupportedTypes DenoiseGadget::process_function(DenoiseSupportedTypes input) const {
        return Core::visit(
            [&,this](auto& image) { return DenoiseSupportedTypes(this->denoise(std::move(image))); }, input);
    }

    template <class T> DenoiseImage<T> DenoiseGadget::denoise(DenoiseImage<T> image) const {
        return DenoiseImage<T>{ std::move(std::get<ISMRMRD::ImageHeader>(image)),
            denoise_function(std::get<hoNDArray<T>>(image)) };
    }

    IsmrmrdImageArray DenoiseGadget::denoise(IsmrmrdImageArray image_array) const {
        auto& input = image_array.data_;
        input       = denoise_function(input);
        return std::move(image_array);
    }


    GADGETRON_GADGET_EXPORT(DenoiseGadget)
}