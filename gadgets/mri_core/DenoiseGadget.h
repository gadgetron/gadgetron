/**
    \brief  Performs denoising with non-local means and non-local Bayes
    \author Original: David Christoffer Hansen
    \author PureGadget Conversion: David Christoffer Hansen
    \test   Tested by: generic_cartesian_cine_denoise.cfg
*/

#pragma once

#include "Gadget.h"
#include "gadgetron_mricore_export.h"
#include "hoNDArray.h"

#include "PureGadget.h"
#include <ismrmrd/ismrmrd.h>
#include <mri_core_data.h>
#include <string>

namespace Gadgetron {

    template <class T> using DenoiseImage = Core::tuple<ISMRMRD::ImageHeader, hoNDArray<T>>;

    using DenoiseSupportedTypes =
        Core::variant<DenoiseImage<float>, DenoiseImage<std::complex<float>>, IsmrmrdImageArray>;

    class EXPORTGADGETSMRICORE DenoiseGadget
        : public Core::PureGadget<DenoiseSupportedTypes, DenoiseSupportedTypes> {

    public:
        using Core::PureGadget<DenoiseSupportedTypes, DenoiseSupportedTypes>::PureGadget;

        DenoiseSupportedTypes process_function(DenoiseSupportedTypes input) const;
        NODE_PROPERTY(image_std, float, "Standard deviation of the noise in the produced image", 1);
        NODE_PROPERTY(search_radius, int, "Standard deviation of the noise in the produced image", 25);
        NODE_PROPERTY(denoiser, std::string, "Type of denoiser - non_local_means or non_local_bayes", "non_local_bayes");

    protected:
        template <class T>
        DenoiseImage<T> denoise(DenoiseImage<T> image) const;
        IsmrmrdImageArray denoise(IsmrmrdImageArray image_array) const;

        template <class T> hoNDArray<T> denoise_function(const hoNDArray<T>&) const;
    };

}
