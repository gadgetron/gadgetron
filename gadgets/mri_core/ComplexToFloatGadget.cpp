#include "ComplexToFloatGadget.h"
#include "hoNDArray_elemwise.h"
#include "log.h"

Gadgetron::ComplexToFloatGadget::ComplexToFloatGadget(
    const Gadgetron::Core::Context& context, const Gadgetron::Core::GadgetProperties& props)
    : PureGadget(context,props) {

    converters = { { ISMRMRD::ISMRMRD_IMTYPE_MAGNITUDE, [](const auto& image) { return abs(image); } },
                   { ISMRMRD::ISMRMRD_IMTYPE_PHASE,     [](const auto& image) { return argument(image); } },
                   { ISMRMRD::ISMRMRD_IMTYPE_REAL,      [](const auto& image) { return real(image); } },
                   { ISMRMRD::ISMRMRD_IMTYPE_IMAG,      [](const auto& image) { return imag(image); } }};
};

Gadgetron::Core::Image<float> Gadgetron::ComplexToFloatGadget::process_function(
    Gadgetron::Core::Image<std::complex<float>> input_image) const {

    auto& header     = std::get<ISMRMRD::ImageHeader>(input_image);
    auto& meta       = std::get<2>(input_image);
    auto& input_data = std::get<hoNDArray<std::complex<float>>>(input_image);

    hoNDArray<float> output_data;

    if (converters.count(header.image_type)) {
        output_data = converters.at(header.image_type)(input_data);
    }
    else {
        GDEBUG_STREAM("Image type not set; defaulting to magnitude image.");
        output_data = converters.at(ISMRMRD::ISMRMRD_IMTYPE_MAGNITUDE)(input_data);
        header.image_type = ISMRMRD::ISMRMRD_IMTYPE_MAGNITUDE;
    }

    return { header, output_data, meta };
}

namespace Gadgetron{
    GADGETRON_GADGET_EXPORT(ComplexToFloatGadget)
}