#include "AugmentImageMetadataGadget.h"

#include "log.h"

Gadgetron::AugmentImageMetadataGadget::AugmentImageMetadataGadget(const Core::Context& context, const Core::GadgetProperties& props)
    : PureGadget(context,props)
{}

Gadgetron::Core::Image<std::complex<float>> Gadgetron::AugmentImageMetadataGadget::process_function(
    Gadgetron::Core::Image<std::complex<float>> input_image) const {

    auto& header     = std::get<ISMRMRD::ImageHeader>(input_image);
    auto& meta       = std::get<2>(input_image);
    auto& input_data = std::get<hoNDArray<std::complex<float>>>(input_image);

    if (!meta.has_value()) {
        meta = ISMRMRD::MetaContainer();
    }

    meta->append("ImageRowDir", header.read_dir[0]);
    meta->append("ImageRowDir", header.read_dir[1]);
    meta->append("ImageRowDir", header.read_dir[2]);

    meta->append("ImageColumnDir", header.phase_dir[0]);
    meta->append("ImageColumnDir", header.phase_dir[1]);
    meta->append("ImageColumnDir", header.phase_dir[2]);

    return { header, input_data, meta };
}

namespace Gadgetron{
    GADGETRON_GADGET_EXPORT(AugmentImageMetadataGadget)
}