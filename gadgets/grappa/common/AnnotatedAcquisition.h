#pragma once

#include "hoNDArray.h"
#include "Types.h"

namespace Gadgetron::Grappa {

    struct ChannelAnnotation {
        uint64_t number_of_combined_channels, number_of_uncombined_channels;
        std::vector<uint64_t> reordering;
    };

    using AnnotatedAcquisition = Core::tuple<ISMRMRD::AcquisitionHeader,
                                             Core::optional<hoNDArray<float>>,
                                             hoNDArray<std::complex<float>>,
                                             Core::optional<ChannelAnnotation>>;
}
