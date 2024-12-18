#pragma once

#include "hoNDArray.h"

#include "mrd/types.h"

namespace Gadgetron::Grappa {

    struct ChannelAnnotation {
        uint64_t number_of_combined_channels, number_of_uncombined_channels;
        std::vector<uint64_t> reordering;
    };

    /** TODO: The ChannelReorderer is no longer used, so all Gadgets using this
     *      AnnotatedAcquisition can be updated to just use mrd::Acquisition.
     */
    using AnnotatedAcquisition = std::tuple<mrd::Acquisition, std::optional<ChannelAnnotation>>;
}
