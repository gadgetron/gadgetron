#pragma once

#include <map>
#include <set>

#include "common/AnnotatedAcquisition.h"

#include "PureGadget.h"
#include "Types.h"

namespace Gadgetron::Grappa {

    class ChannelReorderer : public Core::PureGadget<AnnotatedAcquisition, Core::Acquisition> {
    public:
        ChannelReorderer(const Core::Context &, const std::unordered_map<std::string, std::string> &);

        AnnotatedAcquisition process_function(Core::Acquisition acquisition) const override;

        NODE_PROPERTY(
                uncombined_channels, std::string,
                "Uncombined channels as a comma separated list of channel indices or names (single quoted).", ""
        );

    private:
        const Core::Context context;
        const std::map<std::string, size_t> channel_labels;
        const std::vector<size_t> uncombined_indices;

        std::map<std::string, size_t>
        build_channel_label_map();

        std::vector<size_t>
        parse_uncombined_channels();

        size_t
        parse_uncombined_channel(const std::string &);

        std::vector<uint64_t>
        create_reordering(size_t number_of_channels) const;
    };
}