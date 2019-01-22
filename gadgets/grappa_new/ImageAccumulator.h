#pragma once

#include <chrono>

#include "Channel.h"
#include "Node.h"

namespace Gadgetron::Grappa {

    class ImageAccumulator : public Core::TypedGadgetNode<Core::Acquisition> {
    public:
        ImageAccumulator(const Core::Context &, const std::unordered_map<std::string, std::string> &);

        // TODO: These all belong on the later recon merge node?
        NODE_PROPERTY(image_series, int, "Image series number for output images", 0);
        NODE_PROPERTY(target_coils, int, "Number of target coils for GRAPPA recon", 0);
        NODE_PROPERTY(uncombined_channels, std::string, "Uncombined channels (as a comma separated list of channel indices", "");

        // TODO: This doesn't need to be special? We can handle all of it in 'uncombined_channels'?
        NODE_PROPERTY(uncombined_channels_by_name, std::string, "Uncombined channels (as a comma separated list of channel names", "");

        // TODO: Stop laughing and do something not a damn config flag.
        NODE_PROPERTY(use_gpu, bool, "If true, recon will try to use GPU resources (when available)", false);

        // TODO: Get rid of these pointless things.
        NODE_PROPERTY(device_channels, int, "Number of device channels", 0);

        void process(Core::TypedInputChannel<Core::Acquisition> &in, Core::OutputChannel &out) override;

    private:
        const Core::Context context;
    };
}