#pragma once
#include "Message.h"
#include "PureGadget.h"
#include "Loader.h"
#include "Config.h"

namespace Gadgetron::Main::Nodes {
    class PureStream {
    public:
        PureStream(const Config::PureStream&, const Core::Context&, Loader&);
        Core::Message process_function(Core::Message) const;

    private:
        const std::vector<std::unique_ptr<Core::GenericPureGadget>> pure_gadgets;
    };
}
