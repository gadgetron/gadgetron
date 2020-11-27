#include "PureStream.h"

namespace {
    using namespace Gadgetron::Core;
    using namespace Gadgetron::Server::Connection;

    std::unique_ptr<GenericPureGadget> load_pure_gadget(const Config::Gadget& conf, const Context& context, Loader& loader) {
        auto factory
            = loader.load_factory<Loader::generic_factory<Node>>("gadget_factory_export_", conf.classname, conf.dll);

        auto gadget = factory(context, conf.properties);

        if (dynamic_cast<GenericPureGadget*>(gadget.get())) {
            return std::unique_ptr<GenericPureGadget>(dynamic_cast<GenericPureGadget*>(gadget.release()));
        }

        throw std::runtime_error("Non-pure Gadget \"" + conf.classname + "\" in pure stream.");
    }

    std::vector<std::unique_ptr<GenericPureGadget>> load_pure_gadgets(
        const std::vector<Config::Gadget>& configs,
        const Context& context,
        Loader& loader
    ) {
        std::vector<std::unique_ptr<GenericPureGadget>> gadgets;
        for (auto& config : configs) {
            gadgets.emplace_back(load_pure_gadget(config, context, loader));
        }
        return gadgets;
    }
}

Gadgetron::Server::Connection::Stream::PureStream::PureStream(
    const Gadgetron::Server::Connection::Config::PureStream& conf,
    const Gadgetron::Core::Context& context,
    Loader& loader
) : pure_gadgets{ load_pure_gadgets(conf.gadgets, context, loader) } {}

Gadgetron::Core::Message Gadgetron::Server::Connection::Stream::PureStream::process_function(
    Gadgetron::Core::Message message
) const {
    return std::accumulate(
            pure_gadgets.begin(),
            pure_gadgets.end(),
            std::move(message),
            [](auto& message, auto& gadget) {
                return gadget->process_function(std::move(message));
            }
    );
}
