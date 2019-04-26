#pragma once

#include "Channel.h"
#include "PropertyMixin.h"
#include "Context.h"
#include <boost/dll/alias.hpp>

namespace Gadgetron::Core {

    class Node {
    public:
        virtual ~Node() = default;
        virtual void process(InputChannel& in, OutputChannel& out) = 0;
    };

    class ChannelGadget : public Node, public PropertyMixin {
    public:
            using PropertyMixin::PropertyMixin;
    };

    template<class ...ARGS >
    class TypedChannelGadget : public ChannelGadget {
    public:

        explicit TypedChannelGadget(const GadgetProperties& properties): ChannelGadget(properties) {}

        TypedChannelGadget(const Context& context, const GadgetProperties& properties): ChannelGadget(properties) {

        }

        void process(InputChannel& in, OutputChannel& out)  final  {
            auto typed_input = TypedInputChannel<ARGS...>(in, out);
            this->process(typed_input, out);
        }


        virtual void process(TypedInputChannel<ARGS...> &in, OutputChannel &out) = 0;

    };



}

#define GADGETRON_GADGET_EXPORT(GadgetClass)                                        \
std::unique_ptr<Gadgetron::Core::Node> gadget_factory_##GadgetClass(                \
        const Gadgetron::Core::Context &context,                                    \
        const Gadgetron::Core::GadgetProperties &props                   \
) {                                                                                 \
    return std::make_unique<GadgetClass>(context, props);                           \
}                                                                                   \
                                                                                    \
BOOST_DLL_ALIAS(                                                                    \
        gadget_factory_##GadgetClass,                                               \
        gadget_factory_export_##GadgetClass                                         \
)                                                                                   \


