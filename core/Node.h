#pragma once

#include "Channel.h"
#include <thread>
#include <memory>
#include <future>
#include "log.h"
#include <ismrmrd/xml.h>
#include "PropertyMixin.h"
#include "Context.h"
#include <boost/dll/alias.hpp>

namespace Gadgetron::Core {

    class Node {
    public:
        virtual ~Node() = default;
        virtual void process(InputChannel& in, OutputChannel& out) = 0;
    };

    class GadgetNode : public Node, public PropertyMixin {
    public:
            GadgetNode(
                    const GadgetProperties& properties
            ) : PropertyMixin(properties) {};
            virtual ~GadgetNode() = default;
    };

    template<class ...ARGS >
    class TypedGadgetNode : public GadgetNode {
    public:
        TypedGadgetNode(const GadgetProperties& properties): GadgetNode(properties) {

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


