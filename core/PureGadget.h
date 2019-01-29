#pragma once
#include "Node.h"

namespace Gadgetron::Core {
    class PureGadget : public GadgetNode {
    public:
       PureGadget(const GadgetProperties& props) : GadgetNode(props) {}

    };

   template<class Base,class... ARGS>
class TypedPureGadget : public PureGadget {
public:
    void process(InputChannel &in, OutputChannel &out) override {
        auto self = (const Base*) this;

        auto tin = TypedInputChannel<ARGS...>(in,out);

        for (auto message : tin){
            out.push(self->process_function(message));
        }

    }
};
}
