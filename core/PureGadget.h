#pragma once
#include "Node.h"

namespace Gadgetron::Core {
class PureGadget : public GadgetNode {
public:
    explicit PureGadget(const GadgetProperties& props)
        : GadgetNode(props)
    {
    }

    void process(InputChannel& in, OutputChannel& out) final
    {
        for (auto message : in )
            out.push(this->process_function(std::move(message)));
    }

    virtual Message process_function(Message) const = 0;
};

template <class RETURN, class... ARGS>
class TypedPureGadget : public PureGadget {
public:

    Message process_function(Message message) const override{
        if (!convertible_to<ARGS...>(message)) return message;
        return  process_function(force_unpack<ARGS...>(message));
    }
    virtual RETURN process_function(ARGS&&... args) const = 0;

};
}
