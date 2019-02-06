#pragma once

namespace Gadgetron::Core::Parallel {

    template<class... ARGS>
    TypedBranch<ARGS...>::TypedBranch(const GadgetProperties &props) : Branch(props) {}

    template<class... ARGS>
    void TypedBranch<ARGS...>::process(InputChannel input,
                                       std::map<std::string, OutputChannel> output,
                                       OutputChannel bypass) {
        auto typed_input = TypedInputChannel<ARGS...>(input, bypass);
        process(typed_input, std::move(output));
    }
}
