#pragma once

namespace Gadgetron::Core::Parallel {

    template<class... ARGS>
    TypedBranch<ARGS...>::TypedBranch(const GadgetProperties &props) : Branch(props) {}

    template<class... ARGS>
    void TypedBranch<ARGS...>::process(std::shared_ptr<Channel> input,
                                       const std::map<std::string, std::shared_ptr<Channel>> &output,
                                       std::shared_ptr<Channel> bypass) {

        auto typed_input = TypedInputChannel<ARGS...>(*input, *bypass);

        std::map<std::string, std::shared_ptr<OutputChannel>> typed_output;
        for (auto &pair : output) typed_output[pair.first] = pair.second;

        process(typed_input, typed_output);
    }
}
