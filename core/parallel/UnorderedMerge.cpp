#include "UnorderedMerge.h"

#include <thread>

namespace {

    using namespace Gadgetron::Core;

    void move_input_to_output(GenericInputChannel input, OutputChannel output) {
        std::transform(
                std::begin(input),
                std::end(input),
                std::begin(output),
                [](auto message) { return message; }
        );
    }
}

namespace Gadgetron::Core::Parallel {

    UnorderedMerge::UnorderedMerge(const Context &, const GadgetProperties &props) : Merge(props) {}

    void UnorderedMerge::process(std::map<std::string, GenericInputChannel> input, OutputChannel output) {

        std::vector<std::thread> threads;

        for (auto &pair : input) {
            threads.emplace_back(
                    std::thread(
                            move_input_to_output,
                            split(pair.second),
                            split(output)
                    )
            );
        }

        for (auto &thread : threads) thread.join();
    }

    GADGETRON_MERGE_EXPORT(UnorderedMerge)
}

