
#include <map>
#include "AcquisitionDistributor.h"

namespace {

    using Header = ISMRMRD::AcquisitionHeader;
    const std::unordered_map<std::string, std::function<uint16_t(const Header &)>> function_map = {
            {"kspace_encode_step_1", [](const Header &header) { return header.idx.kspace_encode_step_1; }},
            {"kspace_encode_step_2", [](const Header &header) { return header.idx.kspace_encode_step_2; }},
            {"average",              [](const Header &header) { return header.idx.average; }},
            {"slice",                [](const Header &header) { return header.idx.slice; }},
            {"contrast",             [](const Header &header) { return header.idx.contrast; }},
            {"phase",                [](const Header &header) { return header.idx.phase; }},
            {"repetition",           [](const Header &header) { return header.idx.repetition; }},
            {"set",                  [](const Header &header) { return header.idx.set; }},
            {"segment",              [](const Header &header) { return header.idx.segment; }},
            {"user_0",               [](const Header &header) { return header.idx.user[0]; }},
            {"user_1",               [](const Header &header) { return header.idx.user[1]; }},
            {"user_2",               [](const Header &header) { return header.idx.user[2]; }},
            {"user_3",               [](const Header &header) { return header.idx.user[3]; }},
            {"user_4",               [](const Header &header) { return header.idx.user[4]; }},
            {"user_5",               [](const Header &header) { return header.idx.user[5]; }},
            {"user_6",               [](const Header &header) { return header.idx.user[6]; }},
            {"user_7",               [](const Header &header) { return header.idx.user[7]; }}
    };

    const std::function<uint16_t(const ISMRMRD::AcquisitionHeader&)> &selector_function(const std::string &key) {
        return function_map.at(key);
    }
}


namespace Gadgetron::Core::Distributed {

    void AcquisitionDistributor::process(InputChannel<Acquisition> &input,
            ChannelCreator &creator
    ) {
        std::map<uint16_t, OutputChannel> channels{};

        for (Acquisition acq : input) {
            auto index = selector(std::get<Header>(acq));

            if (!channels.count(index)) channels.emplace(index, creator.create());

            channels.at(index).push(std::move(acq));
        }
    }

    AcquisitionDistributor::AcquisitionDistributor(
            const Gadgetron::Core::Context &context,
            const Gadgetron::Core::GadgetProperties &props
    ) : TypedDistributor(props), selector(selector_function(parallel_dimension)) {}

    GADGETRON_DISTRIBUTOR_EXPORT(AcquisitionDistributor);
}