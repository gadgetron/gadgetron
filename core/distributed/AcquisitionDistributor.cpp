//
// Created by dchansen on 1/16/19.
//

#include "AcquisitionDistributor.h"
#include <map>
namespace {
    std::function<uint16_t(const ISMRMRD::AcquisitionHeader&)> selector_function(const std::string &key) {

        const std::unordered_map<std::string, std::function<uint16_t(const ISMRMRD::AcquisitionHeader&)>> function_map = {
                {"kspace_encode_step_1", [](
                        const ISMRMRD::AcquisitionHeader &header) { return header.idx.kspace_encode_step_1; }},
                {"kspace_encode_step_2", [](
                        const ISMRMRD::AcquisitionHeader &header) { return header.idx.kspace_encode_step_2; }},
                {"average",              [](const ISMRMRD::AcquisitionHeader &header) { return header.idx.average; }},
                {"slice",                [](const ISMRMRD::AcquisitionHeader &header) { return header.idx.slice; }},
                {"contrast",             [](const ISMRMRD::AcquisitionHeader &header) { return header.idx.contrast; }},
                {"phase",                [](const ISMRMRD::AcquisitionHeader &header) { return header.idx.phase; }},
                {"repetition",           [](
                        const ISMRMRD::AcquisitionHeader &header) { return header.idx.repetition; }},
                {"set",                  [](const ISMRMRD::AcquisitionHeader &header) { return header.idx.set; }},
                {"segment",              [](const ISMRMRD::AcquisitionHeader &header) { return header.idx.segment; }},
                {"user_0",               [](const ISMRMRD::AcquisitionHeader &header) { return header.idx.user[0]; }},
                {"user_1",               [](const ISMRMRD::AcquisitionHeader &header) { return header.idx.user[1]; }},
                {"user_2",               [](const ISMRMRD::AcquisitionHeader &header) { return header.idx.user[2]; }},
                {"user_3",               [](const ISMRMRD::AcquisitionHeader &header) { return header.idx.user[3]; }},
                {"user_4",               [](const ISMRMRD::AcquisitionHeader &header) { return header.idx.user[4]; }},
                {"user_5",               [](const ISMRMRD::AcquisitionHeader &header) { return header.idx.user[5]; }},
                {"user_6",               [](const ISMRMRD::AcquisitionHeader &header) { return header.idx.user[6]; }},
                {"user_7",               [](const ISMRMRD::AcquisitionHeader &header) { return header.idx.user[7]; }}
        };

        return function_map.at(key);
    }
}

void Gadgetron::Core::Distributed::AcquisitionDistributor::process(
        Gadgetron::Core::TypedInputChannel<Gadgetron::Core::Acquisition> &input,
        Gadgetron::Core::Distributed::ChannelCreator &creator) {

    std::map<uint16_t,OutputChannel> channels;
    for (Acquisition acq : input) {
        uint16_t  channel_index = 0;
        if (!channels.count(channel_index)) channels.emplace(channel_index,creator.create());
        channels.at(channel_index).push(std::move(acq));
    }


    GDEBUG("AAAAND DONE");
}


Gadgetron::Core::Distributed::AcquisitionDistributor::AcquisitionDistributor(const Gadgetron::Core::Context &context,
                                                                             const Gadgetron::Core::GadgetProperties &props)
        : TypedDistributor(props), selector(selector_function(parallel_dimension)) {

}

namespace Gadgetron::Core::Distributed {
    GADGETRON_DISTRIBUTOR_EXPORT(AcquisitionDistributor);
}