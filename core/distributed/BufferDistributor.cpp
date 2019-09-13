#include "BufferDistributor.h"

void Gadgetron::Core::Distributed::BufferDistributor::process(
        Gadgetron::Core::InputChannel<Gadgetron::IsmrmrdReconData> &input,
        Gadgetron::Core::Distributed::ChannelCreator &creator) {

    std::vector<OutputChannel> channels;
    for (size_t i = 0; i < encoding_spaces; i++) channels.push_back(creator.create());

    for (IsmrmrdReconData reconData : input){

        for (size_t i =0; i < reconData.rbit_.size(); i++)
            channels[i].push(IsmrmrdReconData{{std::move(reconData.rbit_[i])}});
    }
}

Gadgetron::Core::Distributed::BufferDistributor::BufferDistributor(const Gadgetron::Core::Context &context,
                                                                   const Gadgetron::Core::GadgetProperties &props) : TypedDistributor(props), encoding_spaces{context.header.encoding.size()} {
}

namespace Gadgetron::Core::Distributed {
    GADGETRON_DISTRIBUTOR_EXPORT(BufferDistributor);
}
