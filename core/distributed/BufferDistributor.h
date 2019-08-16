#pragma once

#include "Context.h"
#include "Distributor.h"

#include "mri_core_data.h"

namespace Gadgetron::Core::Distributed {
    class BufferDistributor : public TypedDistributor<IsmrmrdReconData> {

    public:
        BufferDistributor(const Context &context, const GadgetProperties &props);

        void process(InputChannel<IsmrmrdReconData> &input, ChannelCreator &creator) override;

    private:
        const size_t encoding_spaces;
    };
}


