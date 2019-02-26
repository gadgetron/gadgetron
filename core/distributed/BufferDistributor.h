#pragma once

#include <Context.h>
#include "Distributor.h"
#include "../../toolboxes/mri_core/mri_core_data.h"

namespace Gadgetron::Core::Distributed {
    class BufferDistributor : public TypedDistributor<IsmrmrdReconData> {

    public:
        BufferDistributor(const Context &context, const GadgetProperties &props);

        void process(TypedInputChannel<IsmrmrdReconData> &input, ChannelCreator &creator) override;

    private:
        const size_t encoding_spaces;
    };


}


