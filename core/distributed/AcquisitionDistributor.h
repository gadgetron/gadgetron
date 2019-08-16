#pragma once

#include "Context.h"
#include "Distributor.h"
#include "Types.h"

namespace Gadgetron::Core::Distributed {
    class AcquisitionDistributor : public TypedDistributor<Core::Acquisition>  {

    public:
        AcquisitionDistributor(const Context& context, const GadgetProperties& props);

        void process(InputChannel<Acquisition> &input, ChannelCreator &creator) override;

        NODE_PROPERTY(parallel_dimension, std::string, "Dimension that data will be parallelized over", "slice");

    private:
        const std::function<uint16_t(const ISMRMRD::AcquisitionHeader&)> &selector;
    };
}



