//
// Created by dchansen on 1/16/19.
//

#pragma once

#include "Context.h"
#include "Distributor.h"
#include "Types.h"
#include <boost/hana.hpp>


namespace Gadgetron::Core::Distributed {
    class AcquisitionDistributor : public TypedDistributor<Core::Acquisition>  {

    public:

        AcquisitionDistributor(const Context& context, const GadgetProperties& props);

        void process(TypedInputChannel<Acquisition> &input, ChannelCreator &creator) override;



        NODE_PROPERTY(parallel_dimensions,std::string, "slice", "Dimension that data will be parallelized over" );

    private:
        const std::function<uint16_t(const ISMRMRD::AcquisitionHeader&)> selector;

    };
}



