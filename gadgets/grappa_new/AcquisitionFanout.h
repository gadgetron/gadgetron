#pragma once

#include <ismrmrd/ismrmrd.h>

#include "parallel/Fanout.h"
#include "hoNDArray.h"
#include "Types.h"

namespace Gadgetron::Grappa {
//    using AcquisitionFanout = Core::Parallel::Fanout<Core::Acquisition>;

    class AcquisitionFanout : public Core::Parallel::Fanout<Core::Acquisition> {
    public:
        AcquisitionFanout(const Core::Context &, const Core::GadgetProperties &);
    };
}