#pragma once

#include "SliceAccumulator.h"

#include "parallel/Fanout.h"

namespace Gadgetron::Grappa {
    using AcquisitionFanout = Core::Parallel::Fanout<Slice>;
}