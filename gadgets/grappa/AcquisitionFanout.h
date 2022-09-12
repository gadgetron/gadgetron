#pragma once

#include <ismrmrd/ismrmrd.h>

#include "SliceAccumulator.h"

#include "parallel/Fanout.h"
#include "Types.h"

namespace Gadgetron::Grappa {
    using AcquisitionFanout = Core::Parallel::Fanout<Slice>;
}