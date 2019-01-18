#include "AcquisitionFanout.h"

#include "parallel/Branch.h"

namespace Gadgetron::Grappa {

    AcquisitionFanout::AcquisitionFanout(
            const Core::Context &context,
            const Core::GadgetProperties &props
    ) : Fanout(context, props) {}

    GADGETRON_BRANCH_EXPORT(AcquisitionFanout);
}
