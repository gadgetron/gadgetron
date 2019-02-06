
#include "slice.h"

#include <ismrmrd/ismrmrd.h>

#include "Types.h"

namespace Gadgetron::Grappa {
    using namespace Gadgetron::Core;

    bool is_last_in_slice(const Acquisition &acquisition) {
        return std::get<ISMRMRD::AcquisitionHeader>(acquisition).isFlagSet(ISMRMRD::ISMRMRD_ACQ_LAST_IN_SLICE);
    }

    size_t slice_of(const Core::Acquisition &acquisition) {
        return std::get<ISMRMRD::AcquisitionHeader>(acquisition).idx.slice;
    }

    size_t line_of(const Core::Acquisition &acquisition) {
        return std::get<ISMRMRD::AcquisitionHeader>(acquisition).idx.kspace_encode_step_1;
    }
}