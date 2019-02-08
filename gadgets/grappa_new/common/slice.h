#pragma once

#include "Types.h"

namespace Gadgetron::Grappa {
    bool is_last_in_slice(const Core::Acquisition &acquisition);
    size_t slice_of(const Core::Acquisition &acquisition);
    size_t line_of(const Core::Acquisition &acquisition);
    size_t samples_in(const Core::Acquisition &acquisition);
}