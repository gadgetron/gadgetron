#pragma once

#include <ismrmrd/ismrmrd.h>
#include <ismrmrd/waveform.h>

#include <boost/optional.hpp>
#include <boost/variant.hpp>
#include <tuple>

#include "hoNDArray.h"

namespace Gadgetron::Core {

    template<class T>
    using optional = boost::optional<T>;
    static const auto none = boost::none;

    template<class... ARGS>
    using variant = boost::variant<ARGS...>;

    template<class... ARGS>
    using tuple = std::tuple<ARGS...>;


    using Acquisition = tuple<ISMRMRD::AcquisitionHeader, optional<hoNDArray<float>>, hoNDArray<std::complex<float>>>;
    using Waveform    = tuple<ISMRMRD::WaveformHeader, hoNDArray<uint32_t>>;
}


#include "Types.hpp"