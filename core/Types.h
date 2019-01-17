#pragma once

#include "ismrmrd/ismrmrd.h"
#include "hoNDArray.h"
#include <boost/optional.hpp>
#include <boost/variant.hpp>
#include <tuple>
namespace Gadgetron::Core {

    template<class T>
    using optional = boost::optional<T>;

    template<class... ARGS>
    using variant = boost::variant<ARGS...>;

    template<class... ARGS>
    using tuple = std::tuple<ARGS...>;

    using Acquisition = tuple<ISMRMRD::AcquisitionHeader,optional<hoNDArray<float>>, hoNDArray<std::complex<float>>>;

}


#include "Types.hpp"