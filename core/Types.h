#pragma once

#include <ismrmrd/ismrmrd.h>
#include <ismrmrd/meta.h>
#include <ismrmrd/waveform.h>

#include <tuple>
#include "hoNDArray.h"
#include "TypeTraits.h"
#include <optional>
#include <variant>

namespace Gadgetron::Core {
    template<class T>
    using optional = std::optional<T>;
    constexpr auto none = std::nullopt;

    template<class... TYPES>
    using variant= std::variant<TYPES...>;
    using std::visit;
    using std::holds_alternative;
    using std::get;


    template<class... ARGS>
    using tuple = std::tuple<ARGS...>;

/// An Acquisition consists of a data header, the kspace data itself and optionally an array of kspace trajectories
    using Acquisition = tuple<ISMRMRD::AcquisitionHeader,  hoNDArray<std::complex<float>>,optional<hoNDArray<float>>>;

/// A Waveform consiste of a header, followed by the raw Waveform data. See the MRD documentation page for more details
    using Waveform    = tuple<ISMRMRD::WaveformHeader, hoNDArray<uint32_t>>;

    ///An image consists of a header, an array of image data and optionally some metadata
    template<class T>
    using Image = tuple<ISMRMRD::ImageHeader,  hoNDArray<T>, optional<ISMRMRD::MetaContainer>>;

    using AnyImage =
            variant<
                    Image<short>,
                    Image<unsigned short>,
                    Image<float>,
                    Image<double>,
                    Image<std::complex<float>>,
                    Image<std::complex<double>>
            >;
}


#include "Types.hpp"