#pragma once

#include "complext.h"
#include <boost/type_traits.hpp>

namespace Gadgetron {

    template <class T> constexpr bool enable_operator(){ return false; };
    template <> constexpr bool enable_operator<float>(){ return true; };
    template <> constexpr bool enable_operator<double>(){ return true; };

    template <> constexpr bool enable_operator<complext<float>>(){ return true; };
    template <> constexpr bool enable_operator<complext<double>>(){ return true; };

    template<class T> constexpr bool enable_operator_v = enable_operator<T>();

}
