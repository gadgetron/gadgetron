#pragma once

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


}