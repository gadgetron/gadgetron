#pragma once

#include <boost/python.hpp>
namespace bp = boost::python;

namespace Gadgetron {

/// Interface for registering C++ <-> NumPy type converters.
/// A static function on the `python_converter` struct allows
/// for partial template specialization.
template <typename T>
struct python_converter {
    static void create() { }
};

/// User-facing interface for registering a Python converter for specific types
template <typename ...TS>
void register_converter(void) {
    // Register the converter for each template type
    (python_converter<TS>::create(), ...);
}

/// Internal interface for registering a C++ type with a Python converter
template <typename T, typename C>
void register_with_boost() {
    bp::type_info info = bp::type_id<T>();
    const bp::converter::registration* reg = bp::converter::registry::query(info);
    if (nullptr == reg || nullptr == reg->m_to_python) {
        bp::converter::registry::push_back(&C::convertible, &C::construct, bp::type_id<T>());
        bp::to_python_converter<T, C>();
    }
}

}

#include "patchlevel.h"
#include "python_tuple_converter.h"
#include "python_optional_converter.h"
#include "python_vector_converter.h"
#include "python_hoNDArray_converter.h"
#include "python_mrd_converters.h"