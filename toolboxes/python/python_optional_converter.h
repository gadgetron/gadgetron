#pragma once

#include "python_converters.h"

namespace Gadgetron {

namespace Python {

template <typename T>
struct optional_converter {
    static void* convertible(PyObject* obj) {
        return obj;
    }

    static void construct(PyObject* obj, bp::converter::rvalue_from_python_stage1_data* data) {
        void* storage = ((bp::converter::rvalue_from_python_storage<std::optional<T>>*)data)->storage.bytes;
        data->convertible = storage;

        if (obj == Py_None) {
            new (storage) std::optional<T>();
        } else {
            new (storage) std::optional<T>(bp::extract<T>(bp::object(bp::handle<>(bp::borrowed(obj)))));
        }
    }

    static PyObject* convert(const std::optional<T>& opt) {
        if (!opt) {
            Py_RETURN_NONE;
        }
        return bp::incref(bp::object(*opt).ptr());
    }
};

} // namespace Python


template <typename T>
struct python_converter<std::optional<T>> {
    static void create()
    {
        // Register converter for inner type
        register_converter<T>();

        // Register optional_converter
        register_with_boost<std::optional<T>, Python::optional_converter<T>>();
    }
};

} // namespace Gadgetron