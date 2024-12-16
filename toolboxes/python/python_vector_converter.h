#pragma once

#include "python_converters.h"
#include "log.h"

#include <vector>

namespace Gadgetron {

namespace Python {

template <typename T>
struct vector_converter {
    static void* convertible(PyObject* obj) {
        return obj;
    }

    static void construct(PyObject* obj, bp::converter::rvalue_from_python_stage1_data* data) {
        try
        {
            void* storage = ((bp::converter::rvalue_from_python_storage<std::vector<T> >*)data)->storage.bytes;
            data->convertible = storage;

            bp::list pyVec((bp::handle<>(bp::borrowed(obj))));
            auto length = bp::len(pyVec);

            std::vector<T>* vec = new (storage) std::vector<T>(length);
            for (size_t n = 0; n < length; n++) {
                bp::object item = pyVec[n];
                (*vec)[n] = bp::extract<T>(item);
            }
        }
        catch (const bp::error_already_set&)
        {
            std::string err = pyerr_to_string();
            GERROR(err.c_str());
            throw std::runtime_error(err);
        }
    }

    static PyObject* convert(const std::vector<T>& vec)
    {
        try
        {
            auto pyVec = bp::list();
            for (size_t n = 0; n < vec.size(); n++) {
                auto item = bp::object(vec[n]);
                pyVec.append(item);
            }
            // increment the reference count so it exists after `return`
            return bp::incref(pyVec.ptr());
        }
        catch (const bp::error_already_set&)
        {
            std::string err = pyerr_to_string();
            GERROR(err.c_str());
            throw std::runtime_error(err);
        }
    }

};

} // namespace Python

/// Partial specialization of `python_converter` for std::vector
template <typename T>
struct python_converter<std::vector<T> > {
    static void create()
    {
        register_converter<T>();

        register_with_boost<std::vector<T>, Python::vector_converter<T>>();
    }
};

} // namespace Gadgetron