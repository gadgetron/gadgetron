#ifndef GADGETRON_PYTHON_VECTOR_CONVERTER_H
#define GADGETRON_PYTHON_VECTOR_CONVERTER_H

#include "python_toolbox.h"
#include "python_numpy_wrappers.h"
#include "log.h"

#include <vector>

#include <boost/python.hpp>
namespace bp = boost::python;

namespace Gadgetron {

template <typename T>
struct vector_to_numpy_array {
    static PyObject* convert(const std::vector<T>& vec)
    {
        std::vector<npy_intp> dims(1);
        dims[0] = vec.size();

        // TODO: This probably only works for types that map to NumPy types
        // so e.g. a std::vector<std::string> shouldn't work
        PyObject* obj = NumPyArray_SimpleNew(dims.size(), &dims[0], get_numpy_type<T>());
        if (sizeof(T) != NumPyArray_ITEMSIZE(obj)) {
            GERROR("sizeof(T): %d, ITEMSIZE: %d\n", sizeof(T), NumPyArray_ITEMSIZE(obj));
            throw std::runtime_error("vector_to_numpy_array: "
                    "python object and std::vector data type sizes do not match");
        }

        // Copy data... this is safe right? or use a for-loop
        memcpy(NumPyArray_DATA(obj), &vec[0], vec.size() * sizeof(T));

        // increment the reference count so it exists after `return`
        return bp::incref(obj);
    }
};

template <typename T>
struct vector_from_numpy_array {
    vector_from_numpy_array() {
        // actually register this converter with Boost
        bp::converter::registry::push_back(
                &convertible,
                &construct,
                bp::type_id<std::vector<T> >());
    }

    static void* convertible(PyObject* obj) {
        if (sizeof(T) != NumPyArray_ITEMSIZE(obj)) {
            GERROR("sizeof(T): %d, ITEMSIZE: %d\n", sizeof(T), NumPyArray_ITEMSIZE(obj));
            return NULL;
        }
        return obj;
    }

    static void construct(PyObject* obj, bp::converter::rvalue_from_python_stage1_data* data) {
        void* storage = ((bp::converter::rvalue_from_python_storage<hoNDArray<T> >*)data)->storage.bytes;
        data->convertible = storage;

        size_t length = NumPyArray_SIZE(obj);
        std::vector<T>* vec = new (storage) std::vector<T>(length);
        memcpy(&(*vec)[0], NumPyArray_DATA(obj), sizeof(T) * length);
    }
};

/// Create and register vector converter as necessary
template <typename T> void create_vector_converter() {
    bp::type_info info = bp::type_id<std::vector<T> >();
    const bp::converter::registration* reg = bp::converter::registry::query(info);
    // only register if not already registered!
    if (nullptr == reg || nullptr == (*reg).m_to_python) {
        bp::to_python_converter<std::vector<T>, vector_to_numpy_array<T> >();
        vector_from_numpy_array<T>();
    }
}

/// Partial specialization of `python_converter` for std::vector
template <typename T>
struct python_converter<std::vector<T> > {
    static void create()
    {
        // ensure NumPy C-API is initialized
        initialize_numpy();
        // register std::vector converter
        create_vector_converter<T>();
    }
};

}

#endif /* GADGETRON_PYTHON_VECTOR_CONVERTER_H */
