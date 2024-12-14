#pragma once

#include "python_converters.h"
#include "python_toolbox.h"
#include "python_numpy_wrappers.h"

#include "hoNDArray.h"
#include "log.h"

namespace Gadgetron {

namespace Python {

/// Used for converting a NumPy array to/from an hoNDArray
template <typename T>
struct hoNDArray_converter
{
    /// Returns NULL if the NumPy array is not convertible
    static void* convertible(PyObject* obj) {
        if (get_numpy_type<T>() == NPY_OBJECT) {
            return obj;
        } else if (sizeof(T) != NumPyArray_ITEMSIZE(obj)) {
            GERROR("sizeof(T): %d, ITEMSIZE: %d\n", sizeof(T), NumPyArray_ITEMSIZE(obj));
            return NULL;
        }
        return obj;
    }

    /// Construct an hoNDArray in-place
    static void construct(PyObject* obj, bp::converter::rvalue_from_python_stage1_data* data) {
        void* storage = ((bp::converter::rvalue_from_python_storage<hoNDArray<T> >*)data)->storage.bytes;
        data->convertible = storage;

        size_t ndim = NumPyArray_NDIM(obj);
        std::vector<size_t> dims(ndim);
        for (size_t i = 0; i < ndim; i++) {
            dims[i] = NumPyArray_DIM(obj, ndim - i - 1);
        }

        // Placement-new of hoNDArray in memory provided by Boost
        hoNDArray<T>* arr = new (storage) hoNDArray<T>(dims);

        auto dtype = get_numpy_type<T>();
        if (dtype == NPY_OBJECT) {
            PyObject** pyobjects = (PyObject**) NumPyArray_DATA(obj);
            std::transform(pyobjects, pyobjects + NumPyArray_SIZE(obj), arr->begin(),
                    [](PyObject* obj) -> T { return bp::extract<T>(bp::object(bp::borrowed(obj))); });

        } else {
            memcpy(arr->get_data_ptr(), NumPyArray_DATA(obj), sizeof(T) * arr->get_number_of_elements());
        }
    }

    static PyObject* convert(const hoNDArray<T>& arr) {
        size_t ndim = arr.get_number_of_dimensions();
        std::vector<npy_intp> dims2(ndim);
        for (size_t i = 0; i < ndim; i++) {
            dims2[i] = static_cast<npy_intp>(arr.get_size(ndim - i - 1));
        }

        auto dtype = get_numpy_type<T>();
        PyObject *obj = NumPyArray_EMPTY(dims2.size(), dims2.data(), dtype, false);

        if (dtype == NPY_OBJECT) {
            std::transform(
                arr.begin(), arr.end(), (PyObject**)NumPyArray_DATA(obj),
                [](T const& item) -> PyObject* { return bp::incref(bp::object(item).ptr()); }
            );
        } else if (sizeof(T) == NumPyArray_ITEMSIZE(obj)) {
            memcpy(NumPyArray_DATA(obj), arr.get_data_ptr(),
                    arr.get_number_of_elements() * sizeof(T));
        } else {
            GERROR("sizeof(T): %d, ITEMSIZE: %d\n", sizeof(T), NumPyArray_ITEMSIZE(obj));
            throw std::runtime_error("hondarray_to_numpy_array: "
                    "python object and array data type sizes do not match");
        }

        return obj;
    }
};

} // namespace Python


/// Partial specialization of `python_converter` for hoNDArray
template <typename T>
struct python_converter<hoNDArray<T> > {
    static void create()
    {
        // Ensure NumPy C-API is initialized
        initialize_numpy();

        // Register converter for inner type
        register_converter<T>();

        // Register hoNDArray converter
        register_with_boost<hoNDArray<T>, Python::hoNDArray_converter<T>>();
    }
};

} // namespace Gadgetron