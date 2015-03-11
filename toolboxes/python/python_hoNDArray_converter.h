#ifndef GADGETRON_PYTHON_HONDARRAY_CONVERTER_H
#define GADGETRON_PYTHON_HONDARRAY_CONVERTER_H

#include "python_toolbox.h"
#include "python_numpy_wrappers.h"

#include "hoNDArray.h"
#include "log.h"

#include <boost/python.hpp>
namespace bp = boost::python;

namespace Gadgetron {

/// return the enumerated numpy type for a given C++ type
template <typename T> int get_numpy_type();
template <> inline int get_numpy_type< char >() { return NPY_INT8; }
template <> inline int get_numpy_type< unsigned char >() { return NPY_UINT8; }
template <> inline int get_numpy_type< short >() { return NPY_INT16; }
template <> inline int get_numpy_type< unsigned short >() { return NPY_UINT16; }
template <> inline int get_numpy_type< int >() { return NPY_INT32; }
template <> inline int get_numpy_type< unsigned int >() { return NPY_UINT32; }
template <> inline int get_numpy_type< long >() { return NPY_INT64; }
template <> inline int get_numpy_type< unsigned long >() { return NPY_UINT64; }
template <> inline int get_numpy_type< float >() { return NPY_FLOAT32; }
template <> inline int get_numpy_type< double >() { return NPY_FLOAT64; }
template <> inline int get_numpy_type< std::complex<float> >() { return NPY_COMPLEX64; }
template <> inline int get_numpy_type< std::complex<double> >() { return NPY_COMPLEX128; }

/// Used for making a NumPy array from and hoNDArray
template <typename T>
struct hoNDArray_to_numpy_array {
    static PyObject* convert(const hoNDArray<T>& arr) {
        size_t ndim = arr.get_number_of_dimensions();
        std::vector<npy_intp> dims2(ndim);
        for (size_t i = 0; i < ndim; i++) {
            dims2[ndim-i-1] = static_cast<npy_intp>(arr.get_size(i));
        }
        PyObject *obj = NumPyArray_SimpleNew(dims2.size(), &dims2[0], get_numpy_type<T>());
        if (sizeof(T) != NumPyArray_ITEMSIZE(obj)) {
            GERROR("sizeof(T): %d, ITEMSIZE: %d\n", sizeof(T), NumPyArray_ITEMSIZE(obj));
            throw std::runtime_error("hondarray_to_bp_object: "
                    "python object and array data type sizes do not match");
        }

        // Copy data
        memcpy(NumPyArray_DATA(obj), arr.get_data_ptr(),
                arr.get_number_of_elements() * sizeof(T));

        // increment the reference count so it exists after `return`
        return bp::incref(obj);
    }
};

/// Used for making an hoNDArray from a NumPy array
template <typename T>
struct hoNDArray_from_numpy_array {
    hoNDArray_from_numpy_array() {
        // actually register this converter with Boost
        bp::converter::registry::push_back(
                &convertible,
                &construct,
                bp::type_id<hoNDArray<T> >());
    }

    /// Returns NULL if the NumPy array is not convertible
    static void* convertible(PyObject* obj) {
        if (sizeof(T) != NumPyArray_ITEMSIZE(obj)) {
            GERROR("sizeof(T): %d, ITEMSIZE: %d\n", sizeof(T), NumPyArray_ITEMSIZE(obj));
            return NULL;
        }
        return obj;
    }

    /// Construct an hoNDArray in-place
    static void construct(PyObject* obj, bp::converter::rvalue_from_python_stage1_data* data) {
        void* storage = ((bp::converter::rvalue_from_python_storage<hoNDArray<T> >*)data)->storage.bytes;
        size_t ndim = NumPyArray_NDIM(obj);
        std::vector<size_t> dims(ndim);
        for (size_t i = 0; i < ndim; i++) {
            dims[ndim - i - 1] = NumPyArray_DIM(obj, i);
        }

        // Placement-new of hoNDArray in memory provided by Boost
        hoNDArray<T>* arr = new (storage) hoNDArray<T>(dims);
        memcpy(arr->get_data_ptr(), NumPyArray_DATA(obj),
                sizeof(T) * arr->get_number_of_elements());
        data->convertible = storage;
    }
};

/// Create and register hoNDArray converter as necessary
template <typename T> void create_hoNDArray_converter() {
    bp::type_info info = bp::type_id<hoNDArray<T> >();
    const bp::converter::registration* reg = bp::converter::registry::query(info);
    // only register if not already registered!
    if (nullptr == reg || nullptr == (*reg).m_to_python) {
        bp::to_python_converter<hoNDArray<T>, hoNDArray_to_numpy_array<T> >();
        hoNDArray_from_numpy_array<T>();
    }
}

/// Partial specialization of `python_converter` for hoNDArray
template <typename T>
struct python_converter<hoNDArray<T> > {
    static void create()
    {
        // ensure NumPy C-API is initialized
        initialize_numpy();
        // register hoNDArray converter
        create_hoNDArray_converter<T>();
    }
};

}

#endif // GADGETRON_PYTHON_HONDARRAY_CONVERTER_H
