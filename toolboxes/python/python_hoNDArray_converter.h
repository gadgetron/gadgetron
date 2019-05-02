#ifndef GADGETRON_PYTHON_HONDARRAY_CONVERTER_H
#define GADGETRON_PYTHON_HONDARRAY_CONVERTER_H

#include "python_toolbox.h"
#include "python_numpy_wrappers.h"
#if (NPY_API_VERSION <= 6) && !defined(NPY_ARRAY_IN_FARRAY)
  // work-around for NumPy 1.6 (or earlier?)
  #define NPY_ARRAY_IN_FARRAY NPY_IN_FARRAY
#endif

#include "hoNDArray.h"
#include "log.h"

#include <boost/python.hpp>
namespace bp = boost::python;

namespace Gadgetron {

// -------------------------------------------------------------------------------
/// Used for making a NumPy array from and hoNDArray
template <typename T>
struct hoNDArray_to_numpy_array {
    static PyObject* convert(const hoNDArray<T>& arr) {
        size_t ndim = arr.get_number_of_dimensions();
        std::vector<npy_intp> dims2(ndim);
        for (size_t i = 0; i < ndim; i++) {
            dims2[i] = static_cast<npy_intp>(arr.get_size(i));
        }
        PyObject *obj = NumPyArray_EMPTY(dims2.size(), dims2.data(), get_numpy_type<T>(),true);
        if (sizeof(T) != NumPyArray_ITEMSIZE(obj)) {
            GERROR("sizeof(T): %d, ITEMSIZE: %d\n", sizeof(T), NumPyArray_ITEMSIZE(obj));
            throw std::runtime_error("hondarray_to_numpy_array: "
                    "python object and array data type sizes do not match");
        }

        // Copy data
        memcpy(NumPyArray_DATA(obj), arr.get_data_ptr(),
                arr.get_number_of_elements() * sizeof(T));

        // increment the reference count so it exists after `return`
        return obj;
    }
};

// -------------------------------------------------------------------------------
/// ISMRMRD::AcquisitionHeader
template <>
struct hoNDArray_to_numpy_array<ISMRMRD::AcquisitionHeader> {
    static PyObject* convert(const hoNDArray<ISMRMRD::AcquisitionHeader>& arr) {
        size_t ndim = arr.get_number_of_dimensions();
        std::vector<npy_intp> dims2(ndim);
        for (size_t i = 0; i < ndim; i++) {
            dims2[i] = static_cast<npy_intp>(arr.get_size(i));
        }
        PyObject *obj = NumPyArray_EMPTY(dims2.size(), dims2.data(), get_numpy_type<ISMRMRD::AcquisitionHeader>(),1);

        std::vector<PyObject*> pyobjects;
        for (auto & acq : arr){
          pyobjects.push_back(bp::incref(bp::object(acq).ptr()));
        }

        // Copy data
        memcpy(NumPyArray_DATA(obj), pyobjects.data(),
                pyobjects.size()* sizeof(PyObject*));

        // increment the reference count so it exists after `return`
        return obj;
    }
};

// -------------------------------------------------------------------------------
/// ISMRMRD::ImageHeader
template <>
struct hoNDArray_to_numpy_array<ISMRMRD::ImageHeader>
{
    static PyObject* convert(const hoNDArray<ISMRMRD::ImageHeader>& arr)
    {
        size_t ndim = arr.get_number_of_dimensions();
        std::vector<npy_intp> dims2(ndim);
        for (size_t i = 0; i < ndim; i++)
        {
            dims2[i] = static_cast<npy_intp>(arr.get_size(i));
        }
        PyObject *obj = NumPyArray_EMPTY(dims2.size(), dims2.data(), get_numpy_type<ISMRMRD::ImageHeader>(), 1);

        std::vector<PyObject*> pyobjects;
        for (auto & acq : arr)
        {
            pyobjects.push_back(bp::incref(bp::object(acq).ptr()));
        }

        // Copy data
        memcpy(NumPyArray_DATA(obj), pyobjects.data(),
            pyobjects.size() * sizeof(PyObject*));

        // increment the reference count so it exists after `return`
        return obj;
    }
};

// ===========================================================================================================

// -------------------------------------------------------------------------------
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
    static void construct(PyObject* obj_orig, bp::converter::rvalue_from_python_stage1_data* data) {
        void* storage = ((bp::converter::rvalue_from_python_storage<hoNDArray<T> >*)data)->storage.bytes;
        data->convertible = storage;

        PyObject* obj =  NumPyArray_FromAny(obj_orig, nullptr, 1, 36,  NPY_ARRAY_IN_FARRAY, nullptr);
        size_t ndim = NumPyArray_NDIM(obj);
        std::vector<size_t> dims(ndim);
        for (size_t i = 0; i < ndim; i++) {
            dims[i] = NumPyArray_DIM(obj, i);
        }
        // Placement-new of hoNDArray in memory provided by Boost
        hoNDArray<T>* arr = new (storage) hoNDArray<T>(dims);
        memcpy(arr->get_data_ptr(), NumPyArray_DATA(obj),
                sizeof(T) * arr->get_number_of_elements());
        bp::decref(obj);
    }
};

// --------------------------------------------------------------------------------
/// Used for making an hoNDArray from a NumPy array
template <>
struct hoNDArray_from_numpy_array<ISMRMRD::AcquisitionHeader> {
    hoNDArray_from_numpy_array() {
        // actually register this converter with Boost
        bp::converter::registry::push_back(
                &convertible,
                &construct,
                bp::type_id<hoNDArray<ISMRMRD::AcquisitionHeader> >());
    }

    /// Returns NULL if the NumPy array is not convertible.
    //Or well.. it should
    static void* convertible(PyObject* obj) {
        return obj;
    }

    /// Construct an hoNDArray in-place
    static void construct(PyObject* obj_orig, bp::converter::rvalue_from_python_stage1_data* data) {
        void* storage = ((bp::converter::rvalue_from_python_storage<hoNDArray<ISMRMRD::AcquisitionHeader> >*)data)->storage.bytes;
        data->convertible = storage;
        //Ensure fortran byte order
        PyObject* obj =  NumPyArray_FromAny(obj_orig, nullptr, 1, 36,  NPY_ARRAY_IN_FARRAY, nullptr);
        size_t ndim = NumPyArray_NDIM(obj);
        std::vector<size_t> dims(ndim);
        for (size_t i = 0; i < ndim; i++) {
            dims[i] = NumPyArray_DIM(obj, i);
        }

        // Placement-new of hoNDArray in memory provided by Boost
        hoNDArray<ISMRMRD::AcquisitionHeader>* arr = new (storage) hoNDArray<ISMRMRD::AcquisitionHeader>(dims);

        size_t elements = arr->get_number_of_elements();
        auto data_ptr = arr->get_data_ptr();
        PyObject** pyobjects = (PyObject**) NumPyArray_DATA(obj);

        for (size_t i = 0; i < elements; i++){
          data_ptr[i] = bp::extract<ISMRMRD::AcquisitionHeader>(bp::object(bp::borrowed(pyobjects[i])));
        }

        bp::decref(obj);

    }
};

// --------------------------------------------------------------------------------
template <>
struct hoNDArray_from_numpy_array<ISMRMRD::ImageHeader>
{
    hoNDArray_from_numpy_array()
    {
        bp::converter::registry::push_back(
            &convertible,
            &construct,
            bp::type_id<hoNDArray<ISMRMRD::ImageHeader> >());
    }

    static void* convertible(PyObject* obj) {
        return obj;
    }

    static void construct(PyObject* obj_orig, bp::converter::rvalue_from_python_stage1_data* data)
    {
        void* storage = ((bp::converter::rvalue_from_python_storage<hoNDArray<ISMRMRD::ImageHeader> >*)data)->storage.bytes;
        data->convertible = storage;
        //Ensure fortran byte order
        PyObject* obj = NumPyArray_FromAny(obj_orig, nullptr, 1, 36, NPY_ARRAY_IN_FARRAY, nullptr);
        size_t ndim = NumPyArray_NDIM(obj);
        std::vector<size_t> dims(ndim);
        for (size_t i = 0; i < ndim; i++)
        {
            dims[i] = NumPyArray_DIM(obj, i);
        }

        hoNDArray<ISMRMRD::ImageHeader>* arr = new (storage) hoNDArray<ISMRMRD::ImageHeader>(dims);

        size_t elements = arr->get_number_of_elements();
        auto data_ptr = arr->get_data_ptr();
        PyObject** pyobjects = (PyObject**)NumPyArray_DATA(obj);

        for (size_t i = 0; i < elements; i++)
        {
            data_ptr[i] = bp::extract<ISMRMRD::ImageHeader>(bp::object(bp::borrowed(pyobjects[i])));
        }

        bp::decref(obj);
    }
};

// --------------------------------------------------------------------------------
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
