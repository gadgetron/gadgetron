#ifndef GADGETRON_PYTHON_VECTOR_CONVERTER_H
#define GADGETRON_PYTHON_VECTOR_CONVERTER_H

#include "python_toolbox.h"
#include "python_numpy_wrappers.h"
#include "log.h"

#include "ismrmrd/ismrmrd.h"
#include "ismrmrd/meta.h"

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

// -----------------------------------------------------------------------------------

template <>
struct vector_to_numpy_array<std::string>
{
    static PyObject* convert(const std::vector<std::string>& vec)
    {
        try
        {
            std::vector<npy_intp> dims(1);
            dims[0] = vec.size();

            auto pyVecStr = bp::list();
            for (size_t n=0; n<vec.size(); n++)
            {
                auto curr_str = bp::object(vec[n]);
                pyVecStr.append(curr_str);

            }
            // increment the reference count so it exists after `return`
            return bp::incref(pyVecStr.ptr());
        }
        catch (const bp::error_already_set&)
        {
            std::string err = pyerr_to_string();
            GERROR(err.c_str());
            GERROR_STREAM("Exceptions happened in vector_to_numpy_array<std::string> ... ");
            throw std::runtime_error(err);
        }
    }
};

template <>
struct vector_from_numpy_array<std::string>
{
    vector_from_numpy_array()
    {
        bp::converter::registry::push_back(
            &convertible,
            &construct,
            bp::type_id<std::vector<std::string> >());
    }

    static void* convertible(PyObject* obj)
    {
        return obj;
    }

    static void construct(PyObject* obj, bp::converter::rvalue_from_python_stage1_data* data)
    {
        try
        {
            void* storage = ((bp::converter::rvalue_from_python_storage<std::vector<std::string> >*)data)->storage.bytes;
            data->convertible = storage;

            bp::list pyVecStr((bp::handle<>(bp::borrowed(obj))));
            auto length = bp::len(pyVecStr);

            std::vector<std::string>* vec = new (storage) std::vector<std::string>(length);
            for (size_t n=0; n<length; n++)
            {
                bp::object curr_str = pyVecStr[n];
                (*vec)[n] = bp::extract<std::string>(curr_str);
            }
        }
        catch (const bp::error_already_set&)
        {
            std::string err = pyerr_to_string();
            GERROR(err.c_str());
            GERROR_STREAM("Exceptions happened in vector_from_numpy_array<std::string> ... ");
            throw std::runtime_error(err);
        }
    }
};

// -----------------------------------------------------------------------------------

template <>
struct vector_to_numpy_array<ISMRMRD::MetaContainer>
{
    static PyObject* convert(const std::vector<ISMRMRD::MetaContainer>& vec)
    {
        try
        {
            std::vector<npy_intp> dims(1);
            dims[0] = vec.size();

            auto pyVecStr = bp::list();
            for (size_t n = 0; n<vec.size(); n++)
            {
                std::stringstream str;
                ISMRMRD::serialize( const_cast<ISMRMRD::MetaContainer&>(vec[n]), str);
                auto curr_str = bp::object(str.str());
                pyVecStr.append(curr_str);

            }
            // increment the reference count so it exists after `return`
            return bp::incref(pyVecStr.ptr());
        }
        catch (const bp::error_already_set&)
        {
            std::string err = pyerr_to_string();
            GERROR(err.c_str());
            GERROR_STREAM("Exceptions happened in vector_to_numpy_array<ISMRMRD::MetaContainer> ... ");
            throw std::runtime_error(err);
        }
    }
};

template <>
struct vector_from_numpy_array<ISMRMRD::MetaContainer>
{
    vector_from_numpy_array()
    {
        bp::converter::registry::push_back(
            &convertible,
            &construct,
            bp::type_id<std::vector<ISMRMRD::MetaContainer> >());
    }

    static void* convertible(PyObject* obj)
    {
        return obj;
    }

    static void construct(PyObject* obj, bp::converter::rvalue_from_python_stage1_data* data)
    {
        try
        {
            void* storage = ((bp::converter::rvalue_from_python_storage<std::vector<ISMRMRD::MetaContainer> >*)data)->storage.bytes;
            data->convertible = storage;

            bp::list pyVecStr((bp::handle<>(bp::borrowed(obj))));
            auto length = bp::len(pyVecStr);

            std::vector<ISMRMRD::MetaContainer>* vec = new (storage) std::vector<ISMRMRD::MetaContainer>(length);
            for (size_t n = 0; n<length; n++)
            {
                bp::object curr_str = pyVecStr[n];
                std::string str = bp::extract<std::string>(curr_str);
                ISMRMRD::deserialize(str.c_str(), (*vec)[n]);
            }
        }
        catch (const bp::error_already_set&)
        {
            std::string err = pyerr_to_string();
            GERROR(err.c_str());
            GERROR_STREAM("Exceptions happened in vector_from_numpy_array<ISMRMRD::MetaContainer> ... ");
            throw std::runtime_error(err);
        }
    }
};

// -----------------------------------------------------------------------------------

template <>
struct vector_to_numpy_array<ISMRMRD::Waveform>
{
    static PyObject* convert(const std::vector<ISMRMRD::Waveform>& vec)
    {
        try
        {
            std::vector<npy_intp> dims(1);
            dims[0] = vec.size();

            auto pyVecWav = bp::list();
            for (size_t n = 0; n<vec.size(); n++)
            {
                auto curr_wav = bp::object(vec[n]);
                bp::incref(curr_wav.ptr());
                pyVecWav.append(curr_wav);

            }
            // increment the reference count so it exists after `return`
            return bp::incref(pyVecWav.ptr());
        }
        catch (const bp::error_already_set&)
        {
            std::string err = pyerr_to_string();
            GERROR(err.c_str());
            GERROR_STREAM("Exceptions happened in vector_to_numpy_array<ISMRMRD::Waveform> ... ");
            throw std::runtime_error(err);
        }
    }
};

template <>
struct vector_from_numpy_array<ISMRMRD::Waveform>
{
    vector_from_numpy_array()
    {
        bp::converter::registry::push_back(
            &convertible,
            &construct,
            bp::type_id<std::vector<ISMRMRD::Waveform> >());
    }

    static void* convertible(PyObject* obj)
    {
        return obj;
    }

    static void construct(PyObject* obj, bp::converter::rvalue_from_python_stage1_data* data)
    {
        try
        {
            void* storage = ((bp::converter::rvalue_from_python_storage<std::vector<ISMRMRD::Waveform> >*)data)->storage.bytes;
            data->convertible = storage;

            bp::list pyVecWav((bp::handle<>(bp::borrowed(obj))));
            auto length = bp::len(pyVecWav);

            std::vector<ISMRMRD::Waveform>* vec = new (storage) std::vector<ISMRMRD::Waveform>(length);
            for (size_t n = 0; n<length; n++)
            {
                bp::object curr_wav = pyVecWav[n];
                ISMRMRD::Waveform wav = bp::extract<ISMRMRD::Waveform>(curr_wav);
                (*vec)[n] = wav;
            }
        }
        catch (const bp::error_already_set&)
        {
            std::string err = pyerr_to_string();
            GERROR(err.c_str());
            GERROR_STREAM("Exceptions happened in vector_from_numpy_array<ISMRMRD::Waveform> ... ");
            throw std::runtime_error(err);
        }
    }
};

// -----------------------------------------------------------------------------------
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

template <>
struct python_converter<std::vector<std::string> >
{
    static void create()
    {
        create_vector_converter<std::string>();
    }
};

template <>
struct python_converter<std::vector<ISMRMRD::MetaContainer> >
{
    static void create()
    {
        create_vector_converter<ISMRMRD::MetaContainer>();
    }
};

template <>
struct python_converter<std::vector<ISMRMRD::Waveform> >
{
    static void create()
    {
        create_vector_converter<ISMRMRD::Waveform>();
    }
};

}

#endif /* GADGETRON_PYTHON_VECTOR_CONVERTER_H */
