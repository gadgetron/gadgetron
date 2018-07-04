#ifndef GADGETRON_PYTHON_ISMRMRD_CONVERTER_H
#define GADGETRON_PYTHON_ISMRMRD_CONVERTER_H

#include "python_toolbox.h" // for pyerr_to_string()
#include "ismrmrd/ismrmrd.h"
#include "ismrmrd/meta.h"

#include <boost/python.hpp>
#include <ismrmrd/waveform.h>

namespace bp = boost::python;

namespace Gadgetron {

    namespace {
       template<class T> T* get_ctypes_ptr(bp::object& obj){
           static bp::object ctypes_addressof = bp::import("ctypes").attr("addressof");
           T* py_ptr = reinterpret_cast<T*>(uintptr_t(bp::extract<uintptr_t>(ctypes_addressof(obj))));
           return py_ptr;
       }

       template<class T> bp::object python_module();

       template<> bp::object python_module<ISMRMRD::AcquisitionHeader>(){
           static bp::object module = bp::import("ismrmrd");
           return module.attr("AcquisitionHeader");
       }
      template<> bp::object python_module<ISMRMRD::ImageHeader>(){
           static bp::object module = bp::import("ismrmrd");
           return module.attr("ImageHeader");
       }
              template<> bp::object python_module<ISMRMRD::WaveformHeader>(){
           static bp::object module = bp::import("ismrmrd");
           return module.attr("WaveformHeader");
       }
    }
// -------------------------------------------------------------------------------------------------------
// ISMRMRD::AcquisitionHeader

template<class T> struct Header_to_PythonHeader {
    static PyObject* convert(const T& head) {
        try {
            GILLock lock;
            bp::object module = bp::import("ismrmrd");
            bp::object header_class = python_module<T>();

            bp::object pyhead = header_class();
            auto py_ptr = get_ctypes_ptr<T>(pyhead);
             *py_ptr = head;
            // increment the reference count so it exists after `return`
            return bp::incref(pyhead.ptr());
        } catch (const bp::error_already_set&) {
            std::string err = pyerr_to_string();
            GERROR(err.c_str());
            throw std::runtime_error(err);
        }
    }
};

template<class T> struct Header_from_PythonHeader {
    Header_from_PythonHeader() {
        // actually register this converter with Boost
        bp::converter::registry::push_back(
                &convertible,
                &construct,
                bp::type_id<T>());
    }

    /// Returns NULL if the Python AcquisitionHeader is not convertible
    static void* convertible(PyObject* obj) {
        return obj;
    }

    /// Construct an ISMRMRD::AcquisitionHeader in-place
    static void construct(PyObject* obj, bp::converter::rvalue_from_python_stage1_data* data) {
        void* storage = ((bp::converter::rvalue_from_python_storage<T>*)data)->storage.bytes;

        // Placement-new of ISMRMRD::AcquisitionHeader in memory provided by Boost
        auto head = new (storage) T;
        data->convertible = storage;

        try {
            bp::object pyhead((bp::handle<>(bp::borrowed(obj))));

            auto py_ptr = get_ctypes_ptr<T>(pyhead);
            *head = *py_ptr;
        } catch (const bp::error_already_set&) {
            std::string err = pyerr_to_string();
            GERROR(err.c_str());
            throw std::runtime_error(err);
        }
    }
};

// ---------------------
struct MetaContainer_to_PythonMetaContainer
{
    static PyObject* convert(const ISMRMRD::MetaContainer& meta)
    {
        try
        {
            bp::object module = bp::import("ismrmrd");
            bp::object pymeta = module.attr("Meta")();

            std::stringstream str;
            ISMRMRD::serialize(const_cast<ISMRMRD::MetaContainer&>(meta), str);

            pymeta = boost::python::object(str.str());

            // increment the reference count so it exists after `return`
            return bp::incref(pymeta.ptr());
        }
        catch (const bp::error_already_set&)
        {
            std::string err = pyerr_to_string();
            GERROR(err.c_str());
            throw std::runtime_error(err);
        }
    }
};

struct MetaContainer_from_PythonMetaContainer
{
    MetaContainer_from_PythonMetaContainer()
    {
        bp::converter::registry::push_back(
            &convertible,
            &construct,
            bp::type_id<ISMRMRD::MetaContainer>());
    }

    static void* convertible(PyObject* obj)
    {
        return obj;
    }

    /// Construct an ISMRMRD::MetaContainer in-place
    static void construct(PyObject* obj, bp::converter::rvalue_from_python_stage1_data* data)
    {
        void* storage = ((bp::converter::rvalue_from_python_storage<ISMRMRD::MetaContainer>*)data)->storage.bytes;

        ISMRMRD::MetaContainer* meta = new (storage) ISMRMRD::MetaContainer;
        data->convertible = storage;

        try {
            bp::object pyMeta((bp::handle<>(bp::borrowed(obj))));
            std::string meta_str = bp::extract<std::string>(pyMeta);
            ISMRMRD::deserialize(meta_str.c_str(), *meta);
        }
        catch (const bp::error_already_set&)
        {
            std::string err = pyerr_to_string();
            GERROR(err.c_str());
            throw std::runtime_error(err);
        }
    }
};

// -------------------------------------------------------------------------------------------------------
/// Create and register WaveformHeader converter as necessary
inline void create_ismrmrd_WaveformHeader_converter() {
    bp::type_info info = bp::type_id<ISMRMRD::WaveformHeader>();
    const bp::converter::registration* reg = bp::converter::registry::query(info);
    // only register if not already registered!
    if (nullptr == reg || nullptr == (*reg).m_to_python) {
        bp::to_python_converter<ISMRMRD::WaveformHeader,Header_to_PythonHeader<ISMRMRD::WaveformHeader>>();
        Header_from_PythonHeader<ISMRMRD::WaveformHeader>();
    }
}
/// Create and register AcquisitionHeader converter as necessary
inline void create_ismrmrd_AcquisitionHeader_converter() {
    bp::type_info info = bp::type_id<ISMRMRD::AcquisitionHeader>();
    const bp::converter::registration* reg = bp::converter::registry::query(info);
    // only register if not already registered!
    if (nullptr == reg || nullptr == (*reg).m_to_python) {
        bp::to_python_converter<ISMRMRD::AcquisitionHeader, Header_to_PythonHeader<ISMRMRD::AcquisitionHeader>>();
        Header_from_PythonHeader<ISMRMRD::AcquisitionHeader>();
    }
}

// -------------------------------------------------------------------------------------------------------
/// Create and register ImageHeader converter as necessary
inline void create_ismrmrd_ImageHeader_converter() {
    bp::type_info info = bp::type_id<ISMRMRD::ImageHeader>();
    const bp::converter::registration* reg = bp::converter::registry::query(info);
    // only register if not already registered!
    if (nullptr == reg || nullptr == (*reg).m_to_python) {
        bp::to_python_converter<ISMRMRD::ImageHeader, Header_to_PythonHeader<ISMRMRD::ImageHeader>>();
        Header_from_PythonHeader<ISMRMRD::ImageHeader>();
    }
}

// -------------------------------------------------------------------------------------------------------
/// Create and register MetaContainer converter as necessary
inline void create_ismrmrd_MetaContainer_converter()
{
    bp::type_info info = bp::type_id<ISMRMRD::MetaContainer>();
    const bp::converter::registration* reg = bp::converter::registry::query(info);
    // only register if not already registered!
    if (nullptr == reg || nullptr == (*reg).m_to_python)
    {
        bp::to_python_converter<ISMRMRD::MetaContainer, MetaContainer_to_PythonMetaContainer>();
        MetaContainer_from_PythonMetaContainer();
    }
}

// -------------------------------------------------------------------------------------------------------
/// Partial specialization of `python_converter` for ISMRMRD::AcquisitionHeader
template<> struct python_converter<ISMRMRD::WaveformHeader> {
    static void create()
    {
        create_ismrmrd_WaveformHeader_converter();
    }
};
// -------------------------------------------------------------------------------------------------------
/// Partial specialization of `python_converter` for ISMRMRD::AcquisitionHeader
template<> struct python_converter<ISMRMRD::AcquisitionHeader> {
    static void create()
    {
        create_ismrmrd_AcquisitionHeader_converter();
    }
};

// -------------------------------------------------------------------------------------------------------
/// Partial specialization of `python_converter` for ISMRMRD::ImageHeader
template<> struct python_converter<ISMRMRD::ImageHeader> {
    static void create()
    {
        create_ismrmrd_ImageHeader_converter();
    }
};

// -------------------------------------------------------------------------------------------------------
/// Partial specialization of `python_converter` for ISMRMRD::MetaContainer
template<> struct python_converter<ISMRMRD::MetaContainer> {
    static void create()
    {
        create_ismrmrd_MetaContainer_converter();
    }
};

}

#endif /* GADGETRON_PYTHON_ISMRMRD_CONVERTER_H */
