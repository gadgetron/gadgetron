#ifndef GADGETRON_PYTHON_ISMRMRD_CONVERTER_H
#define GADGETRON_PYTHON_ISMRMRD_CONVERTER_H

#include "python_toolbox.h" // for pyerr_to_string()
#include "ismrmrd/ismrmrd.h"
#include "ismrmrd/meta.h"

#include <boost/python.hpp>
namespace bp = boost::python;

namespace Gadgetron {

    namespace {
       template<class T> T* get_ctypes_ptr(bp::object& obj){
           static bp::object ctypes_addressof = bp::import("ctypes").attr("addressof");
           T* py_ptr = reinterpret_cast<T*>(uintptr_t(bp::extract<uintptr_t>(ctypes_addressof(obj))));
           return py_ptr;
       }
    }
// -------------------------------------------------------------------------------------------------------
// ISMRMRD::AcquisitionHeader

struct AcquisitionHeader_to_PythonAcquisitionHeader {
    static PyObject* convert(const ISMRMRD::AcquisitionHeader& head) {
        try {
            GILLock lock;
            bp::object module = bp::import("ismrmrd");
            bp::object acqheader_class = module.attr("AcquisitionHeader");

            bp::object pyhead = acqheader_class();
            auto py_ptr = get_ctypes_ptr<ISMRMRD::AcquisitionHeader>(pyhead);
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

struct AcquisitionHeader_from_PythonAcquisitionHeader {
    AcquisitionHeader_from_PythonAcquisitionHeader() {
        // actually register this converter with Boost
        bp::converter::registry::push_back(
                &convertible,
                &construct,
                bp::type_id<ISMRMRD::AcquisitionHeader>());
    }

    /// Returns NULL if the Python AcquisitionHeader is not convertible
    static void* convertible(PyObject* obj) {
        return obj;
    }

    /// Construct an ISMRMRD::AcquisitionHeader in-place
    static void construct(PyObject* obj, bp::converter::rvalue_from_python_stage1_data* data) {
        void* storage = ((bp::converter::rvalue_from_python_storage<ISMRMRD::AcquisitionHeader>*)data)->storage.bytes;

        // Placement-new of ISMRMRD::AcquisitionHeader in memory provided by Boost
        ISMRMRD::AcquisitionHeader* head = new (storage) ISMRMRD::AcquisitionHeader;
        data->convertible = storage;

        try {
            bp::object pyhead((bp::handle<>(bp::borrowed(obj))));

            auto py_ptr = get_ctypes_ptr<ISMRMRD::AcquisitionHeader>(pyhead);
            *head = *py_ptr;
        } catch (const bp::error_already_set&) {
            std::string err = pyerr_to_string();
            GERROR(err.c_str());
            throw std::runtime_error(err);
        }
    }
};

// -------------------------------------------------------------------------------------------------------
// ISMRMRD::ImageHeader

struct ImageHeader_to_PythonImageHeader {
    static PyObject* convert(const ISMRMRD::ImageHeader& head) {
        try {
            GILLock lock;
            bp::object module = bp::import("ismrmrd");
            bp::object pyhead = module.attr("ImageHeader")();

            pyhead.attr("version") = head.version;
            pyhead.attr("data_type") = head.data_type;
            pyhead.attr("flags") = head.flags;
            pyhead.attr("measurement_uid") = head.measurement_uid;
            for (int i = 0; i < ISMRMRD::ISMRMRD_POSITION_LENGTH; i++) {
                pyhead.attr("matrix_size")[i] = head.matrix_size[i];
            }
            for (int i = 0; i < ISMRMRD::ISMRMRD_POSITION_LENGTH; i++) {
                pyhead.attr("field_of_view")[i] = head.field_of_view[i];
            }
            pyhead.attr("channels") = head.channels;
            for (int i = 0; i < ISMRMRD::ISMRMRD_POSITION_LENGTH; i++) {
                pyhead.attr("position")[i] = head.position[i];
            }
            for (int i = 0; i < ISMRMRD::ISMRMRD_DIRECTION_LENGTH; i++) {
                pyhead.attr("read_dir")[i] = head.read_dir[i];
            }
            for (int i = 0; i < ISMRMRD::ISMRMRD_DIRECTION_LENGTH; i++) {
                pyhead.attr("phase_dir")[i] = head.phase_dir[i];
            }
            for (int i = 0; i < ISMRMRD::ISMRMRD_DIRECTION_LENGTH; i++) {
                pyhead.attr("slice_dir")[i] = head.slice_dir[i];
            }
            for (int i = 0; i < ISMRMRD::ISMRMRD_POSITION_LENGTH; i++) {
                pyhead.attr("patient_table_position")[i] = head.patient_table_position[i];
            }
            pyhead.attr("average") = head.average;
            pyhead.attr("slice") = head.slice;
            pyhead.attr("contrast") = head.contrast;
            pyhead.attr("phase") = head.phase;
            pyhead.attr("repetition") = head.repetition;
            pyhead.attr("set") = head.set;
            pyhead.attr("acquisition_time_stamp") = head.acquisition_time_stamp;
            for (int i = 0; i < ISMRMRD::ISMRMRD_PHYS_STAMPS; i++) {
                pyhead.attr("physiology_time_stamp")[i] = head.physiology_time_stamp[i];
            }
            pyhead.attr("image_type") = head.image_type;
            pyhead.attr("image_index") = head.image_index;
            pyhead.attr("image_series_index") = head.image_series_index;
            for (int i = 0; i < ISMRMRD::ISMRMRD_USER_INTS; i++) {
                pyhead.attr("user_int")[i] = head.user_int[i];
            }
            for (int i = 0; i < ISMRMRD::ISMRMRD_USER_FLOATS; i++) {
                pyhead.attr("user_float")[i] = head.user_float[i];
            }
            pyhead.attr("attribute_string_len") = head.attribute_string_len;

            // increment the reference count so it exists after `return`
            return bp::incref(pyhead.ptr());
        } catch (const bp::error_already_set&) {
            std::string err = pyerr_to_string();
            GERROR(err.c_str());
            throw std::runtime_error(err);
        }
    }
};

struct ImageHeader_from_PythonImageHeader {
    ImageHeader_from_PythonImageHeader() {
        // actually register this converter with Boost
        bp::converter::registry::push_back(
                &convertible,
                &construct,
                bp::type_id<ISMRMRD::ImageHeader>());
    }

    /// Returns NULL if the Python ImageHeader is not convertible
    static void* convertible(PyObject* obj) {
        return obj;
    }

    /// Construct an ISMRMRD::ImageHeader in-place
    static void construct(PyObject* obj, bp::converter::rvalue_from_python_stage1_data* data) {
        void* storage = ((bp::converter::rvalue_from_python_storage<ISMRMRD::ImageHeader>*)data)->storage.bytes;

        // Placement-new of ISMRMRD::ImageHeader in memory provided by Boost
        ISMRMRD::ImageHeader* head = new (storage) ISMRMRD::ImageHeader;
        data->convertible = storage;

        try {
            bp::object pyhead((bp::handle<>(bp::borrowed(obj))));

            head->version = bp::extract<uint16_t>(pyhead.attr("version"));
            head->data_type = bp::extract<uint16_t>(pyhead.attr("data_type"));
            head->flags = bp::extract<uint64_t>(pyhead.attr("flags"));
            head->measurement_uid = bp::extract<uint32_t>(pyhead.attr("measurement_uid"));
            for (int i = 0; i < ISMRMRD::ISMRMRD_POSITION_LENGTH; i++) {
                head->matrix_size[i] = bp::extract<uint16_t>(pyhead.attr("matrix_size")[i]);
            }
            for (int i = 0; i < ISMRMRD::ISMRMRD_POSITION_LENGTH; i++) {
                head->field_of_view[i] = bp::extract<float>(pyhead.attr("field_of_view")[i]);
            }
            head->channels = bp::extract<uint16_t>(pyhead.attr("channels"));
            for (int i = 0; i < ISMRMRD::ISMRMRD_POSITION_LENGTH; i++) {
                head->position[i] = bp::extract<float>(pyhead.attr("position")[i]);
            }
            for (int i = 0; i < ISMRMRD::ISMRMRD_DIRECTION_LENGTH; i++) {
                head->read_dir[i] = bp::extract<float>(pyhead.attr("read_dir")[i]);
            }
            for (int i = 0; i < ISMRMRD::ISMRMRD_DIRECTION_LENGTH; i++) {
                head->phase_dir[i] = bp::extract<float>(pyhead.attr("phase_dir")[i]);
            }
            for (int i = 0; i < ISMRMRD::ISMRMRD_DIRECTION_LENGTH; i++) {
                head->slice_dir[i] = bp::extract<float>(pyhead.attr("slice_dir")[i]);
            }
            for (int i = 0; i < ISMRMRD::ISMRMRD_POSITION_LENGTH; i++) {
                head->patient_table_position[i] = bp::extract<float>(pyhead.attr("patient_table_position")[i]);
            }
            head->average = bp::extract<uint16_t>(pyhead.attr("average"));
            head->slice = bp::extract<uint16_t>(pyhead.attr("slice"));
            head->contrast = bp::extract<uint16_t>(pyhead.attr("contrast"));
            head->phase = bp::extract<uint16_t>(pyhead.attr("phase"));
            head->repetition = bp::extract<uint16_t>(pyhead.attr("repetition"));
            head->set = bp::extract<uint16_t>(pyhead.attr("set"));
            head->acquisition_time_stamp = bp::extract<uint32_t>(pyhead.attr("acquisition_time_stamp"));
            for (int i = 0; i < ISMRMRD::ISMRMRD_PHYS_STAMPS; i++) {
                head->physiology_time_stamp[i] = bp::extract<uint32_t>(pyhead.attr("physiology_time_stamp")[i]);
            }
            head->image_type = bp::extract<uint16_t>(pyhead.attr("image_type"));
            head->image_index = bp::extract<uint16_t>(pyhead.attr("image_index"));
            head->image_series_index = bp::extract<uint16_t>(pyhead.attr("image_series_index"));
            for (int i = 0; i < ISMRMRD::ISMRMRD_USER_INTS; i++) {
                head->user_int[i] = bp::extract<int32_t>(pyhead.attr("user_int")[i]);
            }
            for (int i = 0; i < ISMRMRD::ISMRMRD_USER_FLOATS; i++) {
                head->user_float[i] = bp::extract<float>(pyhead.attr("user_float")[i]);
            }
            head->attribute_string_len = bp::extract<uint32_t>(pyhead.attr("attribute_string_len"));
        } catch (const bp::error_already_set&) {
            std::string err = pyerr_to_string();
            GERROR(err.c_str());
            throw std::runtime_error(err);
        }
    }
};

// -------------------------------------------------------------------------------------------------------
// ISMRMRD::MetaContainer

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

/// Create and register AcquisitionHeader converter as necessary
inline void create_ismrmrd_AcquisitionHeader_converter() {
    bp::type_info info = bp::type_id<ISMRMRD::AcquisitionHeader>();
    const bp::converter::registration* reg = bp::converter::registry::query(info);
    // only register if not already registered!
    if (nullptr == reg || nullptr == (*reg).m_to_python) {
        bp::to_python_converter<ISMRMRD::AcquisitionHeader, AcquisitionHeader_to_PythonAcquisitionHeader>();
        AcquisitionHeader_from_PythonAcquisitionHeader();
    }
}

// -------------------------------------------------------------------------------------------------------
/// Create and register ImageHeader converter as necessary
inline void create_ismrmrd_ImageHeader_converter() {
    bp::type_info info = bp::type_id<ISMRMRD::ImageHeader>();
    const bp::converter::registration* reg = bp::converter::registry::query(info);
    // only register if not already registered!
    if (nullptr == reg || nullptr == (*reg).m_to_python) {
        bp::to_python_converter<ISMRMRD::ImageHeader, ImageHeader_to_PythonImageHeader>();
        ImageHeader_from_PythonImageHeader();
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
