#ifndef GADGETRON_PYTHON_ISMRMRD_CONVERTER_H
#define GADGETRON_PYTHON_ISMRMRD_CONVERTER_H

#include "python_toolbox.h" // for pyerr_to_string()
#include "ismrmrd/ismrmrd.h"

#include <boost/python.hpp>
namespace bp = boost::python;

namespace Gadgetron {

// -------------------------------------------------------------------------------------------------------
// ISMRMRD::AcquisitionHeader

struct AcquisitionHeader_to_PythonAcquisitionHeader {
    static PyObject* convert(const ISMRMRD::AcquisitionHeader& head) {
        try {
            bp::object module = bp::import("ismrmrd");
            bp::object pyhead = module.attr("AcquisitionHeader")();

            pyhead.attr("version") = head.version;
            pyhead.attr("flags") = head.flags;
            pyhead.attr("measurement_uid") = head.measurement_uid;
            pyhead.attr("scan_counter") = head.scan_counter;
            pyhead.attr("acquisition_time_stamp") = head.acquisition_time_stamp;
            for (int i = 0; i < ISMRMRD::ISMRMRD_PHYS_STAMPS; i++) {
                pyhead.attr("physiology_time_stamp")[i] = head.physiology_time_stamp[i];
            }
            pyhead.attr("number_of_samples") = head.number_of_samples;
            pyhead.attr("available_channels") = head.available_channels;
            pyhead.attr("active_channels") = head.active_channels;

            for (int i = 0; i < ISMRMRD::ISMRMRD_CHANNEL_MASKS; i++) {
                pyhead.attr("channel_mask")[i] = head.channel_mask[i];
            }

            pyhead.attr("discard_pre") = head.discard_pre;
            pyhead.attr("discard_post") = head.discard_post;
            pyhead.attr("center_sample") = head.center_sample;
            pyhead.attr("encoding_space_ref") = head.encoding_space_ref;
            pyhead.attr("trajectory_dimensions") = head.trajectory_dimensions;
            pyhead.attr("sample_time_us") = head.sample_time_us;

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

            pyhead.attr("idx").attr("kspace_encode_step_1") = head.idx.kspace_encode_step_1;
            pyhead.attr("idx").attr("kspace_encode_step_2") = head.idx.kspace_encode_step_2;
            pyhead.attr("idx").attr("average") = head.idx.average;
            pyhead.attr("idx").attr("slice") = head.idx.slice;
            pyhead.attr("idx").attr("contrast") = head.idx.contrast;
            pyhead.attr("idx").attr("phase") = head.idx.phase;
            pyhead.attr("idx").attr("repetition") = head.idx.repetition;
            pyhead.attr("idx").attr("set") = head.idx.set;
            pyhead.attr("idx").attr("segment") = head.idx.segment;
            for (int i = 0; i < ISMRMRD::ISMRMRD_USER_INTS; i++) {
                pyhead.attr("idx").attr("user")[i] = head.idx.user[i];
            }

            for (int i = 0; i < ISMRMRD::ISMRMRD_USER_INTS; i++) {
                pyhead.attr("user_int")[i] = head.user_int[i];
            }
            for (int i = 0; i < ISMRMRD::ISMRMRD_USER_FLOATS; i++) {
                pyhead.attr("user_float")[i] = head.user_float[i];
            }

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

            head->version = bp::extract<uint16_t>(pyhead.attr("version"));
            head->flags = bp::extract<uint64_t>(pyhead.attr("flags"));
            head->measurement_uid = bp::extract<uint32_t>(pyhead.attr("measurement_uid"));
            head->scan_counter = bp::extract<uint32_t>(pyhead.attr("scan_counter"));
            head->acquisition_time_stamp = bp::extract<uint32_t>(pyhead.attr("acquisition_time_stamp"));
            for (int i = 0; i < ISMRMRD::ISMRMRD_PHYS_STAMPS; i++) {
                head->physiology_time_stamp[i] = bp::extract<uint32_t>(pyhead.attr("physiology_time_stamp")[i]);
            }
            head->number_of_samples = bp::extract<uint16_t>(pyhead.attr("number_of_samples"));
            head->available_channels = bp::extract<uint16_t>(pyhead.attr("available_channels"));
            head->active_channels = bp::extract<uint16_t>(pyhead.attr("active_channels"));
            for (int i = 0; i < ISMRMRD::ISMRMRD_CHANNEL_MASKS; i++) {
                head->channel_mask[i] = bp::extract<uint64_t>(pyhead.attr("channel_mask")[i]);
            }
            head->discard_pre = bp::extract<uint16_t>(pyhead.attr("discard_pre"));
            head->discard_post = bp::extract<uint16_t>(pyhead.attr("discard_post"));
            head->center_sample = bp::extract<uint16_t>(pyhead.attr("center_sample"));
            head->encoding_space_ref = bp::extract<uint16_t>(pyhead.attr("encoding_space_ref"));
            head->trajectory_dimensions = bp::extract<uint16_t>(pyhead.attr("trajectory_dimensions"));
            head->sample_time_us = bp::extract<float>(pyhead.attr("sample_time_us"));
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

            head->idx.kspace_encode_step_1 = bp::extract<uint16_t>(pyhead.attr("idx").attr("kspace_encode_step_1"));
            head->idx.kspace_encode_step_2 = bp::extract<uint16_t>(pyhead.attr("idx").attr("kspace_encode_step_2"));
            head->idx.average = bp::extract<uint16_t>(pyhead.attr("idx").attr("average"));
            head->idx.slice = bp::extract<uint16_t>(pyhead.attr("idx").attr("slice"));
            head->idx.contrast = bp::extract<uint16_t>(pyhead.attr("idx").attr("contrast"));
            head->idx.phase = bp::extract<uint16_t>(pyhead.attr("idx").attr("phase"));
            head->idx.repetition = bp::extract<uint16_t>(pyhead.attr("idx").attr("repetition"));
            head->idx.set = bp::extract<uint16_t>(pyhead.attr("idx").attr("set"));
            head->idx.segment = bp::extract<uint16_t>(pyhead.attr("idx").attr("segment"));
            for (int i = 0; i < ISMRMRD::ISMRMRD_USER_INTS; i++) {
                head->idx.user[i] = bp::extract<uint16_t>(pyhead.attr("idx").attr("user")[i]);
            }

            for (int i = 0; i < ISMRMRD::ISMRMRD_USER_INTS; i++) {
                head->user_int[i] = bp::extract<int32_t>(pyhead.attr("user_int")[i]);
            }
            for (int i = 0; i < ISMRMRD::ISMRMRD_USER_FLOATS; i++) {
                head->user_float[i] = bp::extract<float>(pyhead.attr("user_float")[i]);
            }
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


/// Partial specialization of `python_converter` for ISMRMRD::AcquisitionHeader
template<> struct python_converter<ISMRMRD::AcquisitionHeader> {
    static void create()
    {
        create_ismrmrd_AcquisitionHeader_converter();
    }
};

/// Partial specialization of `python_converter` for ISMRMRD::ImageHeader
template<> struct python_converter<ISMRMRD::ImageHeader> {
    static void create()
    {
        create_ismrmrd_ImageHeader_converter();
    }
};

}

#endif /* GADGETRON_PYTHON_ISMRMRD_CONVERTER_H */
