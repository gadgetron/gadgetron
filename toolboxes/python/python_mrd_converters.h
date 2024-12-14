#pragma once

#include "python_converters.h"
#include "python_toolbox.h"

#include "mrd/types.h"

namespace Gadgetron {
namespace Python {

struct EncodingCounters_converter {
    static void* convertible(PyObject* obj) {
        return obj;
    }

    /// Construct an mrd::EncodingCounters in-place
    static void construct(PyObject* obj, bp::converter::rvalue_from_python_stage1_data* data) {
        void* storage = ((bp::converter::rvalue_from_python_storage<mrd::EncodingCounters>*)data)->storage.bytes;

        // Placement-new of mrd::EncodingCounters in memory provided by Boost
        auto head = new (storage) mrd::EncodingCounters;
        data->convertible = storage;

        try {
            bp::object pycounters((bp::handle<>(bp::borrowed(obj))));

            head->kspace_encode_step_1 = bp::extract<std::optional<uint32_t>>(pycounters.attr("kspace_encode_step_1"));
            head->kspace_encode_step_2 = bp::extract<std::optional<uint32_t>>(pycounters.attr("kspace_encode_step_2"));
            head->average = bp::extract<std::optional<uint32_t>>(pycounters.attr("average"));
            head->slice = bp::extract<std::optional<uint32_t>>(pycounters.attr("slice"));
            head->contrast = bp::extract<std::optional<uint32_t>>(pycounters.attr("contrast"));
            head->phase = bp::extract<std::optional<uint32_t>>(pycounters.attr("phase"));
            head->repetition = bp::extract<std::optional<uint32_t>>(pycounters.attr("repetition"));
            head->set = bp::extract<std::optional<uint32_t>>(pycounters.attr("set"));
            head->segment = bp::extract<std::optional<uint32_t>>(pycounters.attr("segment"));
            head->user = bp::extract<std::vector<uint32_t>>(pycounters.attr("user"));
        } catch (const bp::error_already_set&) {
            std::string err = pyerr_to_string();
            GERROR(err.c_str());
            throw std::runtime_error(err);
        }
    }

    static PyObject* convert(const mrd::EncodingCounters& e) {
        try {
            GILLock lock;
            bp::object module = bp::import("mrd");
            bp::object cls = module.attr("EncodingCounters");
            bp::object pycounters = cls();

            pycounters.attr("kspace_encode_step_1") = e.kspace_encode_step_1;
            pycounters.attr("kspace_encode_step_2") = e.kspace_encode_step_2;
            pycounters.attr("average") = e.average;
            pycounters.attr("slice") = e.slice;
            pycounters.attr("contrast") = e.contrast;
            pycounters.attr("phase") = e.phase;
            pycounters.attr("repetition") = e.repetition;
            pycounters.attr("set") = e.set;
            pycounters.attr("segment") = e.segment;
            pycounters.attr("user") = e.user;

            // increment the reference count so it exists after `return`
            return bp::incref(pycounters.ptr());
        } catch (const bp::error_already_set&) {
            std::string err = pyerr_to_string();
            GERROR(err.c_str());
            throw std::runtime_error(err);
        }
    }
};

template <typename T, size_t... Dims>
struct FixedNDArray_converter {
    static void* convertible(PyObject* obj) {
        if (sizeof(T) != NumPyArray_ITEMSIZE(obj)) {
            GERROR("sizeof(T): %d, ITEMSIZE: %d\n", sizeof(T), NumPyArray_ITEMSIZE(obj));
            return NULL;
        }
        return obj;
    }

    static void construct(PyObject* obj, bp::converter::rvalue_from_python_stage1_data* data) {
        void* storage = ((bp::converter::rvalue_from_python_storage<yardl::FixedNDArray<T, Dims...> >*)data)->storage.bytes;
        data->convertible = storage;

        size_t length = NumPyArray_SIZE(obj);
        auto arr = new (storage) yardl::FixedNDArray<T, Dims...>();
        // memcpy(&(*vec)[0], NumPyArray_DATA(obj), sizeof(T) * arr->size());
        memcpy(arr->data(), NumPyArray_DATA(obj), sizeof(T) * arr->size());
    }

    static PyObject* convert(const yardl::FixedNDArray<T, Dims...>& arr) {
        try {
            std::vector<npy_intp> reverse_dims{Dims...};
            std::vector<npy_intp> dims(reverse_dims.rbegin(), reverse_dims.rend());

            PyObject* obj = NumPyArray_SimpleNew(dims.size(), &dims[0], get_numpy_type<T>());
            if (sizeof(T) != NumPyArray_ITEMSIZE(obj)) {
                GERROR("sizeof(T): %d, ITEMSIZE: %d\n", sizeof(T), NumPyArray_ITEMSIZE(obj));
                throw std::runtime_error("FixedNDArray_to_Python: python object and array data type sizes do not match");
            }

            memcpy(NumPyArray_DATA(obj), arr.data(), arr.size() * sizeof(T));

            // increment the reference count so it exists after `return`
            return bp::incref(obj);
        } catch (const bp::error_already_set&) {
            std::string err = pyerr_to_string();
            GERROR(err.c_str());
            throw std::runtime_error(err);
        }
    }
};

struct AcquisitionHeader_converter {
    static void* convertible(PyObject* obj) {
        return obj;
    }

    /// Construct an mrd::AcquisitionHeader in-place
    static void construct(PyObject* obj, bp::converter::rvalue_from_python_stage1_data* data) {
        void* storage = ((bp::converter::rvalue_from_python_storage<mrd::AcquisitionHeader>*)data)->storage.bytes;

        // Placement-new of mrd::AcquisitionHeader in memory provided by Boost
        auto head = new (storage) mrd::AcquisitionHeader;
        data->convertible = storage;

        try {
            bp::object pyhead((bp::handle<>(bp::borrowed(obj))));

            head->flags = bp::extract<uint64_t>(pyhead.attr("flags"))();
            head->idx = bp::extract<mrd::EncodingCounters>(pyhead.attr("idx"));
            head->measurement_uid = bp::extract<uint32_t>(pyhead.attr("measurement_uid"));
            head->scan_counter = bp::extract<std::optional<uint32_t>>(pyhead.attr("scan_counter"));
            head->acquisition_time_stamp = bp::extract<std::optional<uint32_t>>(pyhead.attr("acquisition_time_stamp"));
            head->physiology_time_stamp = bp::extract<std::vector<uint32_t>>(pyhead.attr("physiology_time_stamp"));
            head->channel_order = bp::extract<std::vector<uint32_t>>(pyhead.attr("channel_order"));
            head->discard_pre = bp::extract<std::optional<uint32_t>>(pyhead.attr("discard_pre"));
            head->discard_post = bp::extract<std::optional<uint32_t>>(pyhead.attr("discard_post"));
            head->center_sample = bp::extract<std::optional<uint32_t>>(pyhead.attr("center_sample"));
            head->encoding_space_ref = bp::extract<std::optional<uint32_t>>(pyhead.attr("encoding_space_ref"));
            head->sample_time_us = bp::extract<std::optional<float>>(pyhead.attr("sample_time_us"));
            head->position = bp::extract<yardl::FixedNDArray<float, 3>>(pyhead.attr("position"));
            head->read_dir = bp::extract<yardl::FixedNDArray<float, 3>>(pyhead.attr("read_dir"));
            head->phase_dir = bp::extract<yardl::FixedNDArray<float, 3>>(pyhead.attr("phase_dir"));
            head->slice_dir = bp::extract<yardl::FixedNDArray<float, 3>>(pyhead.attr("slice_dir"));
            head->patient_table_position = bp::extract<yardl::FixedNDArray<float, 3>>(pyhead.attr("patient_table_position"));
            head->user_int = bp::extract<std::vector<int32_t>>(pyhead.attr("user_int"));
            head->user_float = bp::extract<std::vector<float>>(pyhead.attr("user_float"));
        } catch (const bp::error_already_set&) {
            std::string err = pyerr_to_string();
            GERROR(err.c_str());
            throw std::runtime_error(err);
        }
    }

    static PyObject* convert(const mrd::AcquisitionHeader& head) {
        try {
            GILLock lock;
            bp::object module = bp::import("mrd");
            bp::object cls = module.attr("AcquisitionHeader");
            bp::object pyhead = cls();

            pyhead.attr("flags") = head.flags.Value();
            pyhead.attr("idx") = head.idx;
            pyhead.attr("measurement_uid") = head.measurement_uid;
            pyhead.attr("scan_counter") = head.scan_counter;
            pyhead.attr("acquisition_time_stamp") = head.acquisition_time_stamp;
            pyhead.attr("physiology_time_stamp") = head.physiology_time_stamp;
            pyhead.attr("channel_order") = head.channel_order;
            pyhead.attr("discard_pre") = head.discard_pre;
            pyhead.attr("discard_post") = head.discard_post;
            pyhead.attr("center_sample") = head.center_sample;
            pyhead.attr("encoding_space_ref") = head.encoding_space_ref;
            pyhead.attr("sample_time_us") = head.sample_time_us;
            pyhead.attr("position") = head.position;
            pyhead.attr("read_dir") = head.read_dir;
            pyhead.attr("phase_dir") = head.phase_dir;
            pyhead.attr("slice_dir") = head.slice_dir;
            pyhead.attr("patient_table_position") = head.patient_table_position;
            pyhead.attr("user_int") = head.user_int;
            pyhead.attr("user_float") = head.user_float;

            // increment the reference count so it exists after `return`
            return bp::incref(pyhead.ptr());
        } catch (const bp::error_already_set&) {
            std::string err = pyerr_to_string();
            GERROR(err.c_str());
            throw std::runtime_error(err);
        }
    }
};

struct ImageHeader_converter {
    static void* convertible(PyObject* obj) {
        return obj;
    }

    /// Construct an mrd::ImageHeader in-place
    static void construct(PyObject* obj, bp::converter::rvalue_from_python_stage1_data* data) {
        void* storage = ((bp::converter::rvalue_from_python_storage<mrd::ImageHeader>*)data)->storage.bytes;

        // Placement-new of mrd::ImageHeader in memory provided by Boost
        auto head = new (storage) mrd::ImageHeader;
        data->convertible = storage;

        try {
            bp::object pyhead((bp::handle<>(bp::borrowed(obj))));

            head->flags = bp::extract<uint64_t>(pyhead.attr("flags"))();
            head->measurement_uid = bp::extract<uint32_t>(pyhead.attr("measurement_uid"));
            head->field_of_view = bp::extract<yardl::FixedNDArray<float, 3>>(pyhead.attr("field_of_view"));
            head->position = bp::extract<yardl::FixedNDArray<float, 3>>(pyhead.attr("position"));
            head->col_dir = bp::extract<yardl::FixedNDArray<float, 3>>(pyhead.attr("col_dir"));
            head->line_dir = bp::extract<yardl::FixedNDArray<float, 3>>(pyhead.attr("line_dir"));
            head->slice_dir = bp::extract<yardl::FixedNDArray<float, 3>>(pyhead.attr("slice_dir"));
            head->patient_table_position = bp::extract<yardl::FixedNDArray<float, 3>>(pyhead.attr("patient_table_position"));
            head->average = bp::extract<std::optional<uint32_t>>(pyhead.attr("average"));
            head->slice = bp::extract<std::optional<uint32_t>>(pyhead.attr("slice"));
            head->contrast = bp::extract<std::optional<uint32_t>>(pyhead.attr("contrast"));
            head->phase = bp::extract<std::optional<uint32_t>>(pyhead.attr("phase"));
            head->repetition = bp::extract<std::optional<uint32_t>>(pyhead.attr("repetition"));
            head->set = bp::extract<std::optional<uint32_t>>(pyhead.attr("set"));
            head->acquisition_time_stamp = bp::extract<std::optional<uint32_t>>(pyhead.attr("acquisition_time_stamp"));
            head->physiology_time_stamp = bp::extract<std::vector<uint32_t>>(pyhead.attr("physiology_time_stamp"));
            head->image_type = static_cast<mrd::ImageType>(bp::extract<int32_t>(pyhead.attr("image_type"))());
            head->image_index = bp::extract<std::optional<uint32_t>>(pyhead.attr("image_index"));
            head->image_series_index = bp::extract<std::optional<uint32_t>>(pyhead.attr("image_series_index"));
            head->user_int = bp::extract<std::vector<int32_t>>(pyhead.attr("user_int"));
            head->user_float = bp::extract<std::vector<float>>(pyhead.attr("user_float"));
        } catch (const bp::error_already_set&) {
            std::string err = pyerr_to_string();
            GERROR(err.c_str());
            throw std::runtime_error(err);
        }
    }

    static PyObject* convert(const mrd::ImageHeader& head) {
        try {
            GILLock lock;
            bp::object module = bp::import("mrd");

            bp::tuple args;
            bp::dict kwargs;
            kwargs["image_type"] = static_cast<int32_t>(head.image_type);

            bp::object pyhead = module.attr("ImageHeader")(*args, **kwargs);

            pyhead.attr("flags") = head.flags.Value();
            pyhead.attr("measurement_uid") = head.measurement_uid;
            pyhead.attr("field_of_view") = head.field_of_view;
            pyhead.attr("position") = head.position;
            pyhead.attr("col_dir") = head.col_dir;
            pyhead.attr("line_dir") = head.line_dir;
            pyhead.attr("slice_dir") = head.slice_dir;
            pyhead.attr("patient_table_position") = head.patient_table_position;
            pyhead.attr("average") = head.average;
            pyhead.attr("slice") = head.slice;
            pyhead.attr("contrast") = head.contrast;
            pyhead.attr("phase") = head.phase;
            pyhead.attr("repetition") = head.repetition;
            pyhead.attr("set") = head.set;
            pyhead.attr("acquisition_time_stamp") = head.acquisition_time_stamp;
            pyhead.attr("physiology_time_stamp") = head.physiology_time_stamp;
            pyhead.attr("image_index") = head.image_index;
            pyhead.attr("image_series_index") = head.image_series_index;
            pyhead.attr("user_int") = head.user_int;
            pyhead.attr("user_float") = head.user_float;

            // increment the reference count so it exists after `return`
            return bp::incref(pyhead.ptr());
        } catch (const bp::error_already_set&) {
            std::string err = pyerr_to_string();
            GERROR(err.c_str());
            throw std::runtime_error(err);
        }
    }
};

struct ImageMeta_converter
{
    static void* convertible(PyObject* obj)
    {
        return obj;
    }

    /// Construct an mrd::ImageMeta in-place
    static void construct(PyObject* obj, bp::converter::rvalue_from_python_stage1_data* data)
    {
        void* storage = ((bp::converter::rvalue_from_python_storage<mrd::ImageMeta>*)data)->storage.bytes;

        mrd::ImageMeta* meta = new (storage) mrd::ImageMeta;
        data->convertible = storage;

        try {
            bp::dict pyMeta((bp::handle<>(bp::borrowed(obj))));

            bp::list pyKeys(pyMeta.keys());
            for (size_t i = 0; i < bp::len(pyKeys); i++) {
                auto pyKey = pyKeys[i];
                std::string key = bp::extract<std::string>(pyKey);
                bp::list values = bp::extract<bp::list>(pyMeta[key]);

                std::vector<mrd::ImageMetaValue> vals;
                for (int i = 0; i < bp::len(values); i++) {
                    bp::object val = values[i].attr("value");
                    if (bp::extract<std::string>(val).check()) {
                        vals.push_back(bp::extract<std::string>(val));
                    } else if (bp::extract<int64_t>(val).check()) {
                        vals.push_back(bp::extract<int64_t>(val));
                    } else if (bp::extract<double>(val).check()) {
                        vals.push_back(bp::extract<double>(val));
                    } else {
                        GERROR("ImageMeta_from_Python: unknown variant type\n");
                        throw std::runtime_error("ImageMeta_from_Python: unknown variant type");
                    }
                }

                (*meta)[key] = vals;
            }
        }
        catch (const bp::error_already_set&)
        {
            std::string err = pyerr_to_string();
            GERROR(err.c_str());
            throw std::runtime_error(err);
        }
    }

    static PyObject* convert(const mrd::ImageMeta& meta)
    {
        try
        {
            bp::object module = bp::import("mrd");
            bp::object pymeta = module.attr("ImageMeta")();

            for (const auto& pair : meta) {
                const std::string& key = pair.first;
                const std::vector<mrd::ImageMetaValue>& values = pair.second;

                bp::list pyvalues;
                for (const auto& value : values) {
                    if (std::holds_alternative<std::string>(value)) {
                        auto str = std::get<std::string>(value);
                        pyvalues.append(module.attr("ImageMetaValue").attr("String")(str));
                    } else if (std::holds_alternative<int64_t>(value)) {
                        auto l = std::get<int64_t>(value);
                        pyvalues.append(module.attr("ImageMetaValue").attr("Int64")(l));
                    } else if (std::holds_alternative<double>(value)) {
                        auto d = std::get<double>(value);
                        pyvalues.append(module.attr("ImageMetaValue").attr("Float64")(d));
                    } else {
                        GERROR("ImageMeta_to_Python: unknown variant type\n");
                        throw std::runtime_error("ImageMeta_to_Python: unknown variant type");
                    }
                }

                pymeta[key] = pyvalues;
            }

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

struct Waveform_converter
{
    static void* convertible(PyObject* obj)
    {
        return obj;
    }

    /// Construct an mrd::Waveform in-place
    static void construct(PyObject* obj, bp::converter::rvalue_from_python_stage1_data* data)
    {
        void* storage = ((bp::converter::rvalue_from_python_storage<mrd::WaveformUint32>*)data)->storage.bytes;

        mrd::WaveformUint32* wav = new (storage) mrd::WaveformUint32;
        data->convertible = storage;

        try
        {
            bp::object pyWav((bp::handle<>(bp::borrowed(obj))));
            wav->flags = bp::extract<uint64_t>(pyWav.attr("flags"));
            wav->measurement_uid = bp::extract<uint32_t>(pyWav.attr("measurement_uid"));
            wav->scan_counter = bp::extract<uint32_t>(pyWav.attr("scan_counter"));
            wav->time_stamp = bp::extract<uint32_t>(pyWav.attr("time_stamp"));
            wav->sample_time_us = bp::extract<float>(pyWav.attr("sample_time_us"));
            wav->waveform_id = bp::extract<uint32_t>(pyWav.attr("waveform_id"));

            auto data = bp::extract<hoNDArray<uint32_t>>(pyWav.attr("data"));
            wav->data = data();
        }
        catch (const bp::error_already_set&)
        {
            std::string err = pyerr_to_string();
            GERROR(err.c_str());
            throw std::runtime_error(err);
        }
    }

    static PyObject* convert(const mrd::WaveformUint32& wav)
    {
        try
        {
            // Ensure converter is registered for Waveform `data`
            register_converter<hoNDArray<uint32_t>>();

            bp::object module = bp::import("mrd");

            const hoNDArray<uint32_t>& wavdata = wav.data;
            auto data = boost::python::object(wavdata);
            bp::incref(data.ptr());

            bp::tuple args;
            bp::dict kwargs;
            kwargs["data"] = data;

            bp::object pyWav = module.attr("WaveformUint32")(*args, **kwargs);

            pyWav.attr("flags") = wav.flags;
            pyWav.attr("measurement_uid") = wav.measurement_uid;
            pyWav.attr("scan_counter") = wav.scan_counter;
            pyWav.attr("time_stamp") = wav.time_stamp;
            pyWav.attr("sample_time_us") = wav.sample_time_us;
            pyWav.attr("waveform_id") = wav.waveform_id;

            // increment the reference count so it exists after `return`
            return bp::incref(pyWav.ptr());
        }
        catch (const bp::error_already_set&)
        {
            std::string err = pyerr_to_string();
            GERROR(err.c_str());
            throw std::runtime_error(err);
        }
    }
};

class ReconData_converter
{
public:
  /// Returns NULL if the object is not convertible. Or well.... it should
  static void* convertible(PyObject* obj) {
    return obj;
  }

  /// Construct an hoNDArray in-place
  static void construct(PyObject* obj, bp::converter::rvalue_from_python_stage1_data* data) {
    void* storage = ((bp::converter::rvalue_from_python_storage<mrd::ReconData >*)data)->storage.bytes;
    mrd::ReconData* reconData = new (storage) mrd::ReconData;
    data->convertible = storage;

    try {
      bp::object pyRecondata((bp::handle<>(bp::borrowed(obj))));

      bp::list pyBuffers = bp::list(pyRecondata.attr("buffers"));

      auto length = bp::len(pyBuffers);
      for (int i = 0; i < length; i++){
        bp::object pyReconAssembly = pyBuffers[i];

        mrd::ReconAssembly assembly;
        assembly.data = extractReconBuffer(pyReconAssembly.attr("data"));

        auto pyRef = bp::object(pyReconAssembly.attr("ref"));
        if (pyRef.is_none()){
          assembly.ref = std::nullopt;
        } else {
          assembly.ref = extractReconBuffer(pyReconAssembly.attr("ref"));
        }

        reconData->buffers.push_back(assembly);
      }

    }catch (const bp::error_already_set&) {
      std::string err = pyerr_to_string();
      GERROR(err.c_str());
      throw std::runtime_error(err);
    }
  }

  static PyObject* convert(const mrd::ReconData& reconData) {
      GILLock lock;
      bp::object module = bp::import("mrd");

      bp::list pyBuffers;
      for (auto& reconAssembly : reconData.buffers) {
          auto data = ReconBufferToPython(reconAssembly.data);
          auto ref = reconAssembly.ref ? ReconBufferToPython(*reconAssembly.ref) : bp::object();

          auto pyReconAssembly = module.attr("ReconAssembly")();
          pyReconAssembly.attr("data") = data;
          pyReconAssembly.attr("ref") = ref;

          pyBuffers.append(pyReconAssembly);
      }
      bp::incref(pyBuffers.ptr());

      auto pyReconData = module.attr("ReconData")();
      pyReconData.attr("buffers") = pyBuffers;

      // increment the reference count so it exists after `return`
      return bp::incref(pyReconData.ptr());
  }

private:
  static bp::object ReconBufferToPython( const mrd::ReconBuffer & dataBuffer){
    const hoNDArray<std::complex<float>>& data = dataBuffer.data;
    auto pyData = bp::object(data);
    bp::incref(pyData.ptr());

    const hoNDArray<mrd::AcquisitionHeader>& headers = dataBuffer.headers;
    auto pyHeaders = boost::python::object(headers);
    bp::incref(pyHeaders.ptr());

    const hoNDArray<float>& trajectory = dataBuffer.trajectory;
    auto pyTraj = bp::object(trajectory);
    bp::incref(pyTraj.ptr());

    bp::object pyDensity;
    if (dataBuffer.density) {
      const hoNDArray<float>& density = *dataBuffer.density;
      pyDensity = bp::object(density);
      bp::incref(pyDensity.ptr());
    }

    auto pySampling = SamplingDescriptionToPython(dataBuffer.sampling);
    bp::incref(pySampling.ptr());

    bp::object module = bp::import("mrd");
    auto buffer = module.attr("ReconBuffer")();

    buffer.attr("data") = pyData;
    buffer.attr("trajectory") = pyTraj;
    buffer.attr("density") = pyDensity;
    buffer.attr("headers") = pyHeaders;
    buffer.attr("sampling") = pySampling;

    return buffer;
  }

  static bp::object SamplingDescriptionToPython(const mrd::SamplingDescription & sD){
    try {
        bp::object module = bp::import("mrd");
        bp::object pySd = module.attr("SamplingDescription")();

        auto encoded_fov = module.attr("FieldOfViewMm")();
        encoded_fov.attr("x") = sD.encoded_fov.x;
        encoded_fov.attr("y") = sD.encoded_fov.y;
        encoded_fov.attr("z") = sD.encoded_fov.z;
        pySd.attr("encoded_fov") = encoded_fov;

        auto encoded_matrix = module.attr("MatrixSizeType")();
        encoded_matrix.attr("x") = sD.encoded_matrix.x;
        encoded_matrix.attr("y") = sD.encoded_matrix.y;
        encoded_matrix.attr("z") = sD.encoded_matrix.z;
        pySd.attr("encoded_matrix") = encoded_matrix;

        auto recon_fov = module.attr("FieldOfViewMm")();
        recon_fov.attr("x") = sD.recon_fov.x;
        recon_fov.attr("y") = sD.recon_fov.y;
        recon_fov.attr("z") = sD.recon_fov.z;
        pySd.attr("recon_fov") = recon_fov;

        auto recon_matrix = module.attr("MatrixSizeType")();
        recon_matrix.attr("x") = sD.recon_matrix.x;
        recon_matrix.attr("y") = sD.recon_matrix.y;
        recon_matrix.attr("z") = sD.recon_matrix.z;
        pySd.attr("recon_matrix") = recon_matrix;

        auto sampling_limit_e0 = module.attr("LimitType")();
        sampling_limit_e0.attr("minimum") = sD.sampling_limits.kspace_encoding_step_0.minimum;
        sampling_limit_e0.attr("maximum") = sD.sampling_limits.kspace_encoding_step_0.maximum;
        sampling_limit_e0.attr("center") = sD.sampling_limits.kspace_encoding_step_0.center;

        auto sampling_limit_e1 = module.attr("LimitType")();
        sampling_limit_e1.attr("minimum") = sD.sampling_limits.kspace_encoding_step_1.minimum;
        sampling_limit_e1.attr("maximum") = sD.sampling_limits.kspace_encoding_step_1.maximum;
        sampling_limit_e1.attr("center") = sD.sampling_limits.kspace_encoding_step_1.center;

        auto sampling_limit_e2 = module.attr("LimitType")();
        sampling_limit_e2.attr("minimum") = sD.sampling_limits.kspace_encoding_step_2.minimum;
        sampling_limit_e2.attr("maximum") = sD.sampling_limits.kspace_encoding_step_2.maximum;
        sampling_limit_e2.attr("center") = sD.sampling_limits.kspace_encoding_step_2.center;

        auto sampling_limits = module.attr("SamplingLimits")();
        sampling_limits.attr("kspace_encoding_step_0") = sampling_limit_e0;
        sampling_limits.attr("kspace_encoding_step_1") = sampling_limit_e1;
        sampling_limits.attr("kspace_encoding_step_2") = sampling_limit_e2;
        pySd.attr("sampling_limits") = sampling_limits;

        return pySd;
    } catch (bp::error_already_set const&) {
        std::string err = pyerr_to_string();
        GERROR(err.c_str());
        throw std::runtime_error(err);
    }
  }

  static mrd::ReconBuffer extractReconBuffer(bp::object pyReconBuffer){
    mrd::ReconBuffer result;

    result.data = bp::extract<hoNDArray<std::complex<float>>>(pyReconBuffer.attr("data"))();
    result.trajectory = bp::extract<hoNDArray<float>>(pyReconBuffer.attr("trajectory"))();
    result.density = bp::extract<std::optional<hoNDArray<float>>>(pyReconBuffer.attr("density"))();
    result.headers = bp::extract<hoNDArray<mrd::AcquisitionHeader>>(pyReconBuffer.attr("headers"))();

    auto pySampling = pyReconBuffer.attr("sampling");
    mrd::SamplingDescription sampling;
    sampling.encoded_fov.x = bp::extract<float>(pySampling.attr("encoded_fov").attr("x"));
    sampling.encoded_fov.y = bp::extract<float>(pySampling.attr("encoded_fov").attr("y"));
    sampling.encoded_fov.z = bp::extract<float>(pySampling.attr("encoded_fov").attr("z"));
    sampling.encoded_matrix.x = bp::extract<uint32_t>(pySampling.attr("encoded_matrix").attr("x"));
    sampling.encoded_matrix.y = bp::extract<uint32_t>(pySampling.attr("encoded_matrix").attr("y"));
    sampling.encoded_matrix.z = bp::extract<uint32_t>(pySampling.attr("encoded_matrix").attr("z"));
    sampling.recon_fov.x = bp::extract<float>(pySampling.attr("recon_fov").attr("x"));
    sampling.recon_fov.y = bp::extract<float>(pySampling.attr("recon_fov").attr("y"));
    sampling.recon_fov.z = bp::extract<float>(pySampling.attr("recon_fov").attr("z"));
    sampling.recon_matrix.x = bp::extract<uint32_t>(pySampling.attr("recon_matrix").attr("x"));
    sampling.recon_matrix.y = bp::extract<uint32_t>(pySampling.attr("recon_matrix").attr("y"));
    sampling.recon_matrix.z = bp::extract<uint32_t>(pySampling.attr("recon_matrix").attr("z"));

    auto pySLro = pySampling.attr("sampling_limits").attr("kspace_encoding_step_0");
    sampling.sampling_limits.kspace_encoding_step_0.minimum = bp::extract<uint32_t>(pySLro.attr("minimum"));
    sampling.sampling_limits.kspace_encoding_step_0.center = bp::extract<uint32_t>(pySLro.attr("center"));
    sampling.sampling_limits.kspace_encoding_step_0.maximum = bp::extract<uint32_t>(pySLro.attr("maximum"));
    auto pySLe1 = pySampling.attr("sampling_limits").attr("kspace_encoding_step_1");
    sampling.sampling_limits.kspace_encoding_step_1.minimum = bp::extract<uint32_t>(pySLe1.attr("minimum"));
    sampling.sampling_limits.kspace_encoding_step_1.center = bp::extract<uint32_t>(pySLe1.attr("center"));
    sampling.sampling_limits.kspace_encoding_step_1.maximum = bp::extract<uint32_t>(pySLe1.attr("maximum"));
    auto pySLe2 = pySampling.attr("sampling_limits").attr("kspace_encoding_step_2");
    sampling.sampling_limits.kspace_encoding_step_2.minimum = bp::extract<uint32_t>(pySLe2.attr("minimum"));
    sampling.sampling_limits.kspace_encoding_step_2.center = bp::extract<uint32_t>(pySLe2.attr("center"));
    sampling.sampling_limits.kspace_encoding_step_2.maximum = bp::extract<uint32_t>(pySLe2.attr("maximum"));

    result.sampling = sampling;

    return result;
  }
};

struct ImageArray_converter
{
    /// Returns NULL if the object is not convertible. Or well.... it should
    static void* convertible(PyObject* obj)
    {
        return obj;
    }

    /// Construct an hoNDArray in-place
    static void construct(PyObject* obj, bp::converter::rvalue_from_python_stage1_data* data)
    {
        void* storage = ((bp::converter::rvalue_from_python_storage<mrd::ImageArray >*)data)->storage.bytes;
        mrd::ImageArray* reconData = new (storage) mrd::ImageArray;
        data->convertible = storage;

        try
        {
            bp::object pyImageArray((bp::handle<>(bp::borrowed(obj))));

            reconData->data = bp::extract<hoNDArray<std::complex<float>>>(pyImageArray.attr("data"))();
            reconData->headers = bp::extract<hoNDArray<mrd::ImageHeader>>(pyImageArray.attr("headers"))();
            reconData->meta = bp::extract<hoNDArray<mrd::ImageMeta>>(pyImageArray.attr("meta"))();
            reconData->waveforms = bp::extract<std::vector<mrd::WaveformUint32>>(pyImageArray.attr("waveforms"));
        }
        catch (const bp::error_already_set&)
        {
            std::string err = pyerr_to_string();
            GERROR(err.c_str());
            throw std::runtime_error(err);
        }
    }

    static PyObject* convert(const mrd::ImageArray& arrayData)
    {
        GILLock lock;
        bp::object module = bp::import("mrd");

        const hoNDArray<std::complex<float>>& data = arrayData.data;
        auto pyData = bp::object(data);
        bp::incref(pyData.ptr());

        const hoNDArray<mrd::ImageHeader>& headers = arrayData.headers;
        auto pyHeaders = bp::object(headers);
        bp::incref(pyHeaders.ptr());

        const hoNDArray<mrd::ImageMeta>& meta = arrayData.meta;
        auto pyMeta = bp::object(meta);
        bp::incref(pyMeta.ptr());

        auto pyWav = bp::object(arrayData.waveforms);
        bp::incref(pyWav.ptr());

        auto imarray = module.attr("ImageArray")();
        imarray.attr("data") = pyData;
        imarray.attr("headers") = pyHeaders;
        imarray.attr("meta") = pyMeta;
        imarray.attr("waveforms") = pyWav;

        // increment the reference count so it exists after `return`
        return bp::incref(imarray.ptr());
    }
};

} // namespace Python

// -------------------------------------------------------------------------------------------------------

template<> struct python_converter<mrd::EncodingCounters> {
    static void create()
    {
        register_converter<std::optional<uint32_t>>();
        register_with_boost<mrd::EncodingCounters, Python::EncodingCounters_converter>();
    }
};

template <typename T, size_t... Dims>
struct python_converter<yardl::FixedNDArray<T, Dims...>> {
    static void create()
    {
        register_with_boost<yardl::FixedNDArray<T, Dims...>, Python::FixedNDArray_converter<T, Dims...>>();
    }
};

// -------------------------------------------------------------------------------------------------------
/// Partial specialization of `python_converter` for mrd::AcquisitionHeader
template<> struct python_converter<mrd::AcquisitionHeader> {
    static void create()
    {
        register_converter<mrd::EncodingCounters>();
        register_converter<yardl::FixedNDArray<float, 3>>();
        register_converter<std::vector<uint32_t>>();
        register_converter<std::vector<int32_t>>();
        register_converter<std::vector<float>>();
        register_converter<std::optional<float>>();

        register_with_boost<mrd::AcquisitionHeader, Python::AcquisitionHeader_converter>();
    }
};

// -------------------------------------------------------------------------------------------------------
/// Partial specialization of `python_converter` for mrd::ImageHeader
template<> struct python_converter<mrd::ImageHeader> {
    static void create()
    {
        register_converter<yardl::FixedNDArray<float, 3>>();
        register_converter<std::optional<uint32_t>>();
        register_converter<std::vector<uint32_t>>();
        register_converter<std::vector<int32_t>>();
        register_converter<std::vector<float>>();

        register_with_boost<mrd::ImageHeader, Python::ImageHeader_converter>();
    }
};

// -------------------------------------------------------------------------------------------------------
/// Partial specialization of `python_converter` for mrd::ImageMeta
template<> struct python_converter<mrd::ImageMeta> {
    static void create()
    {
        register_with_boost<mrd::ImageMeta, Python::ImageMeta_converter>();
    }
};

// -------------------------------------------------------------------------------------------------------
/// Partial specialization of `python_converter` for mrd::Waveform
template<> struct python_converter<mrd::WaveformUint32> {
    static void create()
    {
        register_with_boost<mrd::WaveformUint32, Python::Waveform_converter>();
    }
};

/// Partial specialization of `python_converter` for mrd::ReconData
template <> struct python_converter<mrd::ReconData> {
    static void create()
    {
        register_converter<hoNDArray<std::complex<float>>>();
        register_converter<hoNDArray<float>>();
        register_converter<std::optional<hoNDArray<float>>>();
        register_converter<hoNDArray<mrd::AcquisitionHeader>>();

        register_with_boost<mrd::ReconData, Python::ReconData_converter>();
    }
};

template<> struct python_converter<mrd::ImageArray>
{
    static void create()
    {
        register_converter<hoNDArray<std::complex<float>>>();
        register_converter<hoNDArray<mrd::ImageHeader>>();
        register_converter<hoNDArray<mrd::ImageMeta>>();
        register_converter<std::vector<mrd::WaveformUint32>>();

        register_with_boost<mrd::ImageArray, Python::ImageArray_converter>();
    }
};

} // namespace Gadgetron
