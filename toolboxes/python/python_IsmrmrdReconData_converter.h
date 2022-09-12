#pragma once
#include "python_toolbox.h"
#include "python_numpy_wrappers.h"

#include "hoNDArray.h"
#include "mri_core_data.h"
#include "log.h"

#include <boost/python.hpp>
namespace bp = boost::python;

namespace Gadgetron {

class IsmrmrdReconData_to_python_object {
public:
  static PyObject* convert(const IsmrmrdReconData & reconData) {
//      initialize_python();
      GILLock lock;
    bp::object pygadgetron = bp::import("gadgetron");

    auto pyReconData = bp::list();
    for (auto & reconBit : reconData.rbit_ ){
      auto data = DataBufferedToPython(reconBit.data_);
      auto ref = 	reconBit.ref_ ? DataBufferedToPython(*reconBit.ref_) : bp::object();

      auto pyReconBit = pygadgetron.attr("IsmrmrdReconBit")(data,ref);
      pyReconData.append(pyReconBit);

    }
    // increment the reference count so it exists after `return`
    return bp::incref(pyReconData.ptr());
  }

private:
  static bp::object DataBufferedToPython( const IsmrmrdDataBuffered & dataBuffer){
    bp::object pygadgetron = bp::import("gadgetron");
    auto data = bp::object(dataBuffer.data_);
    auto headers = boost::python::object(dataBuffer.headers_);
    auto trajectory = dataBuffer.trajectory_ ? bp::object(*dataBuffer.trajectory_) : bp::object();
    auto sampling = SamplingDescriptionToPython(dataBuffer.sampling_);
    auto buffer = pygadgetron.attr("IsmrmrdDataBuffered")(data,headers,sampling,trajectory);

    bp::incref(data.ptr());
    bp::incref(headers.ptr());
    bp::incref(trajectory.ptr());
    bp::incref(sampling.ptr());

    return buffer;
  }

  static bp::object SamplingDescriptionToPython(const SamplingDescription & sD){
    bp::object pygadgetron = bp::import("gadgetron");
    bp::object result;
    try {
     auto tmp = pygadgetron.attr("SamplingDescription");
     result = tmp();
    result.attr("encoded_FOV") = bp::make_tuple(sD.encoded_FOV_[0],sD.encoded_FOV_[1],sD.encoded_FOV_[2]);
    result.attr("encoded_matrix") = bp::make_tuple(sD.encoded_matrix_[0],sD.encoded_matrix_[1],sD.encoded_matrix_[2]);
    result.attr("recon_FOV") = bp::make_tuple(sD.recon_FOV_[0],sD.recon_FOV_[1],sD.recon_FOV_[2]);
    result.attr("recon_matrix") = bp::make_tuple(sD.recon_matrix_[0],sD.recon_matrix_[1],sD.recon_matrix_[2]);
 } catch (bp::error_already_set const &){
        std::string err = pyerr_to_string();
        GERROR(err.c_str());
        throw std::runtime_error(err);
    }
     std::vector<bp::object> SL;
      for (int i = 0; i < 3; i++){
        auto sampling_limit = pygadgetron.attr("SamplingLimit")();
        sampling_limit.attr("min") = sD.sampling_limits_[i].min_;
        sampling_limit.attr("max") = sD.sampling_limits_[i].max_;
        sampling_limit.attr("center") = sD.sampling_limits_[i].center_;
        SL.push_back(sampling_limit);
      }
      result.attr("sampling_limits") = bp::make_tuple(SL[0],SL[1],SL[2]);

    return result;
  }
};


/// Used for making an hoNDArray from a NumPy array
struct IsmrmrdReconData_from_python_object {
  IsmrmrdReconData_from_python_object() {
    // actually register this converter with Boost
    bp::converter::registry::push_back(
        &convertible,
        &construct,
        bp::type_id<IsmrmrdReconData>());
  }

  /// Returns NULL if the object is not convertible. Or well.... it should
  static void* convertible(PyObject* obj) {
    return obj;
  }

  /// Construct an hoNDArray in-place
  static void construct(PyObject* obj, bp::converter::rvalue_from_python_stage1_data* data) {
    void* storage = ((bp::converter::rvalue_from_python_storage<IsmrmrdReconData >*)data)->storage.bytes;
    IsmrmrdReconData* reconData = new (storage) IsmrmrdReconData;
    data->convertible = storage;


    try {
      bp::list pyRecondata((bp::handle<>(bp::borrowed(obj))));
      auto length = bp::len(pyRecondata);
      GDEBUG("Recon data length: %i\n",length);
      for (int i = 0; i < length; i++){
        bp::object reconBit = pyRecondata[i];
        IsmrmrdReconBit rBit;
        rBit.data_ = extractDataBuffered(reconBit.attr("data"));
        if (PyObject_HasAttrString(reconBit.ptr(),"ref")){
          rBit.ref_ = extractDataBuffered(reconBit.attr("ref"));
        }
        reconData->rbit_.push_back(rBit);
      }

    }catch (const bp::error_already_set&) {
      std::string err = pyerr_to_string();
      GERROR(err.c_str());
      throw std::runtime_error(err);
    }
  }
  static IsmrmrdDataBuffered extractDataBuffered(bp::object pyDataBuffered){
    IsmrmrdDataBuffered result;

    result.data_ = bp::extract<hoNDArray<std::complex<float>>>(pyDataBuffered.attr("data"));
    if (PyObject_HasAttrString(pyDataBuffered.ptr(),"trajectory"))
      result.trajectory_ = bp::extract<hoNDArray<float>>(pyDataBuffered.attr("trajectory"));

    result.headers_ = bp::extract<hoNDArray<ISMRMRD::AcquisitionHeader>>(pyDataBuffered.attr("headers"));

    auto pySampling = pyDataBuffered.attr("sampling");
    SamplingDescription sampling;
    for (int i = 0; i < 3; i++)
      sampling.encoded_FOV_[i] = bp::extract<float>(pySampling.attr("encoded_FOV")[i]);
    for (int i = 0; i < 3; i++)
      sampling.encoded_matrix_[i] = bp::extract<uint16_t>(pySampling.attr("encoded_matrix")[i]);
    for (int i = 0; i < 3; i++)
      sampling.recon_FOV_[i] = bp::extract<float>(pySampling.attr("recon_FOV")[i]);
    for (int i = 0; i < 3; i++)
      sampling.recon_matrix_[i] = bp::extract<uint16_t>(pySampling.attr("recon_matrix")[i]);
    for (int i = 0; i < 3; i++){
      auto pySL = pySampling.attr("sampling_limits")[i];
      sampling.sampling_limits_[i].min_ = bp::extract<uint16_t>(pySL.attr("min"));
      sampling.sampling_limits_[i].center_ = bp::extract<uint16_t>(pySL.attr("center"));
      sampling.sampling_limits_[i].max_ = bp::extract<uint16_t>(pySL.attr("max"));
    }
    result.sampling_ = sampling;
    return result;
  }


};


/// Partial specialization of `python_converter` for hoNDArray
template<> struct python_converter<IsmrmrdReconData> {
  static void create()
  {
    // register hoNDArray converter
    bp::type_info info = bp::type_id<IsmrmrdReconData >();
    const bp::converter::registration* reg = bp::converter::registry::query(info);
    // only register if not already registered!
    if (nullptr == reg || nullptr == (*reg).m_to_python) {
      bp::to_python_converter<IsmrmrdReconData, IsmrmrdReconData_to_python_object >();
      IsmrmrdReconData_from_python_object();
    }

  }
};

}
