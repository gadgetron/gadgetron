/**
\file   python_IsmrmrdImageArray_converter.h
\brief  Implementation converter from/to python for IsmrmrdImageArray data structure
\author Hui Xue
*/

#pragma once
#include "python_toolbox.h"
#include "python_numpy_wrappers.h"

#include "hoNDArray.h"
#include "mri_core_data.h"
#include "log.h"

#include <boost/python.hpp>
namespace bp = boost::python;

namespace Gadgetron {

    class IsmrmrdImageArray_to_python_object
    {
    public:
        static PyObject* convert(const IsmrmrdImageArray & arrayData)
        {
            bp::object pygadgetron = bp::import("gadgetron");

            auto data = bp::object(arrayData.data_);
            auto pyHeaders = bp::list();

            size_t n;
            for (n=0; n<arrayData.headers_.get_number_of_elements(); n++)
            {
                auto headers = boost::python::object(arrayData.headers_(n));
                pyHeaders.append(headers);
            }

            auto pyMeta = bp::list();
            for (n = 0; n<arrayData.meta_.size(); n++)
            {
                auto meta = boost::python::object(arrayData.meta_[n]);
                pyMeta.append(meta);
            }

            auto buffer = pygadgetron.attr("IsmrmrdImageArray")(data, pyHeaders, pyMeta);

            // increment the reference count so it exists after `return`
            return bp::incref(buffer.ptr());
        }
    };

    /// Used for making an hoNDArray from a NumPy array
    struct IsmrmrdImageArray_from_python_object
    {
        IsmrmrdImageArray_from_python_object()
        {
            // actually register this converter with Boost
            bp::converter::registry::push_back(
                &convertible,
                &construct,
                bp::type_id<IsmrmrdImageArray>());
        }

        /// Returns NULL if the object is not convertible. Or well.... it should
        static void* convertible(PyObject* obj)
        {
            return obj;
        }

        /// Construct an hoNDArray in-place
        static void construct(PyObject* obj, bp::converter::rvalue_from_python_stage1_data* data)
        {

            void* storage = ((bp::converter::rvalue_from_python_storage<IsmrmrdImageArray >*)data)->storage.bytes;
            IsmrmrdImageArray* reconData = new (storage) IsmrmrdImageArray;
            data->convertible = storage;


            try {
                bp::list pyRecondata((bp::handle<>(bp::borrowed(obj))));
                auto length = bp::len(pyRecondata);
                GDEBUG("Recon data length: %i\n", length);
                for (int i = 0; i < length; i++) {
                    bp::object reconBit = pyRecondata[i];
                    IsmrmrdReconBit rBit;
                    rBit.data_ = extractDataBuffered(reconBit.attr("data"));
                    if (PyObject_HasAttrString(reconBit.ptr(), "ref")) {
                        rBit.ref_ = extractDataBuffered(reconBit.attr("ref"));
                    }
                    reconData->rbit_.push_back(rBit);
                }

            }
            catch (const bp::error_already_set&) {
                std::string err = pyerr_to_string();
                GERROR(err.c_str());
                throw std::runtime_error(err);
            }
        }
        static IsmrmrdDataBuffered extractDataBuffered(bp::object pyDataBuffered) {
            IsmrmrdDataBuffered result;

            result.data_ = bp::extract<hoNDArray<std::complex<float>>>(pyDataBuffered.attr("data"));
            if (PyObject_HasAttrString(pyDataBuffered.ptr(), "trajectory"))
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
            for (int i = 0; i < 3; i++) {
                auto pySL = pySampling.attr("sampling_limits")[i];
                sampling.sampling_limits_[i].min_ = bp::extract<uint16_t>(pySL.attr("min"));
                sampling.sampling_limits_[i].center_ = bp::extract<uint16_t>(pySL.attr("center"));
                sampling.sampling_limits_[i].max_ = bp::extract<uint16_t>(pySL.attr("max"));
            }
            return result;
        }


    };



    /// Partial specialization of `python_converter` for hoNDArray
    template<> struct python_converter<IsmrmrdImageArray> {
        static void create()
        {
            // register hoNDArray converter
            bp::type_info info = bp::type_id<IsmrmrdImageArray >();
            const bp::converter::registration* reg = bp::converter::registry::query(info);
            // only register if not already registered!
            if (nullptr == reg || nullptr == (*reg).m_to_python) {
                bp::to_python_converter<IsmrmrdImageArray, IsmrmrdImageArray_to_python_object >();
                IsmrmrdImageArray_from_python_object();
            }

        }
    };

}
